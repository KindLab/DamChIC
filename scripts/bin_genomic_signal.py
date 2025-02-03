"""

Based on KR's process_scchic_bamfile.py

"""


from itertools import groupby
from collections import Counter

import numpy as np
import h5py
from sortedcontainers import SortedDict
from pysam import AlignmentFile
from tqdm import tqdm


class SortNearlySorted:
    """
    Reads on the reverse strand are sorted by 3'-end (leftmost genomic
    coordinate).
    This is impractical for counting UMI-unique reads, as we want to see all
    reads grouped by the 5' ends of the molecule, and reads are not always of
    the same length...

    See also `pos_5p_end`
    """

    def __init__(self, it, key=None, mindist=2500, maxdist=100000):
        self._it = it
        self._key = (key if key is not None else (lambda x: x))

        assert maxdist >= mindist
        self.mindist = int(mindist)
        self.maxdist = int(maxdist)
        self._d = SortedDict()

    def __iter__(self):
        return self

    def __next__(self):
        if not self._d:
            self._fill()

        if self._dist < self.mindist:
            try:
                self._fill()
            except StopIteration:
                pass

        items = self._d[self._lowest]
        if len(items) == 1:
            item = next(iter(self._d.pop(self._lowest)))
            self._update_limits()
        else:
            item = items.pop()

        return item

    def _fill(self):
        while True:
            try:
                item = next(self._it)
            except StopIteration:
                # bubble up
                raise

            k = self._key(item)
            if k in self._d:
                self._d[k].add(item)
            else:
                self._d[k] = {item}

            dist = k - self._d.keys()[0]
            if dist > self.maxdist or dist < 0:
                break

        self._update_limits()
        return

    def _update_limits(self):
        if self._d:
            self._lowest = self._d.keys()[0]
            self._highest = self._d.keys()[-1]
            self._dist = self._highest - self._lowest

        return


def pos_5p_end(read):
    if read.is_reverse:
        return read.pos + read.query_alignment_length
    else:
        return read.pos


def get_umi(read):
    s = read.query_name
    fields = s.split(":")
    assert fields[-2] == "UMI"
    return fields[-1]


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--max-readlength", default=500, required=False)
    ap.add_argument("--binsize", type=int, required=True)
    ap.add_argument("--keepn", type=int, default=1, required=False)
    ap.add_argument("--leading", type=str, default='ACGTN', help="First base of read must be in provided bases.", required=False)
    ap.add_argument("--leading-excl", type=str, default='XX', help="Reads starting with these two bases will be excluded.", required=False)
    ap.add_argument("--outfile", metavar="HDF5_FILE", required=True)
    ap.add_argument("infile", nargs=1)

    args = ap.parse_args()

    is_valid_read = lambda read: (not read.is_unmapped) and (not read.is_secondary) and (read.mapq >= 10)

    bf = AlignmentFile(args.infile[0], 'rb')
    bfiter = filter(is_valid_read, bf)

    progress = tqdm()

    chromgroups = groupby(bfiter, key=lambda read: read.reference_name)
    for chrom, readiter in chromgroups:
        l = bf.get_reference_length(chrom)
        nbins = int(np.ceil(l / args.binsize))
        counts = np.zeros((nbins, 2), dtype="uint32")

        sortedit = SortNearlySorted(readiter, key=pos_5p_end)
        #readdups = groupby(readiter, key=lambda read: (read.is_reverse, pos_5p_end(read), get_umi(read)))
        posgroups = groupby(sortedit, key=lambda read: (read.is_reverse, pos_5p_end(read), read.get_forward_sequence()[:2]))
        for (is_reverse, pos, firstbases), posreaditer in posgroups:
            if firstbases[0] in args.leading and firstbases != args.leading_excl:
                umis = Counter()
                for read in posreaditer:
                    umi = get_umi(read)
                    umis[umi] += 1
                    progress.update(1)

                # TODO filter on Hamming distance
                n_events = len(umis)
                #print(chrom, is_reverse, pos, umis, n_events)
                c = pos // args.binsize
                counts[c, int(is_reverse)] += min(args.keepn, n_events)

        with h5py.File(args.outfile) as f:
            f.create_dataset(chrom, data=counts, compression=3)

    return


if __name__ == "__main__":
    main()
