###PR20220206 VitC differentiation data

## Project variables
PROJECTNAME="DamnChIC"; 
LIMSID="KIN8403";

mkdir -p /hpc/hub_kind/$USER/projects/$PROJECTNAME/$LIMSID;
cd /hpc/hub_kind/$USER/projects/$PROJECTNAME/$LIMSID;
mkdir data;
mkdir data/raw;
mkdir metadata; 
mkdir bin;
mkdir filelists;
mkdir joblog;
mkdir tmp;

## get add_read_prefix.awk and run_preformatted_cmd.sh
tutorial_file_folder="/hpc/hub_kind/koos/tmp/";
tar -xvf ${tutorial_file_folder}/DamIDseq_tutorial_scripts.tar.gz;
cp /hpc/hub_kind/prullens/projects/ESC/KIN3147/bin/run_preformatted_cmd_slurm.sh ./bin/;


## Dowload data from USEQ
curl -u 'MggEKtde6dwajLa:ga3ZnJ38f48c' -H 'X-Requested-With: XMLHttpRequest' 'https://ncie01.op.umcutrecht.nl/public.php/webdav/i46-SK-DamChIC-VitC-Dam-K27-K9-plate8.tar' -o i46-SK-DamChIC-VitC-Dam-K27-K9-plate8.tar


## generate index specific barcode files
ANNOFN="./metadata/KIN8403_anno.tsv";
DAMID_BC_FN="./metadata/damid_v3_set1.barcodes.tsv";
OUT_FMT="./metadata/%s.index%02d.barcodes.tsv";
BC_FMT="BC_%03d";

#append barcodes to index specific barcode file
tail -n+2 $ANNOFN | cut -f1,2 | while read indexnr bcnr; do
    printf -v outfmt "$OUT_FMT" "$LIMSID" $indexnr;
    printf -v bc "$BC_FMT" $bcnr;
    seq=$(grep $bc $DAMID_BC_FN | cut -f2);
    echo -e DamID2_${bc}"\t"${seq} >> $outfmt;
done


## demultiplexing !!!! Paired-End demultiplex !!!! to enable paired-end alignment if desirable 
taskname="demultiplex";
filelist="./filelists/$(date +"%Y%m%d_%H%M").${taskname}.cmds.list";

ANNOFN="./metadata/KIN8403_anno.tsv";
LIMSID=$LIMSID;
INDIR="./tmp/";
OUTDIR="./data/demultiplexed/";
BASE_OUTFMT="${OUTDIR}/%s.index%02d";
BCFN_FMT="./metadata/%s.index%02d.barcodes.tsv";

tail -n+2 $ANNOFN | cut -f1 | sort | uniq | while read indexnr; do
    printf -v bcfile "$BCFN_FMT" "$LIMSID" $indexnr;

    fn1_infiles="";
    fn2_infiles="";
    for fn1 in $(ls ${INDIR}* | grep "i${indexnr}" | grep "R1_001.fastq.gz"); do
        fn2="${fn1/%_R1_001.fastq.gz/_R2_001.fastq.gz}";
        if ! [ -r "$fn2" ]; then
            echo "Error! File not found or not readable: ${fn2}";
            break;
        fi;
        fn1_infiles="${fn1_infiles} \"${fn1}\"";
        fn2_infiles="${fn2_infiles} \"${fn2}\"";
    done;
    
    if ! [ -r "$fn2" ]; then
        echo "Error! File not found or not readable: ${fn2}";
        break;
    fi;

    printf -v outbase "$BASE_OUTFMT" "$LIMSID" $indexnr;
    outfmt="${outbase}.{name}.{readname}.fastq.gz";
    ambiguous_outfmt="${outbase}.ambiguous.{readname}.fastq.gz";
    unmatched_outfmt="${outbase}.unmatched.{readname}.fastq.gz";
    stats_outfile="${outbase}.demultiplex_info.tsv";
    errfile="${outbase}.demultiplex.err";

    cmd="\
        demultiplex \
            -vvv \
            -m 0 \
            -o \"${outfmt}\" \
            --unmatched-outfile \"${unmatched_outfmt}\" \
            --ambiguous-outfile \"${ambiguous_outfmt}\" \
            --infofile \"${stats_outfile}\" \
            \"${bcfile}\" \
            <( cat ${fn1_infiles} | gzip -dc ) <( cat ${fn2_infiles} | gzip -dc ) \
            2> \"${errfile}\" \
    ";

    echo "$cmd";
         
done > "$filelist";

mkdir -p "${OUTDIR}" || true;

sbatch -J $LIMSID.${taskname} -o './joblog/%A.task_%a.out' -e './joblog/%A.task_%a.err' --time=10:00:00 --mem=15G --array=1-8 ./bin/run_preformatted_cmd_slurm.sh "$filelist" $(cat $filelist | wc -l);
#jobid 20706282


## aligning DamID2 data to mm10
taskname="align_damid2_data";
filelist="./filelists/$(date +"%Y%m%d_%H%M").${taskname}.cmds.list";

ANNOFN="./metadata/KIN8403_anno.tsv";
LIMSID=$LIMSID;
FNFMT="./data/demultiplexed/$LIMSID.index%02d.DamID2_BC_%03d.R1.fastq.gz";
OUTFNFMT="./data/aligned/$LIMSID.index%02d.DamID2_BC_%03d.sorted.bam";
OUTDIR="./data/aligned/";
BOWTIE2_INDEX="/hpc/hub_kind/koos/references/mouse/mm10/bowtie2/Mus_musculus.GRCm38.dna.primary_assembly.with_ERCC";
PREFIX="GA";

tail -n+2 $ANNOFN | cut -f1,2 | while read indexnr bcnr; do
    printf -v infiles $FNFMT $indexnr $bcnr;
    
    printf -v outfn "$OUTFNFMT" $indexnr $bcnr;
    bowtie2_errfile="${outfn/%sorted.bam/bowtie2.err}";
    samtools_sort_tmp="${outfn/%sorted.bam/samtools_sort_tmp}";

    cmd="\
        cat $infiles \
        | gzip -dc \
        | awk -f ./bin/add_read_prefix.awk -v PREFIX=\"${PREFIX}\" \
        | bowtie2 \
            --seed 42 \
            --very-sensitive -N 1 \
            -x \"${BOWTIE2_INDEX}\" \
            -U - \
            2> \"${bowtie2_errfile}\" \
        | samtools view -ub - \
        | samtools sort -m 500M -T \"${samtools_sort_tmp}\" -l 9 - \
        > \"${outfn}\" \
        ;";

    echo "$cmd";   
done > "$filelist";    

mkdir -p "$OUTDIR" || true;

sbatch -J $LIMSID.${taskname} -o './joblog/%A.task_%a.out' -e './joblog/%A.task_%a.err' --time=10:00:00 --mem=20G --array=1-32 ./bin/run_preformatted_cmd_slurm.sh "$filelist" $(cat $filelist | wc -l);
#jobid 20710991


## counting DamID2 data to mm10
taskname="count_damid_data_mm10";
filelist="./filelists/$(date +"%Y%m%d_%H%M").${taskname}.cmds.list";

ANNOFN="./metadata/KIN8403_anno.tsv";
INDIR="./data/aligned/";
OUTDIR="./data/counts/";
INFNFMT="${INDIR}/${LIMSID}.index%02d.DamID2_BC_%03d.sorted.bam";
ERRFNFMT="${OUTDIR}/$LIMSID.index%02d.DamID2_BC_%03d.event_counts.err";
POSFILE="/hpc/hub_kind/koos/references/mouse/mm10/posarray/Mus_musculus.GRCm38.dna.primary_assembly.with_ERCC.GATC.posarray.hdf5"
MIN_MAPQ=10;
KEEPN=2;

tail -n+2 $ANNOFN | cut -f1,2 | while read indexnr bcnr; do 
    printf -v fn "$INFNFMT" $indexnr $bcnr;
    printf -v errfn "$ERRFNFMT" $indexnr $bcnr;
    bn=$(basename "$fn");
    b=${bn/%.sorted.bam/};
    
    of="${OUTDIR}/${b}.top_n_${KEEPN}.event_counts.pos.hdf5";
    of_invalid="${OUTDIR}/${b}.invalid_pos.sorted.bam";
    cmd="\
        write_counts_at_pos.py \
            -vvv \
            --umi-hamming-dist 1 \
            --keepn $KEEPN \
            --min-mapq $MIN_MAPQ \
            --unique-event-outfile \"${of}\" \
            --invalid-pos-outfile \"${of_invalid}\" \
            --posfile \"${POSFILE}\" \
            \"${fn}\" \
            2> \"${errfn}\" \
    ;";
    
    echo "$cmd";
    
done > "$filelist";

mkdir -p "$OUTDIR" || true;

sbatch -J $LIMSID.${taskname} -o './joblog/%A.task_%a.out' -e './joblog/%A.task_%a.err' --time=5:00:00 --mem=15G --array=1-16 ./bin/run_preformatted_cmd_slurm.sh "$filelist" $(cat $filelist | wc -l);
#jobid 20725109


##(Re)align InvalidPos as ChIC reads to remove the GA prefix that was added during DamID read alignment 
## to mm10
taskname="AlignDamIDInvalidPos_mm10";
filelist="./filelists/$(date +"%Y%m%d_%H%M").${taskname}.cmds.list";

ANNOFN="./metadata/KIN8403_anno.tsv";
LIMSID=$LIMSID;
FNFMT="./data/counts/${LIMSID}.index%02d.DamID2_BC_%03d.invalid_pos.sorted.bam";
OUTFNFMT="./data/aligned/$LIMSID.index%02d.BC_ChIC_%03d.chic.sorted.bam";
OUTDIR="./data/aligned/";
HISAT2_INDEX="/hpc/hub_kind/koos/references/mouse/mm10/hisat2/Mus_musculus.GRCm38_89.primary_assembly.with_ERCC";


tail -n+2 $ANNOFN | cut -f1,2 | while read indexnr bcnr; do
    printf -v infiles $FNFMT $indexnr $bcnr;
    
    printf -v outfn "$OUTFNFMT" $indexnr $bcnr;
    hisat2_errfile="${outfn/%sorted.bam/hisat2.err}";
    samtools_sort_tmp="${outfn/%sorted.bam/samtools_sort_tmp}";

    cmd="\
        hisat2 \
          --seed 42 \
          -x \"${HISAT2_INDEX}\" \
          -U <(samtools fastq ${infiles} | fastx_trimmer -f 3 -) \
          --no-spliced-alignment \
          --mp '2,0' \
          --sp '4,0' \
          2> \"${hisat2_errfile}\" \
        | samtools view -ub - \
        | samtools sort -m 500M -T \"${samtools_sort_tmp}\" -l 9 - \
        > \"${outfn}\" \
        ;";

    echo "$cmd";  
    
done > "$filelist";

sbatch -J $LIMSID.${taskname} -o './joblog/%A.task_%a.out' -e './joblog/%A.task_%a.err' --time=5:00:00 --mem=10G --array=1-32 ./bin/run_preformatted_cmd_slurm.sh "$filelist" $(cat $filelist | wc -l);
#jobid 20734019


## bin DamID2 data to mm10
taskname="bin_damid2_counts_mm10";
filelist="./filelists/$(date +"%Y%m%d_%H%M").${taskname}.cmds.list";

ANNOFN="./metadata/KIN8403_anno.tsv";
LIMSID=$LIMSID;
MAPFILE="/hpc/hub_kind/koos/references/mouse/mm10/mappability/Mus_musculus.GRCm38.dna.primary_assembly.with_ERCC.GATC.bowtie2_very_sensitive_N1.readlength_64.counts.pos.hdf5"
POSFILE="/hpc/hub_kind/koos/references/mouse/mm10/posarray/Mus_musculus.GRCm38.dna.primary_assembly.with_ERCC.GATC.posarray.hdf5"
INDIR="./data/counts/";
OUTDIR="./data/counts/";
BINSIZES="100000";
FNFMT="${LIMSID}.index%02d.DamID2_BC_%03d.top_n_%d.event_counts.pos.hdf5"
OUTFNFMT="${LIMSID}.index%02d.DamID2_BC_%03d.top_n_%d.event_counts.binsize_%d.hdf5"
KEEPN=2;

tail -n+2 $ANNOFN | cut -f1,2 | while read indexnr barcodenr; do
    for binsize in $BINSIZES; do
        printf -v outfn "${OUTDIR}/${OUTFNFMT}" $indexnr $barcodenr $KEEPN $binsize;
        printf -v fn "${INDIR}/${FNFMT}" $indexnr $barcodenr $KEEPN;
        
        cmd="\
            bin_countfile.py \
                -vvv \
                --keepn $KEEPN \
                --mapfile \"${MAPFILE}\" \
                --posfile \"${POSFILE}\" \
                --binsize ${binsize} \
                --outfile \"${outfn}\" \
                \"${fn}\" \
        ;";
        
        echo "$cmd";
        
    done
    
done > "$filelist";

sbatch -J $LIMSID.${taskname} -o './joblog/%A.task_%a.out' -e './joblog/%A.task_%a.err' --time=1:00:00 --mem=5G --array=1-16 ./bin/run_preformatted_cmd_slurm.sh "$filelist" $(cat $filelist | wc -l);
#jobid 20735179


## bin ChIC data
taskname="BinUMIUniqueGenomicSignal";
filelist="./filelists/$(date +"%Y%m%d_%H%M").${taskname}.cmds.list";

ANNOFN="./metadata/KIN8403_anno.tsv";
LIMSID=$LIMSID;
INDIR="./data/aligned/";
OUTDIR="./data/counts/";
BINSIZES="100000";
FNFMT="${LIMSID}.index%02d.BC_ChIC_%03d.chic.sorted.bam";
OUTFNFMT="${LIMSID}.index%02d.BC_ChIC_%03d.top_n_%d.chic.AT_noTC.event_counts.binsize_%d.hdf5"
KEEPN=1;

tail -n+2 $ANNOFN | cut -f1,2 | while read indexnr barcodenr; do
    for binsize in $BINSIZES; do
        printf -v outfn "${OUTDIR}/${OUTFNFMT}" $indexnr $barcodenr $KEEPN $binsize;
        printf -v fn "${INDIR}/${FNFMT}" $indexnr $barcodenr;
        
        
        cmd="\
            python ./bin/bin_genomic_signal.py \
                --max-readlength 250 \
                --leading AT \
                --leading-excl TC \
                --binsize ${binsize} \
                --keepn $KEEPN \
                --outfile \"${outfn}\" \
                \"${fn}\" \
        ;";                
               
        echo "$cmd";
        
    done
    
done > "$filelist";

sbatch -J $LIMSID.${taskname} -o './joblog/%A.task_%a.out' -e './joblog/%A.task_%a.err' --time=1:00:00 --mem=5G --array=1-16 ./bin/run_preformatted_cmd_slurm.sh "$filelist" $(cat $filelist | wc -l);
#jobid 20736235
