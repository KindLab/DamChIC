###PR20230207
## xci in VitC allelic PE alignment

## Project variables
PROJECTNAME="DamnChIC"; 
LIMSID="KIN8403";

############################################################################################
### Processing allele-specific DamID data using PE information inspired by FR's approach ###
###/hpc/hub_kind/franka/projects/preimplantation/KIN6777/CLUSTERNOTES_F1xCAST###############
############################################################################################


## align DamID2 data - to both parental genotypes - PE
taskname="align_damid_data_all_genotypes_PE";
filelist="./filelists/$(date +"%Y%m%d_%H%M").${taskname}.cmds.list";

REFDATA="\
129S1_SvImJ_snps\t/hpc/hub_kind/koos/references/mouse/129S1_SvImJ_snps/bowtie2/129S1_SvImJ.from_GRCm38_dna_primary_assembly.mgp_v5_snps_filter_PASS.with_contigs.with_ERCC\nCAST_snps\t/hpc/hub_kind/koos/references/mouse/CAST_snps/bowtie2/CAST_EiJ.from_GRCm38_dna_primary_assembly.mgp_v5_snps_filter_PASS.with_ERCC\
"

ANNOFN="./metadata/KIN8403_anno.tsv";
OUTDIR="./data/aligned/";
FNFMT="./data/demultiplexed/${LIMSID}.index%02d.DamID2_BC_%03d.R1.fastq.gz";
OUTFNFMT="${OUTDIR}/${LIMSID}.index%02d.DamID2_BC_%03d.%s.PE.bam";
PREFIX="GA";

tail -n+2 "$ANNOFN" | cut -f1,2 | while read indexnr bcnr; do
    printf -v fn1 "$FNFMT" $indexnr $bcnr;
    fn2=$(echo $fn1 | sed -e 's/R1\.fastq/R2\.fastq/g');

    echo -e $REFDATA | cut -f1,2 | while read genotype bowtie2_index; do

        printf -v outfn "$OUTFNFMT" $indexnr $bcnr "$genotype";
        bowtie2_errfile="${outfn/%bam/bowtie2.err}";
        
        if [ ! -f $outfn ]; then

            cmd="\
                bowtie2 \
                    --seed 42 \
                    --very-sensitive -N 1 \
                    -x \"${bowtie2_index}\" \
                    -1 <(gzip -dc \"${fn1}\"  | awk -f ./bin/add_read_prefix.awk -v PREFIX=\"${PREFIX}\") \
                    -2 \"${fn2}\" \
                    2> \"${bowtie2_errfile}\" \
                | samtools view -b - \
                > \"${outfn}\" \
                ;";

            echo "$cmd";
        
        fi;
        
    done;
done > "$filelist";

sbatch -J $LIMSID.${taskname} -o './joblog/%A.task_%a.out' -e './joblog/%A.task_%a.err' --time=5:00:00 --mem=20G --array=1-16 ./bin/run_preformatted_cmd_slurm.sh "$filelist" $(cat $filelist | wc -l);
#jobid 20717083


## split PE alignments by genotype
taskname="split_genotype_DamID_PE";
filelist="./filelists/$(date +"%Y%m%d_%H%M").${taskname}.cmds.list";

ANNOFN="./metadata/KIN8403_anno.tsv";
gt1="129S1_SvImJ_snps"
gt2="CAST_snps";
FNFMT="./data/aligned/${LIMSID}.index%02d.DamID2_BC_%03d.%s.PE.bam";
OUTFNFMT="./data/aligned/${LIMSID}.index%02d.DamID2_BC_%03d.${gt1}_vs_${gt2}.{name}.PE_R1.readname_sorted.bam";

tail -n+2 "$ANNOFN" | cut -f1,2 | while read indexnr bcnr; do
    printf -v fn1 "$FNFMT" $indexnr $bcnr "$gt1";
    printf -v fn2 "$FNFMT" $indexnr $bcnr "$gt2";
    printf -v outfmt "$OUTFNFMT" $indexnr $bcnr;
    errfn="${outfmt/%\{name\}.PE_R1.readname_sorted.bam/PE_R1.split_genotype.err}";
    
    cmd="\
        python ./bin/split_PE_DamID_alignments.py \
            -vvv \
            --unsorted-input \
            --name \"${gt1}\" --name \"${gt2}\" \
            --outfmt \"${outfmt}\" \
            \"${fn1}\" \"${fn2}\" \
        2> \"${errfn}\" \
    ;";

    echo "$cmd";
done > "$filelist";

sbatch -J $LIMSID.${taskname} -o './joblog/%A.task_%a.out' -e './joblog/%A.task_%a.err' --time=1:00:00 --mem=5G --array=1-32 ./bin/run_preformatted_cmd_slurm.sh "$filelist" $(cat $filelist | wc -l);
#jobid 20779459


## position sorting split genotype reads
taskname="position_sort_split_gt_DamID_alignments_PE";
filelist="./filelists/$(date +"%Y%m%d_%H%M").${taskname}.cmds.list";

ANNOFN="./metadata/KIN8403_anno.tsv";
gt1="129S1_SvImJ_snps"
gt2="CAST_snps";
GENOTYPES="unique_${gt1} unique_${gt2} ambiguous_${gt1}_${gt2}";
INDIR="./data/aligned/";
FNFMT="${INDIR}/${LIMSID}.index%02d.DamID2_BC_%03d.${gt1}_vs_${gt2}.%s.PE_R1.readname_sorted.bam";

tail -n+2 "$ANNOFN" | cut -f1,2 | while read indexnr bcnr; do
    for gt in $GENOTYPES; do
        printf -v fn "$FNFMT" $indexnr $bcnr "$gt";
        samtools_sort_tmp="${fn/%readname_sorted.bam/sort_tmp}";
        outfn="${fn/%readname_sorted.bam/sorted.bam}";

        cmd="\
            samtools view -ub \"${fn}\" \
            | samtools sort -m 2500M -T \"${samtools_sort_tmp}\" -l 9 - \
            > \"${outfn}\" \
        ;";

        echo $cmd;
    done;
done > "$filelist";

sbatch -J $LIMSID.${taskname} -o './joblog/%A.task_%a.out' -e './joblog/%A.task_%a.err' --time=2:00:00 --mem=5G --array=1-32 ./bin/run_preformatted_cmd_slurm.sh "$filelist" $(cat $filelist | wc -l);
#jobid 20802140


#delete original alignment and readname_sorted bam files
#rm data/aligned/*Dam*snps.PE.bam
#rm data/aligned/*Dam*snps.PE_R1.readname_sorted.bam

## process split gt reads to event counts
# using a new script that can consolidate observations across genotypes
taskname="process_damid_alignments_to_counts_PE";
filelist="./filelists/$(date +"%Y%m%d_%H%M").${taskname}.cmds.list";

ANNOFN="./metadata/KIN8403_anno.tsv";
gt1="129S1_SvImJ_snps"
gt2="CAST_snps";
GENOTYPES="unique_${gt1} unique_${gt2} ambiguous_${gt1}_${gt2}";

INDIR="./data/aligned/";
OUTDIR="./data/counts/";
FNFMT="${INDIR}/${LIMSID}.index%02d.DamID2_BC_%03d.${gt1}_vs_${gt2}.%s.PE_R1.sorted.bam";
OUTFMT="${OUTDIR}/${LIMSID}.index%02d.DamID2_BC_%03d.${gt1}_vs_${gt2}.{name}.top_n_{keepn}.event_counts.PE_R1.pos.hdf5";

POSFN1="/hpc/hub_kind/koos/references/mouse/129S1_SvImJ_snps/posarray/129S1_SvImJ.from_GRCm38_dna_primary_assembly.mgp_v5_snps_filter_PASS.with_contigs.with_ERCC.GATC.posarray.hdf5";
POSFN2="/hpc/hub_kind/koos/references/mouse/CAST_snps/posarray/CAST_EiJ.from_GRCm38_dna_primary_assembly.mgp_v5_snps_filter_PASS.with_ERCC.GATC.posarray.hdf5";

MIN_MAPQ=10;
MIN_EDITDIST=2;
UMILEN=6; #exact UMIlen does not matter; script only uses the fact that UMI is present
KEEPN=1;

tail -n+2 "$ANNOFN" | cut -f1,2,6 | while read indexnr bcnr cellcount; do
    printf -v fn1 "$FNFMT" $indexnr $bcnr "unique_${gt1}";
    printf -v fn2 "$FNFMT" $indexnr $bcnr "unique_${gt2}";
    printf -v fn3 "$FNFMT" $indexnr $bcnr "ambiguous_${gt1}_${gt2}";

    printf -v outfn "${OUTFMT}" $indexnr $bcnr;
    outbase="${outfn/%.{\name\}.top_n_\{keepn\}.event_counts.PE_R1.pos.hdf5}";
    invalid_outfn="${outbase}.PE_R1.invalid_pos.sorted.bam"
    count_err="${outbase}.PE_R1.count.err";

    if [ ! -f $count_err ]; then

        cmd="\
            python ./bin/process_damid_reads_and_generate_event_counts_from_split_gt.py \
                -vvv \
                --outfnfmt \"${outfn}\" \
                --genotype-names ${GENOTYPES} \
                --min-mapq $MIN_MAPQ \
                --umi-length $UMILEN \
                --keep-n $KEEPN \
                --min-editdistance $MIN_EDITDIST \
                --pos-file1 \"${POSFN1}\" \
                --pos-file2 \"${POSFN2}\" \
                \"${fn1}\"  \
                \"${fn2}\"  \
                \"${fn3}\"  \
                2> \"${count_err}\" \
        ;";
        echo $cmd;
        
    fi;
    
done > $filelist;

sbatch -J $LIMSID.${taskname} -o './joblog/%A.task_%a.out' -e './joblog/%A.task_%a.err' --time=3:00:00 --mem=15G --array=1-32 ./bin/run_preformatted_cmd_slurm.sh "$filelist" $(cat $filelist | wc -l);
#jobid 23817673


#delete ambigious DamID genotype, which will not be used
#rm data/aligned/*Dam*ambiguous*
#rm data/counts/*ambiguous*hdf5

## bin event counts - 100kb
# NO mappability file is provided, since this would filter out positions that are unmappable with just R1
# This would effectively remove some of the gains that we're hopefully getting from the PE alignments
taskname="bin_damid_event_counts_PE.F1xCAST";
filelist="./filelists/$(date +"%Y%m%d_%H%M").${taskname}.cmds.list";

POSFNS="\
129S1_SvImJ_snps\t/hpc/hub_kind/koos/references/mouse/129S1_SvImJ_snps/posarray/129S1_SvImJ.from_GRCm38_dna_primary_assembly.mgp_v5_snps_filter_PASS.with_contigs.with_ERCC.GATC.posarray.hdf5\nCAST_snps\t/hpc/hub_kind/koos/references/mouse/CAST_snps/posarray/CAST_EiJ.from_GRCm38_dna_primary_assembly.mgp_v5_snps_filter_PASS.with_ERCC.GATC.posarray.hdf5\
"

ANNOFN="./metadata/KIN8403_anno.tsv";
gt1="129S1_SvImJ_snps"
gt2="CAST_snps";
GENOTYPES="unique_${gt1} unique_${gt2}";

BINSIZE="100000"; #only 1 binsize can be provided 
INDIR="./data/counts/";
KEEPN=1;
INFNFMT="${INDIR}/${LIMSID}.index%02d.DamID2_BC_%03d.${gt1}_vs_${gt2}.%s.top_n_${KEEPN}.event_counts.PE_R1.pos.hdf5";

tail -n+2 "$ANNOFN" | cut -f1,2 | while read indexnr bcnr; do
    echo -e $POSFNS | cut -f1,2 | while read genotype posfn; do
        printf -v fn "$INFNFMT" $indexnr $bcnr "unique_${genotype}";
        outfn="${fn/%.event_counts.PE_R1.pos.hdf5/.event_counts.PE_R1.binsize_${BINSIZE}.hdf5}";
        
        if [ ! -f $outfn ]; then
            cmd="\
                bin_countfile.py \
                    -vvv \
                    --posfile \"${posfn}\" \
                    --binsize ${BINSIZE} \
                    --outfile \"${outfn}\" \
                    \"$fn\" \
            ;";
            echo "$cmd";
        fi;
    done;
done > "$filelist";

sbatch -J $LIMSID.${taskname} -o './joblog/%A.task_%a.out' -e './joblog/%A.task_%a.err' --time=2:00:00 --mem=10G --array=1-32 ./bin/run_preformatted_cmd_slurm.sh "$filelist" $(cat $filelist | wc -l);
#jobid 23837882


## extract read names of invalid_pos DamID reads 
taskname="extract_invalid_pos_readnames";
filelist="./filelists/$(date +"%Y%m%d_%H%M").${taskname}.cmds.list";

ANNOFN="./metadata/KIN8403_anno.tsv";
gt1="129S1_SvImJ_snps"
gt2="CAST_snps";
INDIR="./data/counts/";
OUTDIR="./data/aligned/";
INFMT="${INDIR}/${LIMSID}.index%02d.DamID2_BC_%03d.${gt1}_vs_${gt2}.PE_R1.invalid_pos.sorted.bam";
OUTFMT="${OUTDIR}/${LIMSID}.index%02d.DamID2_BC_%03d.${gt1}_vs_${gt2}.PE_R1.invalid_pos.readname.lst";

tail -n+2 "$ANNOFN" | cut -f1,2,6 | while read indexnr bcnr cellcount; do
    printf -v fn "${INFMT}" $indexnr $bcnr;
    printf -v outfn "${OUTFMT}" $indexnr $bcnr;
    
    cmd="\
        samtools view $fn | cut -f1 > $outfn;";

    echo "$cmd";
  
done > "$filelist";

sbatch -J $LIMSID.${taskname} -o './joblog/%A.task_%a.out' -e './joblog/%A.task_%a.err' --time=1:00:00 --mem=5G --array=1-16 ./bin/run_preformatted_cmd_slurm.sh "$filelist" $(cat $filelist | wc -l);
#jobid 20811850


## invalid_pos DamID reads to new fastq.gz using ``seqtk`` 
taskname="invalid_pos_reads_to_fastq";
filelist="./filelists/$(date +"%Y%m%d_%H%M").${taskname}.cmds.list";

ANNOFN="./metadata/KIN8403_anno.tsv"
SEQTK="/hpc/hub_kind/prullens/software/seqtk/seqtk";
gt1="129S1_SvImJ_snps"
gt2="CAST_snps";
INDIR="./data/aligned/";
LST="${INDIR}/${LIMSID}.index%02d.DamID2_BC_%03d.${gt1}_vs_${gt2}.PE_R1.invalid_pos.readname.lst";
INFMT_BASE="./data/demultiplexed/${LIMSID}.index%02d.DamID2_BC_%03d";
OUTFMT_BASE="./data/demultiplexed/${LIMSID}.index%02d.DamID2_BC_%03d.invalid_pos";

tail -n+2 $ANNOFN | cut -f1,2,6 | while read indexnr bcnr cellcount; do
    printf -v lst "${LST}" $indexnr $bcnr;
    
    for read in .R1.fastq.gz .R2.fastq.gz; do
        printf -v fn "${INFMT_BASE}${read}" $indexnr $bcnr;
        printf -v outfn "${OUTFMT_BASE}${read}" $indexnr $bcnr;
        
        cmd="\
            ${SEQTK} seq \
            ${fn} \
            | ${SEQTK} subseq - ${lst} \
            | gzip \
            > \"${outfn}\" \
        ;";

        echo $cmd;
    done;
done > "$filelist";

sbatch -J $LIMSID.${taskname} -o './joblog/%A.task_%a.out' -e './joblog/%A.task_%a.err' --time=2:00:00 --mem=10G --array=1-16 ./bin/run_preformatted_cmd_slurm.sh "$filelist" $(cat $filelist | wc -l);
#jobid 20813256

#delete readnames of ChIC reads
#rm data/aligned/*PE_R1.invalid_pos.readname.lst


## align ChIC data, without GA prefix - to both parental genotypes - PE: DamID invalid_pos reads are used
taskname="align_damid_data_all_genotypes_PE";
filelist="./filelists/$(date +"%Y%m%d_%H%M").${taskname}.cmds.list";

REFDATA="\
129S1_SvImJ_snps\t/hpc/hub_kind/koos/references/mouse/129S1_SvImJ_snps/bowtie2/129S1_SvImJ.from_GRCm38_dna_primary_assembly.mgp_v5_snps_filter_PASS.with_contigs.with_ERCC\nCAST_snps\t/hpc/hub_kind/koos/references/mouse/CAST_snps/bowtie2/CAST_EiJ.from_GRCm38_dna_primary_assembly.mgp_v5_snps_filter_PASS.with_ERCC\
"

ANNOFN="./metadata/KIN8403_anno.tsv";
OUTDIR="./data/aligned/";
FNFMT="./data/demultiplexed/${LIMSID}.index%02d.DamID2_BC_%03d.invalid_pos.R1.fastq.gz";
OUTFNFMT="${OUTDIR}/${LIMSID}.index%02d.BC_ChIC_%03d.%s.chic.PE.bam";
PREFIX="";

tail -n+2 $ANNOFN | cut -f1,2 | while read indexnr bcnr; do
    printf -v fn1 "$FNFMT" $indexnr $bcnr;
    fn2=$(echo $fn1 | sed -e 's/R1\.fastq/R2\.fastq/g');

    echo -e $REFDATA | cut -f1,2 | while read genotype bowtie2_index; do

        printf -v outfn "$OUTFNFMT" $indexnr $bcnr "$genotype";
        bowtie2_errfile="${outfn/%bam/bowtie2.err}";

        if [ ! -f $outfn ]; then
            cmd="\
                bowtie2 \
                    --seed 42 \
                    --very-sensitive -N 1 \
                    -x \"${bowtie2_index}\" \
                    -1 \"${fn1}\" \
                    -2 \"${fn2}\" \
                    2> \"${bowtie2_errfile}\" \
                | samtools view -b - \
                > \"${outfn}\" \
                ;";

            echo "$cmd";
        
        fi;
        
    done;
done > "$filelist";

sbatch -J $LIMSID.${taskname} -o './joblog/%A.task_%a.out' -e './joblog/%A.task_%a.err' --time=15:00:00 --mem=25G --array=1-42 ./bin/run_preformatted_cmd_slurm.sh "$filelist" $(cat $filelist | wc -l);
#jobid 20814264


## split PE alignments by genotype
taskname="split_genotype_ChIC_PE";
filelist="./filelists/$(date +"%Y%m%d_%H%M").${taskname}.cmds.list";

ANNOFN="./metadata/KIN8403_anno.tsv";
gt1="129S1_SvImJ_snps"
gt2="CAST_snps";
FNFMT="./data/aligned/${LIMSID}.index%02d.BC_ChIC_%03d.%s.chic.PE.bam";
OUTFNFMT="./data/aligned/${LIMSID}.index%02d.BC_ChIC_%03d.${gt1}_vs_${gt2}.{name}.chic.PE_R1.readname_sorted.bam";

tail -n+2 $ANNOFN | cut -f1,2 | while read indexnr bcnr; do
    printf -v fn1 "$FNFMT" $indexnr $bcnr "$gt1";
    printf -v fn2 "$FNFMT" $indexnr $bcnr "$gt2";
    printf -v outfmt "$OUTFNFMT" $indexnr $bcnr;
    errfn="${outfmt/%\{name\}.PE_R1.readname_sorted.bam/PE_R1.split_genotype.err}";
    
    cmd="\
        python ./bin/split_PE_DamID_alignments.py \
            -vvv \
            --unsorted-input \
            --name \"${gt1}\" --name \"${gt2}\" \
            --outfmt \"${outfmt}\" \
            \"${fn1}\" \"${fn2}\" \
        2> \"${errfn}\" \
    ;";

    echo "$cmd";
done > "$filelist";

sbatch -J $LIMSID.${taskname} -o './joblog/%A.task_%a.out' -e './joblog/%A.task_%a.err' --time=3:00:00 --mem=5G --array=1-32 ./bin/run_preformatted_cmd_slurm.sh "$filelist" $(cat $filelist | wc -l);
#jobid 20857195


## position sorting split genotype reads
taskname="position_sort_split_gt_DamID_alignments_PE";
filelist="./filelists/$(date +"%Y%m%d_%H%M").${taskname}.cmds.list";

ANNOFN="./metadata/KIN8403_anno.tsv";
gt1="129S1_SvImJ_snps"
gt2="CAST_snps";
GENOTYPES="unique_${gt1} unique_${gt2}";
INDIR="./data/aligned/";
FNFMT="${INDIR}/${LIMSID}.index%02d.BC_ChIC_%03d.${gt1}_vs_${gt2}.%s.chic.PE_R1.readname_sorted.bam";

tail -n+2 $ANNOFN | cut -f1,2 | while read indexnr bcnr; do
    for gt in $GENOTYPES; do
        printf -v fn "$FNFMT" $indexnr $bcnr "$gt";
        samtools_sort_tmp="${fn/%readname_sorted.bam/sort_tmp}";
        outfn="${fn/%readname_sorted.bam/sorted.bam}";

        if [ ! -f $outfn ]; then
            cmd="\
                samtools view -ub \"${fn}\" \
                | samtools sort -m 2500M -T \"${samtools_sort_tmp}\" -l 9 - \
                > \"${outfn}\" \
            ;";

            echo $cmd;
        fi;
    done;
done > "$filelist";

sbatch -J $LIMSID.${taskname} -o './joblog/%A.task_%a.out' -e './joblog/%A.task_%a.err' --time=1:00:00 --mem=8G --array=1-32 ./bin/run_preformatted_cmd_slurm.sh "$filelist" $(cat $filelist | wc -l);
#jobid 20857341


#delete original alignment and readname_sorted bam files
#rm data/aligned/*ChIC*snps.chic.PE.bam
#rm data/aligned/*chic.PE_R1.readname_sorted.bam


## bin ChIC data to parental genomes 
taskname="BinUMIUniqueGenomicSignal_genotypes";
filelist="./filelists/$(date +"%Y%m%d_%H%M").${taskname}.cmds.list";

gt1="129S1_SvImJ_snps";
gt2="CAST_snps";

ANNOFN="./metadata/KIN8403_anno.tsv";
LIMSID=$LIMSID;
INDIR="./data/aligned/";
OUTDIR="./data/counts/";
BINSIZE="100000";
FNFMT="${LIMSID}.index%02d.BC_ChIC_%03d.${gt1}_vs_${gt2}.%s.chic.PE_R1.sorted.bam";
OUTFNFMT="${LIMSID}.index%02d.BC_ChIC_%03d.${gt1}_vs_${gt2}.%s.top_n_%d.chic.PE_R1.AT_noTC.event_counts.binsize_%d.hdf5"
KEEPN=1;

tail -n+2 $ANNOFN | awk '($1!=9)' | cut -f1,2 | while read indexnr barcodenr; do
    for genotype in "unique_${gt1}" "unique_${gt2}"; do
        printf -v fn "${INDIR}/${FNFMT}" $indexnr $barcodenr "$genotype";
        printf -v outfn "${OUTDIR}/${OUTFNFMT}" $indexnr $barcodenr "$genotype" $KEEPN $BINSIZE;

        cmd="\
            python ./bin/bin_genomic_signal.py \
                --leading AT \
                --leading-excl TC \
                --binsize ${BINSIZE} \
                --keepn $KEEPN \
                --outfile \"${outfn}\" \
                \"${fn}\" \
        ;";  
        
        echo $cmd;
        
    done;
done > "$filelist";

sbatch -J $LIMSID.${taskname} -o './joblog/%A.task_%a.out' -e './joblog/%A.task_%a.err' --time=1:00:00 --mem=8G --array=1-32 ./bin/run_preformatted_cmd_slurm.sh "$filelist" $(cat $filelist | wc -l);
#jobid 20857427
