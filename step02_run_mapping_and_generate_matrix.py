# ——————————————————>>>>>>>>>>
# Project Name: Hi-C
# Author: Hua-nan ZHAO
# E-mail: hermanzhaozzzz@gmail.com
# Update log:
#     2022-07-25: start project
# ——————————————————>>>>>>>>>>
import os
import json
# ------------------------------------------------------------------->>>>>>>>>>
# FUNCTIONS
# ------------------------------------------------------------------->>>>>>>>>>
def print_head(SAMPLES, MODE):
    print('----------\nSAMPLES:')
    [print('\t' + i) for i in SAMPLES]
    print('----------\nMODE:')
    print('\t' + MODE)
    print('----------\n\n')

def check_cmd(x):
    return any(
        os.access(os.path.join(path, x), os.X_OK) 
        for path in os.environ["PATH"].split(os.pathsep)
    )

def check_read(x):
    if x == "PE":
        read = ['R1', 'R2']
    elif x == "SE":
        read = ['SE']
    else:
        raise ValueError()
    return read
# ------------------------------------------------------------------->>>>>>>>>>
# SAMPLE INFO
# ------------------------------------------------------------------->>>>>>>>>>
with open('./samples.json') as f:
    dt = json.loads(f.read())

SAMPLES = dt['samples']
MODE = dt['seq_mode']
READ = check_read(MODE)  # Hi-C默认双端PE，这里是一句废话，不过为了以后万一有单端，备用，就不改了

print_head(SAMPLES, MODE)
print(READ)
# ------------------------------------------------------------------->>>>>>>>>>
# RUN INFO
# ------------------------------------------------------------------->>>>>>>>>>
THREAD = dt['thread']
BOWTIE2_GLOBAL_OPTIONS = dt['bowtie2_params_global']
BOWTIE2_LOCAL_OPTIONS = dt['bowtie2_params_local']
LIGATION_SITE = dt['ligation_site']
SORT_RAM_PER_THREAD = dt['sort_ram_per_thread']
MIN_MAPQ = dt['min_mapq']
BIN_SIZES = dt['bin_sizes']
CHR_SIZES = dt['chr_sizes']
DIGEST_BED = dt['digest_bed']
RESTRICTION = dt['restriction_sequence']
DANGLING = dt['dangling_sequence']
# ------------------------------------------------------------------->>>>>>>>>>
# DATABASE INFO
# ------------------------------------------------------------------->>>>>>>>>>
BOWTIE2_INDEX = dt["bowtie2_index"]
# ------------------------------------------------------------------->>>>>>>>>>
# SOFTWARE INFO
# ------------------------------------------------------------------->>>>>>>>>>
# check if cmd exists
assert check_cmd("bowtie2")  # Bowtie 2 version 2.4.5
assert check_cmd("samtools")  # samtools 1.15.1 Using htslib 1.15.1
assert check_cmd("java")  # openjdk version 1.8.0_312
# hicexplorer=3.7.2
assert check_cmd("hicBuildMatrix") 
assert check_cmd("hicCorrectMatrix") 
assert check_cmd("hicConvertFormat") 
assert check_cmd("hicFindTADs") 

# pip install iced                           
# Collecting iced
#   Downloading iced-0.5.10.tar.gz (2.3 MB)
# manually set cmd path
BOWTIE2 = "bowtie2"
SAMTOOLS = "samtools"
JAVA = "java"
hicBuildMatrix = "hicBuildMatrix"
hicCorrectMatrix = "hicCorrectMatrix"
hicConvertFormat = "hicConvertFormat"
hicFindTADs = "hicFindTADs"


# ------------------------------------------------------------------->>>>>>>>>>
# rule all
# ------------------------------------------------------------------->>>>>>>>>>
rule all:
    input:
        expand("../fastq/{sample}_%s.fastq.gz" % READ[0], sample=SAMPLES),
        expand("../bam/{sample}/{sample}.merged_sortn.bwt2pairs.bam", sample=SAMPLES),
        expand("../valid_pairs/{sample}.merged_sortn.bwt2pairs.validPairs", sample=SAMPLES),
        expand("../valid_pairs/{sample}.rm_dup_pairs.allValidPairs", sample=SAMPLES),
        expand("../matrix/{bin_size}/{sample}_{bin_size}_raw.matrix", sample=SAMPLES, bin_size=BIN_SIZES),
        expand("../matrix/{bin_size}/{sample}_{bin_size}_iced.matrix", sample=SAMPLES, bin_size=BIN_SIZES),
        expand("../bam/{sample}/{sample}_R1.mapstat", sample=SAMPLES),
        expand("../bam/{sample}/{sample}_R2.mapstat", sample=SAMPLES),
        expand("../quality_checks/plotMapping_{sample}.pdf", sample=SAMPLES),
        expand("../quality_checks/plotMappingPairing_{sample}.pdf", sample=SAMPLES),
        expand("../quality_checks/plotHiCFragment_{sample}.pdf", sample=SAMPLES),
        expand("../quality_checks/plotHiCContactRanges_{sample}.pdf", sample=SAMPLES),
        expand("../hic_file/{sample}.rm_dup_pairs.allValidPairs.hic", sample=SAMPLES),
        expand("../matrix/{bin_size}/{sample}_{bin_size}_RawMatrix.h5", sample=SAMPLES, bin_size=BIN_SIZES),
        expand("../matrix/{bin_size}/{sample}_{bin_size}_KRjustify_Matrix.h5", sample=SAMPLES, bin_size=BIN_SIZES),
        # expand("../calling_use_hdf5/{bin_size}/{sample}_zscore_matrix.h5", sample=SAMPLES, bin_size=BIN_SIZES),
        
# ------------------------------------------------------------------->>>>>>>>>>
# mapping_global
# ------------------------------------------------------------------->>>>>>>>>>
rule map_read1_global:
    input: 
        "../fastq/{sample}_%s.fastq.gz" % READ[0]
    output: 
        sam=temp("../bam/{sample}/{sample}_R1.bwt2glob.sam"),
        un=temp("../bam/fastq_unmap_bwt2/{sample}_R1.fastq.gz")
    params:
        SE="../fastq/{sample}_SE.fastq.gz",  # 不存在，为以后smk升级所留
        R1="../fastq/{sample}_R1.fastq.gz",  # map read1
        R2="../fastq/{sample}_R2.fastq.gz",
        RG="SM:{sample}_R1"
    log:
        "../bam/{sample}/{sample}_R1.bwt2glob.log"
    shell:
        """
        INPUT={input}
        if [[ $INPUT =~ .*SE.fastq.gz$ ]]; then 
            echo "[FATAL] find SE reads, raw reads should be PE reads in Hi-C protocol" > {log}
        elif [[ $INPUT =~ .*R1.fastq.gz$ ]]; then
            echo "[DEBUG] find PE reads, go on mapping" > {log}
            {BOWTIE2} {BOWTIE2_GLOBAL_OPTIONS} --un {output.un} --rg-id BMG --rg {params.RG} -p {THREAD} -x {BOWTIE2_INDEX} -U {params.R1} -S {output.sam} >> {log} 2>&1
        else
            echo "[FATAL] fastq should be *.SE.fastq.gz or *.R[1,2].fastq.gz" > {log}
        fi
        echo "[DEBUG] bowtie2 mapping done" >> {log}
        """
rule map_read2_global:
    input: 
        "../fastq/{sample}_%s.fastq.gz" % READ[0]
    output: 
        sam=temp("../bam/{sample}/{sample}_R2.bwt2glob.sam"),
        un=temp("../bam/fastq_unmap_bwt2/{sample}_R2.fastq.gz")
    params:
        SE="../fastq/{sample}_SE.fastq.gz",  # 不存在，为以后smk升级所留
        R1="../fastq/{sample}_R1.fastq.gz",
        R2="../fastq/{sample}_R2.fastq.gz",  # map read2
        RG="SM:{sample}_R2"
    log:
        "../bam/{sample}/{sample}_R2.bwt2glob.log"
    shell:
        """
        INPUT={input}
        if [[ $INPUT =~ .*SE.fastq.gz$ ]]; then 
            echo "[FATAL] find SE reads, raw reads should be PE reads in Hi-C protocol" > {log}
        elif [[ $INPUT =~ .*R1.fastq.gz$ ]]; then
            echo "[DEBUG] find PE reads, go on mapping" > {log}
            {BOWTIE2} {BOWTIE2_GLOBAL_OPTIONS} --un {output.un} --rg-id BMG --rg {params.RG} -p {THREAD} -x {BOWTIE2_INDEX} -U {params.R2} -S {output.sam} >> {log} 2>&1
        else
            echo "[FATAL] fastq should be *.SE.fastq.gz or *.R[1,2].fastq.gz" > {log}
        fi
        echo "[DEBUG] bowtie2 mapping done" >> {log}
        """
rule sam2bam_read1_global:
    input:
        "../bam/{sample}/{sample}_R1.bwt2glob.sam"
    output:
        temp("../bam/{sample}/{sample}_R1.bwt2glob.bam")
    log:
        "../bam/{sample}/{sample}_R1.bwt2glob.log"
    shell:
        """
        echo "[DEBUG] start samtools view (sam to bam)" > {log}
        {SAMTOOLS} view -F 4 -hSb -@ {THREAD} {input} > {output}
        echo "[DEBUG] samtools view (sam to bam) done" >> {log}
        """
rule sam2bam_read2_global:
    input:
        "../bam/{sample}/{sample}_R2.bwt2glob.sam"
    output:
        temp("../bam/{sample}/{sample}_R2.bwt2glob.bam")
    log:
        "../bam/{sample}/{sample}_R2.bwt2glob.log"
    shell:
        """
        echo "[DEBUG] start samtools view (sam to bam)" > {log}
        {SAMTOOLS} view -F 4 -hSb -@ {THREAD} {input} > {output}
        echo "[DEBUG] samtools view (sam to bam) done" >> {log}
        """
# ------------------------------------------------------------------->>>>>>>>>>
# trim unmapped reads
# ------------------------------------------------------------------->>>>>>>>>>
rule hicpro_read1_trim:
    input:
        "../bam/fastq_unmap_bwt2/{sample}_R1.fastq.gz"
    output:
        temp("../bam/fastq_unmap_bwt2/{sample}_R1.trimmed.fastq.gz")
    params:
        FQ="../bam/fastq_unmap_bwt2/{sample}_R1.trimmed.fastq"
    log:
        "../bam/fastq_unmap_bwt2/{sample}_R1.trimmed.log"
    shell:
        """
        program/HiC-Pro_3.1.0/scripts/cutsite_trimming --fastq {input} --cutsite AAGCTAGCTT --out {params.FQ} > {log} 2>&1
        echo "[DEBUG] trimming done, start to gzip file" >> {log}
        pigz -p {THREAD} {params.FQ}
        echo "[DEBUG] gzip done" >> {log}
        """
rule hicpro_read2_trim:
    input:
        "../bam/fastq_unmap_bwt2/{sample}_R2.fastq.gz"
    output:
        temp("../bam/fastq_unmap_bwt2/{sample}_R2.trimmed.fastq.gz")
    params:
        FQ="../bam/fastq_unmap_bwt2/{sample}_R2.trimmed.fastq"
    log:
        "../bam/fastq_unmap_bwt2/{sample}_R2.trimmed.log"
    shell:
        """
        program/HiC-Pro_3.1.0/scripts/cutsite_trimming --fastq {input} --cutsite AAGCTAGCTT --out {params.FQ} > {log} 2>&1
        echo "[DEBUG] trimming done, start to gzip file" >> {log}
        pigz -p {THREAD} {params.FQ}
        echo "[DEBUG] gzip done" >> {log}
        """
# ------------------------------------------------------------------->>>>>>>>>>
# mapping_local
# ------------------------------------------------------------------->>>>>>>>>>
rule map_read1_local:
    input: 
        "../bam/fastq_unmap_bwt2/{sample}_R1.trimmed.fastq.gz"
    output: 
        temp("../bam/{sample}/{sample}_R1.unmap_bwt2loc.sam")
    params:
        RG="SM:{sample}_R1_unmap"
    log:
        temp("../bam/{sample}/{sample}_R1.unmap_bwt2loc.log")
    shell:
        """
        echo "[DEBUG] start bowtie2 mapping" > {log}
        {BOWTIE2} {BOWTIE2_LOCAL_OPTIONS} --rg-id BML --rg {params.RG} -p {THREAD} -x {BOWTIE2_INDEX} -U {input} -S {output} >> {log} 2>&1
        echo "[DEBUG] bowtie2 mapping done" >> {log}
        """
rule map_read2_local:
    input: 
        "../bam/fastq_unmap_bwt2/{sample}_R2.trimmed.fastq.gz"
    output: 
        temp("../bam/{sample}/{sample}_R2.unmap_bwt2loc.sam")
    params:
        RG="SM:{sample}_R2_unmap"
    log:
        temp("../bam/{sample}/{sample}_R2.unmap_bwt2loc.log")
    shell:
        """
        echo "[DEBUG] start bowtie2 mapping" > {log}
        {BOWTIE2} {BOWTIE2_LOCAL_OPTIONS} --rg-id BML --rg {params.RG} -p {THREAD} -x {BOWTIE2_INDEX} -U {input} -S {output} >> {log} 2>&1
        echo "[DEBUG] bowtie2 mapping done" >> {log}
        """
rule sam2bam_read1_local:
    input:
        "../bam/{sample}/{sample}_R1.unmap_bwt2loc.sam"
    output:
        temp("../bam/{sample}/{sample}_R1.unmap_bwt2loc.bam")
    log:
        "../bam/{sample}/{sample}_R1.unmap_bwt2loc.log"
    shell:
        """
        echo "[DEBUG] start samtools view (sam to bam)" > {log}
        {SAMTOOLS} view -hSb -@ {THREAD} {input} > {output}
        echo "[DEBUG] samtools view (sam to bam) done" >> {log}
        """
rule sam2bam_read2_local:
    input:
        "../bam/{sample}/{sample}_R2.unmap_bwt2loc.sam"
    output:
        temp("../bam/{sample}/{sample}_R2.unmap_bwt2loc.bam")
    log:
        "../bam/{sample}/{sample}_R2.unmap_bwt2loc.log"
    shell:
        """
        echo "[DEBUG] start samtools view (sam to bam)" > {log}
        {SAMTOOLS} view -hSb -@ {THREAD} {input} > {output}
        echo "[DEBUG] samtools view (sam to bam) done" >> {log}
        """
# ------------------------------------------------------------------->>>>>>>>>>
# merge_global_and_local_bam
# ------------------------------------------------------------------->>>>>>>>>>
rule merge_global_and_local_read1:
    input:
        bam_global = "../bam/{sample}/{sample}_R1.bwt2glob.bam",
        bam_local = "../bam/{sample}/{sample}_R1.unmap_bwt2loc.bam"
    output:
        temp("../bam/{sample}/{sample}_R1.merged.bam")
    log:
        "../bam/{sample}/{sample}_R1.merged.log"
    shell:
        """
        echo "[DEBUG] start samtools merge" > {log}
        {SAMTOOLS} merge -@ {THREAD} -n -f {output} {input.bam_global} {input.bam_local} >> {log} 2>&1
        echo "[DEBUG] samtools merge done" >> {log}
        """
rule merge_global_and_local_read2:
    input:
        bam_global = "../bam/{sample}/{sample}_R2.bwt2glob.bam",
        bam_local = "../bam/{sample}/{sample}_R2.unmap_bwt2loc.bam"
    output:
        temp("../bam/{sample}/{sample}_R2.merged.bam")
    log:
        "../bam/{sample}/{sample}_R2.merged.log"
    shell:
        """
        echo "[DEBUG] start samtools merge" > {log}
        {SAMTOOLS} merge -@ {THREAD} -n -f {output} {input.bam_global} {input.bam_local} >> {log} 2>&1
        echo "[DEBUG] samtools merge done" >> {log}
        """
# ------------------------------------------------------------------->>>>>>>>>>
# sort bam by name
# ------------------------------------------------------------------->>>>>>>>>>
rule sortn_read1_bam:
    input:
        "../bam/{sample}/{sample}_R1.merged.bam"
    output:
        temp("../bam/{sample}/{sample}_R1.merged_sortn.bam")
    params:
        temp_path="../temp_files/{sample}/{sample}_R1"
    log:
        "../bam/{sample}/{sample}_R1.merged_sortn.log"
    shell:
        """
        mkdir -p {params.temp_path}
        echo "[DEBUG] start samtools sort -n" > {log}
        {SAMTOOLS} sort -@ {THREAD} {SORT_RAM_PER_THREAD} -n -T {params.temp_path} -o {output} {input} >> {log} 2>&1
        echo "[DEBUG] samtools sort -n done" >> {log}
        """
rule sortn_read2_bam:
    input:
        "../bam/{sample}/{sample}_R2.merged.bam"
    output:
        temp("../bam/{sample}/{sample}_R2.merged_sortn.bam")
    params:
        temp_path="../temp_files/{sample}/{sample}_R2"
    log:
        "../bam/{sample}/{sample}_R2.merged_sortn.log"
    shell:
        """
        mkdir -p {params.temp_path}
        echo "[DEBUG] start samtools sort -n" > {log}
        {SAMTOOLS} sort -@ {THREAD} {SORT_RAM_PER_THREAD} -n -T {params.temp_path} -o {output} {input} >> {log} 2>&1
        echo "[DEBUG] samtools sort -n done" >> {log}
        """
# ------------------------------------------------------------------->>>>>>>>>>
# form pairs.bam
# ------------------------------------------------------------------->>>>>>>>>>
rule form_pairs_bam:
    input:
        bam_f="../bam/{sample}/{sample}_R1.merged_sortn.bam",
        bam_r="../bam/{sample}/{sample}_R2.merged_sortn.bam"
    output:
        "../bam/{sample}/{sample}.merged_sortn.bwt2pairs.bam"
    log:
        "../bam/{sample}/{sample}.merged_sortn.bwt2pairs.log"
    shell:
        # Usage : python mergeSAM.py
        # -f/--forward <forward read mapped file>
        # -r/--reverse <reverse read mapped file>
        # [-o/--output] <Output file. Default is stdin>
        # [-s/--single] <report singleton>
        # [-m/--multi] <report multiple hits>
        # [-q/--qual] <minimum reads mapping quality>
        # [-t/--stat] <generate a stat file>
        # [-v/--verbose] <Verbose>
        # [-h/--help] <Help>
        """
        echo "[DEBUG] start mergeSAM.py and generate pairs.bam" > {log}
        echo "[DEBUG] attention！！！[E::idx_find_and_load] is fine because of sorting name method" >> {log}
        python program/HiC-Pro_3.1.0/scripts/mergeSAM.py -q {MIN_MAPQ} -t -v -f {input.bam_f} -r {input.bam_r} -o {output} >> {log} 2>&1
        echo "[DEBUG] mergeSAM.py done" >> {log}
        """
# ------------------------------------------------------------------->>>>>>>>>>
# pairs.bam to validPairs
# ------------------------------------------------------------------->>>>>>>>>>
rule form_valid_pairs:
    input:
        "../bam/{sample}/{sample}.merged_sortn.bwt2pairs.bam"
    output:
        "../valid_pairs/{sample}.merged_sortn.bwt2pairs.validPairs"
    params:
        out_path="../valid_pairs/"
    log:
        "../valid_pairs/{sample}.log"
    shell:
        """
        echo "[DEBUG] start forming validPairs" > {log}
        echo "[DEBUG] attention！！！[E::idx_find_and_load] is fine because of sorting name method" >> {log}
        python program/HiC-Pro_3.1.0/scripts/mapped_2hic_dnase.py -v -a -r {input} -o {params.out_path} >> {log} 2>&1
        echo "[DEBUG] forming validPairs done" >> {log}
        """
# ------------------------------------------------------------------->>>>>>>>>>
# rm duplicates of validPairs, form allValidPairs
# ------------------------------------------------------------------->>>>>>>>>>
rule rm_dup_form_allValidPairs:
    input:
        "../valid_pairs/{sample}.merged_sortn.bwt2pairs.validPairs"
    output:
        "../valid_pairs/{sample}.rm_dup_pairs.allValidPairs"
    params:
        awk_begion="""{c1=0;c2=0;s1=0;s2=0}""",
        awk_body="""{print;c1=$2;c2=$5;s1=$3;s2=$6}"""
    log:
        "../valid_pairs/{sample}.rm_dup_pairs.log"
    shell:
        """
        echo "[DEBUG] start merge valid interactions and remove duplicates" > {log}
        LANG=en;
        sort -T tmp -S 50% -k2,2V -k3,3n -k5,5V -k6,6n -m {input} | \
        awk -F"\t" 'BEGIN{params.awk_begion}(c1!=$2 || c2!=$5 || s1!=$3 || s2!=$6){params.awk_body}' \
        > {output}
        echo "[DEBUG] merge valid interactions and remove duplicates done" >> {log}
        """
# ------------------------------------------------------------------->>>>>>>>>>
# generate contact maps at different resolution (bin_size)
# ------------------------------------------------------------------->>>>>>>>>>     
rule generate_contact_maps:
    input:
        "../valid_pairs/{sample}.rm_dup_pairs.allValidPairs"
    output:
        "../matrix/{bin_size}/{sample}_{bin_size}_raw.matrix"
    params:
        out_path="../matrix/{bin_size}/{sample}_{bin_size}_raw",
        bin_size="{bin_size}"
    log:
        "../matrix/{bin_size}/{sample}_{bin_size}_raw.log"
    shell:
        # usage: build_matrix --binsize BINSIZE|--binfile --chrsizes FILE --ifile FILE 
        #     --oprefix PREFIX [--binadjust] [--step STEP] [--binoffset OFFSET]
        #     [--matrix-format asis|upper|lower|complete][--chrA CHR... --chrB CHR...] [--quiet] [--progress] [--detail-progress]
        # 
        # usage: build_matrix --version
        # usage: build_matrix --help
        """
        echo "[DEBUG] start generate contact maps at {params.bin_size} resolution" > {log}
        cat {input} | \
        program/HiC-Pro_3.1.0/scripts/build_matrix \
        --matrix-format upper \
        --binsize {params.bin_size} \
        --chrsizes {CHR_SIZES} --ifile /dev/stdin \
        --oprefix {params.out_path}
        echo "[DEBUG] generate contact maps at {params.bin_size} resolution done" >> {log}
        """
# ------------------------------------------------------------------->>>>>>>>>>
# ice justify for raw matrix
        # #######################################################################
        # ## Normalization
        # #######################################################################
        # # Maximum number of iteration for ICE normalization. 
        # # Default: 100
        # MAX_ITER = 100

        # # Define which pourcentage of bins with low counts should be force to zero. 
        # # Default: 0.02. 
        # # Replace SPARSE_FILTERING
        # FILTER_LOW_COUNT_PERC = 0.02

        # # Define which pourcentage of bins with low counts should be discarded 
        # # before normalization. 
        # # Default: 0
        # FILTER_HIGH_COUNT_PERC = 0

        # # The relative increment in the results before declaring convergence. 
        # # Default: 0.1
        # EPS = 0.1
# ------------------------------------------------------------------->>>>>>>>>>    
rule ice_justify:
    input:
        "../matrix/{bin_size}/{sample}_{bin_size}_raw.matrix"
    output:
        "../matrix/{bin_size}/{sample}_{bin_size}_iced.matrix"
    params:
        bin_size="{bin_size}"
    log:
        "../matrix/{bin_size}/{sample}_{bin_size}_iced.log"
    shell:
        """
        echo "[DEBUG] start ice justification at {params.bin_size} resolution" > {log}
        program/HiC-Pro_3.1.0/scripts/ice \
            --results_filename {output} \
            --filter_low_counts_perc 0.02 \
            --filter_high_counts_perc 0 \
            --max_iter 100 \
            --eps 0.1 \
            --remove-all-zeros-loci \
            --base 0 \
            {input} >> {log} 2>&1
        echo "[DEBUG] ice justification at {params.bin_size} resolution done" >> {log}
        """
# ------------------------------------------------------------------->>>>>>>>>>
# form mapstat
# ------------------------------------------------------------------->>>>>>>>>>
rule form_mapstat:
    input:
        bam_f="../bam/{sample}/{sample}_R1.merged_sortn.bam",
        bam_r="../bam/{sample}/{sample}_R2.merged_sortn.bam",
        bam_gf="../bam/{sample}/{sample}_R1.bwt2glob.bam",
        bam_gr="../bam/{sample}/{sample}_R2.bwt2glob.bam",
        bam_lf="../bam/{sample}/{sample}_R1.unmap_bwt2loc.bam",
        bam_lr="../bam/{sample}/{sample}_R2.unmap_bwt2loc.bam"
    output:
        mapstat_f="../bam/{sample}/{sample}_R1.mapstat",
        mapstat_r="../bam/{sample}/{sample}_R2.mapstat"
    log:
        f="../bam/{sample}/{sample}_R1.map_count.log",
        r="../bam/{sample}/{sample}_R2.map_count.log",
    shell:
        """
        echo "[DEBUG] start formming mapstat" > {log.f}
        
        TOTAL_R1=`{SAMTOOLS} view -c {input.bam_f}`
        MAPPED_R1=`{SAMTOOLS} view -c -F 4 {input.bam_f}`
        GLOBAL_R1=`{SAMTOOLS} view -c -F 4 {input.bam_gf}`
        LOCAL_R1=`{SAMTOOLS} view -c -F 4 {input.bam_lf}`
        echo "total_R1\t" $TOTAL_R1 >> {output.mapstat_f}
        echo "mapped_R1\t" $MAPPED_R1 >> {output.mapstat_f}
        echo "global_R1\t" $GLOBAL_R1 >> {output.mapstat_f}
        echo "local_R1\t" $LOCAL_R1 >> {output.mapstat_f}
        
        echo "[DEBUG] formming mapstat done" >> {log.f}
        
        echo "[DEBUG] start formming mapstat" > {log.r}
        
        TOTAL_R2=`{SAMTOOLS} view -c {input.bam_r}`
        MAPPED_R2=`{SAMTOOLS} view -c -F 4 {input.bam_r}`
        GLOBAL_R2=`{SAMTOOLS} view -c -F 4 {input.bam_gr}`
        LOCAL_R2=`{SAMTOOLS} view -c -F 4 {input.bam_lr}`
        echo "total_R2\t" $TOTAL_R2 >> {output.mapstat_r}
        echo "mapped_R2\t" $MAPPED_R2 >> {output.mapstat_r}
        echo "global_R2\t" $GLOBAL_R2 >> {output.mapstat_r}
        echo "local_R2\t" $LOCAL_R2 >> {output.mapstat_r}
        
        echo "[DEBUG] formming mapstat done" >> {log.r}
        """
rule form_all_valid_pairs_mergestat:
    input:
        vp="../valid_pairs/{sample}.merged_sortn.bwt2pairs.validPairs",
        vp_rmdup="../valid_pairs/{sample}.rm_dup_pairs.allValidPairs"
    output:
        "../valid_pairs/{sample}.rm_dup_pairs.allValidPairs.mergestat"
    params:
        awk_begion=r"""{cis=0;trans=0;sr=0;lr=0}""",
        awk_begion2=r"""{cis=cis+1; d=$6>$3?$6-$3:$3-$6; if (d<=20000){sr=sr+1}else{lr=lr+1}}""",
        awk_begion3=r"""{trans=trans+1}""",
        awk_end=r"""{print "trans_interaction\t"trans"\ncis_interaction\t"cis"\ncis_shortRange\t"sr"\ncis_longRange\t"lr}"""
    shell:
        """
        allcount=$(cat {input.vp} | wc -l)
        allcount_rmdup=$(cat {input.vp_rmdup} | wc -l)
        echo -e "valid_interaction\t"$allcount > {output}
        echo -e "valid_interaction_rmdup\t"$allcount_rmdup >> {output}
        awk 'BEGIN{params.awk_begion} $2 == $5{params.awk_begion2} $2!=$5{params.awk_begion3}END{params.awk_end}' {input} >> {output}
        """
rule quality_checks:
    input:
        mapstat_f="../bam/{sample}/{sample}_R1.mapstat",
        mapstat_r="../bam/{sample}/{sample}_R2.mapstat",
        mergestat="../valid_pairs/{sample}.rm_dup_pairs.allValidPairs.mergestat"
    output:
        "../quality_checks/plotMapping_{sample}.pdf",
        "../quality_checks/plotMappingPairing_{sample}.pdf",
        "../quality_checks/plotHiCFragment_{sample}.pdf",
        "../quality_checks/plotHiCContactRanges_{sample}.pdf",
    params:
        pic_dir="../quality_checks/",
        bwt_dir="../bam/{sample}/",
        hic_dir="../valid_pairs/",
        sample_name="{sample}"
    log:
        plotMapping="../quality_checks/plotMapping_{sample}.log",
        plotMappingPairing="../quality_checks/plotMappingPairing_{sample}.log",
        plotHiCFragment="../quality_checks/plotHiCFragment_{sample}.log",
        plotHiCContactRanges="../quality_checks/plotHiCContactRanges_{sample}.log"
    shell:
        """
        mkdir -p {params.pic_dir}
        
        # plotMapping
        echo "[DEBUG] start Quality checks - Mapping results ..." > {log.plotMapping}
        R CMD BATCH --no-save --no-restore \
            "--args picDir='{params.pic_dir}' bwtDir='{params.bwt_dir}' sampleName='{params.sample_name}' r1tag='_R1' r2tag='_R2'" \
            program/HiC-Pro_3.1.0/scripts/plot_mapping_portion.R {log.plotMapping}
        echo "[DEBUG] Quality checks - Mapping results done" >> {log.plotMapping}
        
        echo "[DEBUG] start Quality checks - Pairing results ..." > {log.plotMappingPairing}
        R CMD BATCH --no-save --no-restore \
            "--args picDir='{params.pic_dir}' bwtDir='{params.bwt_dir}' sampleName='{params.sample_name}' rmMulti='1' rmSingle='1'" \
            program/HiC-Pro_3.1.0/scripts/plot_pairing_portion.R {log.plotMappingPairing}
        echo "[DEBUG] Quality checks - Pairing results done" >> {log.plotMappingPairing}
        
        # plotHiCFragment [RSstat file]
        echo "[DEBUG] start Quality checks - Hi-C processing ..." > {log.plotHiCFragment}
        R CMD BATCH --no-save --no-restore \
            "--args picDir='{params.pic_dir}' hicDir='{params.hic_dir}' sampleName='{params.sample_name}'" \
            program/HiC-Pro_3.1.0/scripts/plot_hic_fragment.R {log.plotHiCFragment}
        echo "[DEBUG] Quality checks - Hi-C processing done" >> {log.plotHiCFragment}
        
        # plotHiCContactRanges
        echo "[DEBUG] start Quality checks - Hi-C contact maps ..." > {log.plotHiCContactRanges}
        R CMD BATCH --no-save --no-restore \
            "--args picDir='{params.pic_dir}' hicDir='{params.hic_dir}' sampleName='{params.sample_name}'" \
            program/HiC-Pro_3.1.0/scripts/plot_hic_contacts.R {log.plotHiCContactRanges}
        echo "[DEBUG] Quality checks - Hi-C contact maps done" >> {log.plotHiCContactRanges}
        """
# ------------------------------------------------------------------->>>>>>>>>>
# convert valid pairs to Juicer .hic file
        # usage : hicpro2juicebox -i VALIDPAIRS -g GSIZE -j JUICERJAR [-r RESFRAG] [-t TEMP] [-o OUT] [-h]
        # Use option -h|--help for more information

        # Generate JuiceBox input file from HiC-Pro results
        # See http://www.aidenlab.org/juicebox/ for details about Juicebox
        # ---------------
        # OPTIONS

        #    -i|--input VALIDPAIRS : allValidPairs file generated by HiC-Pro >= 2.7.5
        #    -g|--gsize GSIZE : genome size file used during HiC-Pro processing
        #    -j|--jar JUICERJAR : path to juicebox_clt.jar file
        #    [-r|--resfrag] RESFRAG : restriction fragment file used by HiC-Pro
        #    [-t|--temp] TEMP : path to tmp folder. Default is current path
        #    [-o|--out] OUT : output path. Default is current path
        #    [-h|--help]: help

        # java -Xms16g -Xmx16g -Xmn8g -jar program/juicer_tools_1.22.01.jar pre -h
        # 如遇问题看一下这个 ![](https://tva1.sinaimg.cn/large/e6c9d24ely1h4lxwl997zj20nc0eztbi.jpg)
# ------------------------------------------------------------------->>>>>>>>>>
rule valid_pairs_to_hic:
    input:
        vp_rmdup="../valid_pairs/{sample}.rm_dup_pairs.allValidPairs"
    output:
        "../hic_file/{sample}.rm_dup_pairs.allValidPairs.hic"
    params:
        temp_path="../temp_files/{sample}/{sample}_juicer",
        hic_dir="../hic_file/"
    log:
        "../hic_file/{sample}.rm_dup_pairs.allValidPairs.hic.log"
    shell:
        """
        mkdir -p {params.temp_path}
        echo "[DEBUG] start valid pairs convertion to hic file ..." > {log}
        program/HiC-Pro_3.1.0/bin/utils/hicpro2juicebox.sh \
            -i {input.vp_rmdup} \
            -g {CHR_SIZES} \
            -j program/juicer_tools_1.22.01.jar \
            -r {DIGEST_BED} \
            -t {params.temp_path} \
            -o {params.hic_dir} >> {log} 2&>1
        echo "[DEBUG] valid pairs convertion to hic file done" >> {log}
        """
# ------------------------------------------------------------------->>>>>>>>>>
# make hicexplorer matrix
        # hicBuildMatrix -h
        # usage: hicBuildMatrix --samFiles two sam files two sam files --outFileName
        #                       FILENAME --QCfolder FOLDER --restrictionCutFile BED file
        #                       [BED file ...] --restrictionSequence RESTRICTIONSEQUENCE
        #                       [RESTRICTIONSEQUENCE ...] --danglingSequence
        #                       DANGLINGSEQUENCE [DANGLINGSEQUENCE ...]
        #                       [--outBam bam file] [--binSize BINSIZE [BINSIZE ...]]
        #                       [--minDistance MINDISTANCE] [--maxDistance MAXDISTANCE]
        #                       [--maxLibraryInsertSize MAXLIBRARYINSERTSIZE]
        #                       [--genomeAssembly GENOMEASSEMBLY]
        #                       [--region CHR:START-END] [--keepSelfLigation]
        #                       [--keepSelfCircles]
        #                       [--minMappingQuality MINMAPPINGQUALITY]
        #                       [--threads THREADS] [--inputBufferSize INPUTBUFFERSIZE]
        #                       [--doTestRun] [--doTestRunLines DOTESTRUNLINES]
        #                       [--skipDuplicationCheck] [--chromosomeSizes txt file]
        #                       [--help] [--version]
        # need about 150GB memory 
# ------------------------------------------------------------------->>>>>>>>>>
# rule build_hicexplorer_matrix:
#     input:
#         bam_f="../bam/{sample}/{sample}_R1.merged_sortn.bam",
#         bam_r="../bam/{sample}/{sample}_R2.merged_sortn.bam"
#     output:
#         "../hdf5_files/{sample}_RawMatrix.h5"
#     params:
#         # bins="20000 40000 100000 200000 500000 1000000"
#         bins=" ".join([str(bin) for bin in BIN_SIZES]),
#         qc_dir="../quality_checks_by_hicexplorer"
#     log:
#         "../hdf5_files/{sample}_RawMatrix.log"
#     shell:
#         """
#         mkdir -p {params.qc_dir}
#         echo "[DEBUG] start to build raw matrix.hdf5 file ..." > {log}
#         {hicBuildMatrix} --samFiles {input.bam_f} {input.bam_r} \
#             --outFileName {output} \
#             --binSize {params.bins} \
#             --threads {THREAD} \
#             --QCfolder {params.qc_dir} \
#             --restrictionCutFile {DIGEST_BED} \
#             --restrictionSequence {RESTRICTION} \
#             --danglingSequence {DANGLING} \
#             --minMappingQuality 20 \
#             --chromosomeSizes {CHR_SIZES} >> {log} 2&>1
#         echo "[DEBUG] build raw matrix.hdf5 file done" >> {log}
#         """
rule valid_pairs_to_hdf5:
    input:
        matrix="../matrix/{bin_size}/{sample}_{bin_size}_raw.matrix",
        bed="../matrix/{bin_size}/{sample}_{bin_size}_raw_abs.bed"
    output:
        "../matrix/{bin_size}/{sample}_{bin_size}_RawMatrix.h5"
    params:
        bin_size="{bin_size}"
    log:
        "../matrix/{bin_size}/{sample}_{bin_size}_RawMatrix.log"
    shell:
        """
        echo "[DEBUG] start to build raw matrix.hdf5 file ..." > {log}
        export HDF5_USE_FILE_LOCKING='FALSE' # for [locking disabled on this file system] err
        {hicConvertFormat} \
            -m {input.matrix} \
            --bedFileHicpro {input.bed} \
            --inputFormat hicpro \
            --outputFormat h5 -o {output} \
            --resolutions {params.bin_size} >> {log} 2&>1
        echo "[DEBUG] build raw matrix.hdf5 file done" >> {log}
        
        

        """
rule hdf5_matrix_correction:
    input:
        "../matrix/{bin_size}/{sample}_{bin_size}_RawMatrix.h5"
    output:
        "../matrix/{bin_size}/{sample}_{bin_size}_KRjustify_Matrix.h5"
    log:
        "../matrix/{bin_size}/{sample}_{bin_size}_KRjustify_Matrix.log"
    shell:
        """
        echo "[DEBUG] start to justify matrix ..." > {log}
        export HDF5_USE_FILE_LOCKING='FALSE' # for [locking disabled on this file system] err
        {hicCorrectMatrix} correct \
            -m {input} \
            --correctionMethod KR \
            -o {output} >> {log} 2>&1
        echo "[DEBUG] justify matrix done" >> {log}
        """
rule calling_use_hdf5:
    input:
        "../matrix/{bin_size}/{sample}_{bin_size}_KRjustify_Matrix.h5"
    output:
        "../calling_use_hdf5/{bin_size}/{sample}_zscore_matrix.h5"
    params:
        out_dir="../calling_use_hdf5/{bin_size}/{sample}"
    log:
        "../calling_use_hdf5/{bin_size}/{sample}.log"
    shell:
        # 内存开销非常大！！
        """
        echo "[DEBUG] start to calling ..." > {log}
        export HDF5_USE_FILE_LOCKING='FALSE' # for [locking disabled on this file system] err
        {hicFindTADs} \
            -m {input} \
            --outPrefix {params.out_dir} \
            --numberOfProcessors {THREAD} \
            --correctForMultipleTesting fdr >> {log} 2>&1
        echo "[DEBUG] calling done" >> {log}
        """
# # plot 
# hicPlotTADs --tracks 20211101-DdCBE_plot_TAD.ini -o out_image/test.pdf --region chrX:60000000-80000000 &
# hicPlotTADs --tracks 20211101-DdCBE_plot_TAD.ini -o out_image/test_2.pdf --region chrX:60000000-80000000 &
# hicPlotTADs --tracks 20211101-DdCBE_plot_TAD.ini -o out_image/test_3.pdf --region chrX:60000000-80000000 &


# call loops 
# java -Xmx32g -jar /home/menghaowei/menghw_HD/software_package/juicer/juicer_tools_1.22.01.jar hiccups \
# --threads 10 \
# -r 25000,10000,5000 \
# -k KR \
# -f 0.1,0.1,0.1 \
# -p 1,2,4 \
# -d 20000,20000,20000 \
# test.allValidPairs.hic \


# java -Xmx64g -jar /home/menghaowei/menghw_HD/software_package/juicer/juicer_tools_1.22.01.jar hiccups \
# --threads 20 \
# -r 5000,10000,25000 \
# -c chr1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X \
# --ignore-sparsity \
# 293T_WT.allValidPairs.NoFrag.hic all_hiccups_loops.NoFrag.hic


# java -Xmx64g -jar /home/menghaowei/menghw_HD/software_package/juicer/juicer_tools_1.22.01.jar hiccups \
# --threads 20 \
# -r 5000,10000,25000 \
# -c chr1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X \
# --ignore-sparsity \
# 293T_WT.allValidPairs.hic all_hiccups_loops 

# # extract data 
# java -Xmx32g -jar /home/menghaowei/menghw_HD/software_package/juicer/juicer_tools_1.22.01.jar dump \
# observed KR 293T_WT.allValidPairs.NoFrag.hic chr1 chr1 BP 100000 test.txt 
	

# # java -Xmx32g -jar /home/menghaowei/menghw_HD/software_package/juicer/juicer_tools_1.22.01.jar hiccups \
# # -r 5000 -c chr1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X --ignore-sparsity \
# # 293T_WT.allValidPairs.NoFrag.hic all_hiccups_loops

# # # on Lab iMac
# # java -Xmx32g -jar /Users/meng/juicer_tools_1.22.01.jar hiccups \
# # -r 5000,10000,25000 -c chr1,2,3,4,5,6,7,8,9,10 293T_WT.allValidPairs.hic all_hiccups_loops


# # cuda install 
# https://developer.nvidia.com/cuda-10.2-download-archive?target_os=MacOSX&target_arch=x86_64&target_version=1013&target_type=dmglocal



#### call TAD
# hicFindTADs -m 293T_WT.KRCorrectMatrix.h5 --outPrefix 293T_WT.KR --numberOfProcessors 16 --correctForMultipleTesting fdr > 293T_WT.hicFindTADs.log 2>&1 &  
# find TAD 
# srun -T 24 hicFindTADs -m 293T_WT.KR.CorrectMatrix.25Kbp.h5 --outPrefix 293T_WT.KR.25Kbp --numberOfProcessors 24 --correctForMultipleTesting fdr > 293T_WT.hicFindTADs.25Kbp.log 2>&1 &  

# srun -T 24 hicFindTADs -m 293T_WT.KR.CorrectMatrix.50Kbp.h5 --outPrefix 293T_WT.KR.50Kbp --numberOfProcessors 24 --correctForMultipleTesting fdr > 293T_WT.KR.hicFindTADs.50Kbp.log 2>&1 &  

# srun -T 24 hicFindTADs -m 293T_WT.ICE.CorrectMatrix.50Kbp.h5 --outPrefix 293T_WT.ICE.50Kbp --numberOfProcessors 24 --correctForMultipleTesting fdr > 293T_WT.hicFindTADs.50Kbp.log 2>&1 &  

# hicPlotTADs --tracks plot_TAD.DdCBEOnly.ini -o out_image/293T_WT.KR.BinSize5Kbp.chr12_96M_126M.DdCBEOnly.pdf --region chr12:96000000-126000000 &

# hicPlotTADs --tracks plot_TAD.DdCBEOnly.25Kbp.ini -o out_image/293T_WT.KR.BinSize25Kbp.chr12_96M_126M.DdCBEOnly.pdf --region chr12:96000000-126000000 &

# hicPlotTADs --tracks plot_TAD.DdCBEOnly.25Kbp.ini -o out_image/293T_WT.KR.BinSize25Kbp.chr12_96M_116M.DdCBEOnly.pdf --region chr12:96000000-116000000 &

# hicPlotTADs --tracks plot_TAD.Reds.DdCBEOnly.25Kbp.ini -o out_image/293T_WT.KR.BinSize25Kbp.chr12_96M_126M.DdCBEOnly.Reds.pdf --region chr12:96000000-126000000 &

# hicPlotTADs --tracks plot_TAD.DdCBEOnly.25Kbp.ini -o out_image/293T_WT.KR.BinSize25Kbp.chr12_10M_130M.DdCBEOnly.pdf --region chr12:10000000-130000000 &


# # plot TAD score 
# sort -k1,1 -k2,2n 293T_WT.KR_score.bedgraph > 293T_WT.KR_score.sort.bedgraph &

# sort -k1,1 -k2,2n 293T_WT.KR.25Kbp_score.bedgraph > 293T_WT.KR.25Kbp_score.sort.bedgraph &

# # bedgraph to bigwig
# bedGraphToBigWig 293T_WT.KR_score.sort.bedgraph ~/menghw_HD/reference/hg38.only_chrom.sizes 293T_WT.KR_score.bigwig &

# bedGraphToBigWig 293T_WT.KR.25Kbp_score.sort.bedgraph ~/menghw_HD/reference/hg38.only_chrom.sizes 293T_WT.KR.25Kbp_score.bigwig &

# # region in /home/menghaowei/menghw_HD/DdCBE_project/08.DdCBE_merge_all/region_cor_analysis/peak_region

# # plot TAD
# cd /home/menghaowei/menghw_HD/DdCBE_project/10.hic_data/01.hicpro_WT/hicexp_result

# computeMatrix reference-point -S \
# ~/menghw_HD/DdCBE_project/10.hic_data/01.hicpro_WT/hicexp_result/01.TAD_info/293T_WT.KR.25Kbp_score.bigwig \
# ~/menghw_HD/DdCBE_project/10.hic_data/01.hicpro_WT/hicexp_result/01.TAD_info/293T_WT.KR_score.bigwig \
# -R \
# ~/menghw_HD/DdCBE_project/08.DdCBE_merge_all/region_cor_analysis/peak_region/20210312-Motif-ND6-DepSite.merge.bed \
# ~/menghw_HD/DdCBE_project/08.DdCBE_merge_all/region_cor_analysis/peak_region/20210312-Motif-ND6-Indep.merge.bed \
# ~/menghw_HD/DdCBE_project/08.DdCBE_merge_all/region_cor_analysis/peak_region/20210319-hg38_random.50bp.sort.bed \
# --referencePoint center \
# --beforeRegionStartLength 2000000 \
# --afterRegionStartLength 2000000 \
# --skipZeros \
# --binSize 25000 \
# -o deeptools_result/20210421-TAD_score.profile.mat.gz \
# --samplesLabel  TAD.25Kbp  TAD.5Kbp \
# --numberOfProcessors 10 & 

# # plot
# plotHeatmap -m deeptools_result/20210421-TAD_score.profile.mat.gz \
# --colorMap Purples \
# -out deeptools_result/20210421-TAD_score.profile.heatmap.pdf &

# # count site in TAD boundary
# Dep 3 / 72
# Indep 51 / 417 
# Random 16 / 500

# # count site in TAD boundary +- 1 pixel
# Dep 8 / 72
# Indep 134 / 417 
# Random 58 / 500
