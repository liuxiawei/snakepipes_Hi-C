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
# pip install iced                           
# Collecting iced
#   Downloading iced-0.5.10.tar.gz (2.3 MB)
# manually set cmd path
BOWTIE2 = "bowtie2"
SAMTOOLS = "samtools"



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
        [[ "{input}" =~ .*SE.fastq.gz$ ]] &&
        # if true
        echo "[FATAL] find SE reads, raw reads should be PE reads in Hi-C protocol!" > {log} 2>&1 ||
        # if false
        echo "[DEBUG] find PE reads, go on mapping!" > {log} 2>&1
        {BOWTIE2} {BOWTIE2_GLOBAL_OPTIONS} --un {output.un} --rg-id BMG --rg {params.RG} -p {THREAD} -x {BOWTIE2_INDEX} -U {params.R1} -S {output.sam} >> {log} 2>&1
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
        [[ "{input}" =~ .*SE.fastq.gz$ ]] &&
        # if true
        echo "[FATAL] find SE reads, raw reads should be PE reads in Hi-C protocol!" > {log} 2>&1 ||
        # if false
        echo "[DEBUG] find PE reads, go on mapping!" > {log} 2>&1
        {BOWTIE2} {BOWTIE2_GLOBAL_OPTIONS} --un {output.un} --rg-id BMG --rg {params.RG} -p {THREAD} -x {BOWTIE2_INDEX} -U {params.R2} -S {output.sam} >> {log} 2>&1
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
        temp_path="../temp_files/merged_bam_splits_for_sortn/{sample}_R1"
    log:
        "../bam/{sample}/{sample}_R1.merged_sortn.log"
    shell:
        """
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
        temp_path="../temp_files/merged_bam_splits_for_sortn/{sample}_R2"
    log:
        "../bam/{sample}/{sample}_R2.merged_sortn.log"
    shell:
        """
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
        # [TODO] 这里到底是0-base还是1-base有待考证，先用0试试
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