# biopython<=1.76 [Bio.Alphabet was used in hicexplorer] pip install biopython==1.76

ls_sample = ['dixon_2M/SRR400264_00', 'dixon_2M_2/SRR400264_01']
ls_bin = [5000, 10000, 20000, 40000, 100000, 150000, 500000, 1000000]
THREADS = 6



restrictionCutFile = "/lustre1/chengqiyi_pkuhpc/zhaohn/0.apps/HiC-Pro_installed/HiC-Pro_3.1.0/annotation/HindIII_resfrag_hg38.bed"
restrictionSequence = "AAGCTT" # HindIII
danglingSequence = "AGCT" # HindIII
chromosomeSizes = "/lustre1/chengqiyi_pkuhpc/zhaohn/0.apps/HiC-Pro_installed/HiC-Pro_3.1.0/annotation/chrom_hg38.sizes"
        
rule all:
    input:
        expand("../out_dir/hic_results/matrix/{sample}_bin-{bins}.h5", sample=ls_sample, bins=ls_bin)
rule bam2h5:
    input:
        "../out_dir/bowtie_results/bwt2/{sample}_R1_genome_ucsc_hg38.fa.bowtie2_index.bwt2merged.bam",
        "../out_dir/bowtie_results/bwt2/{sample}_R2_genome_ucsc_hg38.fa.bowtie2_index.bwt2merged.bam"
    output:
        h5="../out_dir/hic_results/matrix/{sample}_bin-{bins}.h5",
        qc="../out_dir/hic_results/matrix/{sample}_bin-{bins}.qc"
    params:
        bins = lambda wildcards, output: output[0].split("_bin-")[1].replace(".h5", "")
    shell:
        """
        hicBuildMatrix --samFiles {input[0]} {input[1]} \
            --outFileName {output.h5} \
            --QCfolder {output.qc} \
            --restrictionCutFile {restrictionCutFile} \
            --restrictionSequence {restrictionSequence} \
            --danglingSequence {danglingSequence} \
            --chromosomeSizes {chromosomeSizes} \
            --binSize {params.bins} \
            --threads {THREADS} \
            --inputBufferSize 400000 \
            # --doTestRun
            
        """

# rule h5_normalization:
#     input:
#         pass
#     output:
#         pass
#     shell:
#         "hicNormalize -m {h5} --normalize norm_range -o {h5nm}"
# rule h5_correction:
#     input:
#         pass
#     output:
#         pass
#     shell:
#         "hicCorrectMatrix correct -m {h5nm} --correctionMethod KR -o {h5kr}"