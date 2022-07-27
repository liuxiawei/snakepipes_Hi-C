# snakepipes_Hi-C
**须知**：本仓库还在构建中，暂时只作参考！！

参考和引用了一些[HiC Pro](https://github.com/nservant/HiC-Pro)的代码

---
## 环境
```shell
conda install python r-base bowtie2 samtools iced r-ggplot2 r-rcolorbrewer
conda install -c bioconda java-jdk hicexplorer

# 我用的版本
# python=3.9.13
# R=4.0.5
# bowtie2=2.4.5
# samtools=1.15.1
# iced=0.5.10
# java-jdk=1.8 # java openjdk version "1.8.0_312"
# hicexplorer=3.7.2
```

## 用法

### step 0 测序质量控制
使用 [snakepipes_fastqc-multiqc](https://github.com/hermanzhaozzzz/snakepipes_fastqc-multiqc)进行质量控制

### step 1 运行Snakemake Pipeline，生成Hi-C contact matrix
- **回贴Hi-C reads以及生成RAW矩阵ICE校正矩阵**
- **validPairs convert to .hic file(Juicer)**

```shell
cd HiC
git clone https://github.com/hermanzhaozzzz/snakepipes_fastqc-multiqc.git
cd snakepipes_fastqc-multiqc

# then

# use jupyterlab or runipy to run step01_generate_samples.ipynb
# get samples.json and check it

# then

# dry run, rm -n to run pipeline
snakemake -pr -j 8 -s step02_run_mapping_and_generate_matrix.py -n

# output as below
# HiC|⇒ tree . -L 1
# .
# ├── bam
# ├── fastq
# ├── hic_file
# ├── matrix
# ├── qc
# ├── quality_checks
# ├── snakepipes_fastqc-multiqc
# ├── snakepipes_Hi-C
# ├── temp_files
# └── valid_pairs
```
### step 2 Convert ValidPairs to Juicer .hic¶