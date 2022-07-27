# snakepipes_Hi-C
**须知**：本仓库还在构建中，暂时只作参考！！

---
## 环境
```shell
conda install python r-base bowtie2 samtools iced r-ggplot2 r-rcolorbrewer 

# 我用的版本
# python=3.9.13
# R=4.0.5
# bowtie2=2.4.5
# samtools=1.15.1
# iced=0.5.10
```

## 用法

### step 0 测序质量控制
使用 [snakepipes_fastqc-multiqc](https://github.com/hermanzhaozzzz/snakepipes_fastqc-multiqc)进行质量控制

### step 1 运行Snakemake Pipeline，生成Hi-C contact matrix
**回贴Hi-C reads以及生成RAW矩阵ICE校正矩阵**

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
# ├── matrix
# ├── qc
# ├── quality_checks
# ├── snakepipes_fastqc-multiqc
# ├── snakepipes_Hi-C
# ├── temp_files
# └── valid_pairs
```
![](https://tva1.sinaimg.cn/large/e6c9d24ely1h4lrqeisrlj20kv0cogng.jpg)

![](https://tva1.sinaimg.cn/large/e6c9d24ely1h4lrq238cij20kc0onjvm.jpg)

![](https://tva1.sinaimg.cn/large/e6c9d24ely1h4lrqxfl53j20kn0n1juq.jpg)

![](https://tva1.sinaimg.cn/large/e6c9d24ely1h4lrr7zxxyj20kc0my0wt.jpg)

### step 2