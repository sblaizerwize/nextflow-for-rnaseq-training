# RNA-seq analysis of data from the DIY Transcriptomics course using Kallisto and AWS resources

This pipeline performs an RNA-seq analysis on skin data provided by the [DIY Transcriptomics course](https://diytranscriptomics.com/) (approximately 30 GB), which includes 10 patient samples. The goal is to use Kallisto as the pseudoaligner to generate `abundance.tsv` files for each sample. 

Here is a summary of the pipeline. 

```
nextflow run . -profile batch   --reads 's3://diy-transcriptomics-rna-seq-data/data/single-end.csv'   --transcriptome 's3://diy-transcriptomics-rna-seq-data/data/homo-genome.index'   --outdir 's3://nextflow-rnaseq-kallisto/results'   -with-report report.html   -with-timeline timeline.html   -with-trace trace.txt

 N E X T F L O W   ~  version 25.04.6

Launching `./main.nf` [exotic_monod] DSL2 - revision: b4a11516b4

        R N A S E Q - N F   P I P E L I N E
        ===================================
        transcriptome: s3://diy-transcriptomics-rna-seq-data/data/homo-genome.index
        reads        : s3://diy-transcriptomics-rna-seq-data/data/single-end.csv
        outdir       : s3://nextflow-rnaseq-kallisto/results

executor >  awsbatch (31)
[f7/7cfb36] FASTQC (10)        | 10 of 10 ✔
[fb/718ace] TRIM_GALORE (10)   | 10 of 10 ✔
[98/0b596c] KALLISTO_QUAN (10) | 10 of 10 ✔
[44/cba344] MULTIQC            | 1 of 1 ✔
Waiting for file transfers to complete (1 files)
Completed at: 30-Jul-2025 23:25:58
Duration    : 41m 52s
CPU hours   : 25.9
Succeeded   : 31
```