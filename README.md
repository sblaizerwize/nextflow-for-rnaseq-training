# Nextflow for RNA-seq Training on AWS

This repository provides a cloud-based implementation of RNA-seq workflows for training and experimentation purposes. It combines content from the official [Nextflow for RNA-seq training course](https://training.nextflow.io/latest/nf4_science/rnaseq/) with real-world RNA-seq data from the [DIY Transcriptomics course](https://diytranscriptomics.com/), executed entirely on AWS infrastructure.

The goal is to demonstrate how to run scalable and reproducible RNA-seq pipelines using Nextflow, while leveraging AWS compute resources for performance, flexibility, and cost optimization.

Visit the following link to consult the principal results after running the three presented pipelines: https://sblaizerwize.github.io/

We would like to thank AWS and the University of Navarra for providing cloud resources (1,000 AWS credits) to support this project.

## Table of Contents
- [Motivation](#introduction)
- [Prerequisites](#prerequisites)
- [Getting Started](#getting-started)
    - [Integrating AWS Batch with Nextflow](#integrating-aws-batch-with-nextflow)  
    - [Uploading RNA-seq data to an S3 bucket](#uploading-rna-seq-data-to-an-s3-bucket)
- [Running a Pipeline from this Repository](#running-a-pipeline-from-this-repository)  

---
## Motivation
Nextflow is an open-source tool used to run bioinformatics pipelines for analyzing multiple samples in parallel. As of 2025, Nextflow has established a solid community that actively contributes to its well-documented pipelines and templates for building custom workflows. It has become a popular and practical tool among bioinformaticians.

One of Nextflow’s key strengths is its portability. Whether you're running workflows on your laptop, a high-performance computing (HPC) cluster, or cloud platforms like AWS, Nextflow abstracts the complexity of the compute environment so you can focus on your analysis. 

This project began after we received $1,000 in AWS Cloud credits. Our goal was to lay the foundation for running a Bulk RNA-seq analysis using AWS infrastructure. To avoid discarding valuable insights that could benefit newcomers to Nextflow, we created this practical guide for running RNA-seq analyses on the AWS Cloud.

The content of this repository is structured as follows:

| **Directory** | **Description** |
| ------ | ------ |
| `nextflow-for-rnaseq-aws` | Includes a reproduction of the [Nextflow for RNA-seq](https://training.nextflow.io/latest/nf4_science/rnaseq/) pipeline, adapted for AWS execution. |
| `rnaseq-aws-diy-transcriptomics` | Implements an RNA-seq pipeline using skin tissue data (10 samples) from the [DIY Transcriptomics](https://diytranscriptomics.com/) course. It builds a real-world workflow using **HISAT2** for the alignment of reads, and explores AWS insfrastructure setup.  |
| `rnaseq-aws-diy-kallisto` | Adapts the `rnaseq-aws-diy-transcriptomics` workflow to use **Kallisto**, a pseudoaligner, for lightweight RNA-seq analysis. This pipeline generates expression quantification (`abundance.tsv`) files for the same dataset.  |

---
## Prerequisites
To run this repository's content on a macOS machine, follow these instructions:

- Install [Visual Studio Code](https://code.visualstudio.com/download)
- Install [Nextflow](https://www.nextflow.io/docs/latest/install.html)
- Install [Docker Desktop](https://docs.docker.com/desktop/setup/install/mac-install/).
- Dowload this repository [nextflow-for-rnaseq-training](https://github.com/sblaizerwize/nextflow-for-rnaseq-training/archive/refs/heads/main.zip) 

[⇧ back to top](#table-of-contents)

---
## Getting Started

### Integrating AWS Batch with Nextflow
Before running any pipeline from this repository, you must configure AWS to enable AWS Batch integration with Nextflow. We assume you already have an AWS account with administrative privileges.

To set this up, we strongly recommend following the instructions provided in LW Pembleton’s [Nextflow on AWS Batch](https://lpembleton.rbind.io/posts/nextflow-on-aws-batch/) blog post. It's comprehensive and aligns well with the structure of this project.

> **New to Nextflow?** We suggest completing the "Hello Nextflow" and "Nextflow for RNA-seq" tutorials before diving into this repository. It’ll save you time and confusion.


### Uploading RNA-seq data to an S3 bucket
Once AWS Batch is integrated with Nextflow, you need to upload the RNA-seq data to be analyzed to AWS:

1. Create a new S3 bucket named `nextflow-for-rnaseq`. Make sure the bucket region matches your AWS Batch infrastructure region. 

2. Download the RNA-seq data from the [nextflow-io/training](https://github.com/nextflow-io/training) repository, specifically from the [nf4-science/rnaseq/data/reads](https://github.com/nextflow-io/training/tree/master/nf4-science/rnaseq/data/reads) directory. 

3. Upload the downloaded data to the `nextflow-for-rnaseq` bucket using the AWS CLI or AWS Console. We recommend the AWS CLI if you've already configured credentials locally.

Your S3 bucket should have the following structure:
```
.
├── data
│   ├── aligned
│   │   └── genome_index.tar.gz
│   ├── reads
│   │   ├── ENCSR000COQ1_1.fastq.gz
│   │   ├── ENCSR000COQ1_2.fastq.gz
│   │   ├── ENCSR000COQ2_1.fastq.gz
│   │   ├── ENCSR000COQ2_2.fastq.gz
│   │   ├── ENCSR000COR1_1.fastq.gz
│   │   ├── ENCSR000COR1_2.fastq.gz
│   │   ├── ENCSR000COR2_1.fastq.gz
│   │   ├── ENCSR000COR2_2.fastq.gz
│   │   ├── ENCSR000CPO1_1.fastq.gz
│   │   ├── ENCSR000CPO1_2.fastq.gz
│   │   ├── ENCSR000CPO2_1.fastq.gz
│   │   └── ENCSR000CPO2_2.fastq.gz
│   └── single-end.csv
├── results
└── test_env
```

Here is a summary of the S3 bucket content:

| **Directory** | **Description** |
|---------------|-----------------|
| `data`        | Contains input files for the pipeline, including RNA-seq FASTQ files (raw reads) and the reference genome index (e.g., for HISAT2 or Kallisto). |
| `results`     | Output directory where processed results are stored, such as alignment files (e.g., BAM), quantification results (e.g., counts or TPMs), and quality reports (e.g., FastQC, MultiQC). |
| `test_env`    | Working directory used by Nextflow to store intermediate files, temporary data, and execution logs for each process run. This may include `.command.log`, `.command.err`, and subdirectories for task-specific runs. |

```
> ⚠️ **Note:** The `single-end.csv` file contains references to the single-end samples located in the `nextflow-for-rnaseq` bucket. This approach simplifies pipeline execution.
```

Here is the content of the CSV file:
```
sample_id,fastq_path
ENCSR000COQ1,s3://nextflow-for-rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,s3://nextflow-for-rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,s3://nextflow-for-rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,s3://nextflow-for-rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,s3://nextflow-for-rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,s3://nextflow-for-rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

---
## Running a Pipeline from this Repository  
Once you've set up AWS Batch and your S3 bucket, you're ready to run the `nextflow-for-rnaseq-aws` pipeline. Keep in mind that the provided instructions can be adapted for the other two pipelines. 

1. Open the repository in Visual Studio Code.  
2. Navigate into the `nextflow-for-rnaseq-aws` folder.  
3. Open the `nextflow.config` file.  
4. In the `profiles` section, locate the `batch` profile and update the following variables to match your AWS setup:

   ```groovy
   process.queue = 'nf-queue'
   workDir = 's3://nextflow-for-rnaseq/test_env'
   aws.region = 'eu-west-1'
   ```

5. Open a terminal in VS Code and run:
    ```
    nextflow run . -profile batch \
    --reads 's3://nextflow-for-rnaseq/data/single-end.csv' \
    --transcriptome 's3://nextflow-for-rnaseq/data/aligned/genome_index.tar.gz' \
    --outdir 's3://nextflow-for-rnaseq/results'
    ```

    > **Important:** Ensure your AWS credentials are correctly configured before executing the pipeline. You can verify your credentials using `aws configure list`:

    ```
    aws configure list
      Name                    Value             Type    Location
      ----                    -----             ----    --------
    profile                <not set>             None    None
    access_key     ****************HMMF shared-credentials-file
    secret_key     ****************6OfY shared-credentials-file
    region                eu-west-1      config-file    ~/.aws/config
    ```

6. You will get the following response:

    ```
     N E X T F L O W   ~  version 25.04.6

    Launching `./main.nf` [stoic_kare] DSL2 - revision: 58bf98848c

            R N A S E Q - N F   P I P E L I N E
            ===================================
            transcriptome: s3://nextflow-for-rnaseq/data/aligned/genome_index.tar.gz
            reads        : s3://nextflow-for-rnaseq/data/single-end.csv
            outdir       : s3://nextflow-for-rnaseq/results
    executor >  awsbatch (19)
    [f4/558877] FASTQC (4)       [100%] 6 of 6 ✔
    [11/11a646] TRIM_GALORE (5)  [100%] 6 of 6 ✔
    [30/a4f38c] HISAT2_ALIGN (6) [100%] 6 of 6 ✔
    [f8/8513cf] MULTIQC          [100%] 1 of 1 ✔
    Completed at: 18-Jul-2025 12:25:54
    Duration    : 3m 31s
    CPU hours   : 0.1
    Succeeded   : 19
    ```

And that's it! You can find more details on how to run each pipeline in its corresponding folder.


[⇧ back to top](#table-of-contents)


