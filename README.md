# Nextflow for RNA-seq Training on AWS

This repository provides a cloud-based implementation of RNA-seq workflows for training and experimentation purposes. It combines content from the official [Nextflow for RNA-seq training course](https://training.nextflow.io/latest/nf4_science/rnaseq/) with real-world RNA-seq data from the [DIY Transcriptomics course](https://diytranscriptomics.com/), executed entirely on AWS infrastructure.

The goal is to demonstrate how to run scalable and reproducible RNA-seq pipelines using Nextflow, while leveraging AWS compute resources for performance, flexibility, and cost optimization.

Included in this repository:

- In the nextflow-for-rnaseq-aws folder, a reproduction of the Nextflow RNA-seq training pipeline, adapted for AWS execution.
- In the rnaseq-aws-diy-transcriptomics, an RNA-seq workflow using **HISAT2** to align reads from ~30 GB of skin tissue data (10 samples) and explore infrastructure setup, instance selection, and queue-based workload distribution in AWS.
- In the rnaseq-aws-diy-kallisto folder, a lightweight RNA-seq analysis using **Kallisto** as a pseudoaligner to quickly generate expression quantification (`abundance.tsv`) files for the same dataset.

This setup serves as both a training resource and a practical guide for running RNA-seq analyses in the cloud.

