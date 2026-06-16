# DAP-seq processing

Nextflow workflow for processing and analysis scripts comparing and benchmarking DAP-seq datasets aginst publically available DAP-seq 
data.

The repository contains a general DAP-seq processing workflow that can be run in either paired-end or single-end mode, together with dataset-specific configuration files and downstream R scripts for comparing sample-level statistics.

## Repository structure

```text
DAP-seq_processing/
├── analysis/        # R scripts and rendered reports for downstream comparisons
├── config/          # Dataset-specific Nextflow configuration files
├── pipeline/
│   ├── nf/          # Main Nextflow workflow
│   └── bash/        # Helper scripts
└── samplesheets/    # Example/input samplesheets
```


## Workflow overview

The Nextflow pipeline performs the main processing and QC steps for DAP-seq data:

1. Raw read QC with FastQC and MultiQC
2. Adapter and quality trimming with Trimmomatic
3. Trimmed read QC with FastQC and MultiQC
4. Alignment with Bowtie2
5. BAM filtering with SAMtools
6. Duplicate removal with SAMtools markdup
7. Peak calling with MACS3
8. FRiP calculation
9. Peak overlap with promoter annotations
10. Selection of top signal peaks
11. Motif analysis with MEME-ChIP
12. Collection of per-sample summary statistics

The workflow supports both paired-end and single-end data. The read layout is set in the config file using:

    params.read_layout = 'PE'

or

    params.read_layout = 'SE'

## Input samplesheet

The pipeline expects a CSV samplesheet with the following columns:

    sample,role,fastq_1,fastq_2

For paired-end data:

    DAL1_1,tf,/path/to/DAL1_1_R1.fastq.gz,/path/to/DAL1_1_R2.fastq.gz
    DAL1_Input,input,/path/to/input_R1.fastq.gz,/path/to/input_R2.fastq.gz

For single-end data, leave `fastq_2` empty:

    DAL1_1,tf,/path/to/DAL1_1.fastq.gz,
    DAL1_Input,input,/path/to/input.fastq.gz,

The `role` column should use the labels defined in the config file, normally:

    tf
    input

## Running the pipeline

Example run:

    nextflow run pipeline/nf/main.nf \
      -c config/nextflow_aspen.config \
      -profile slurm

## Configuration

Dataset-specific settings are kept in separate config files under `config/`.

These include:

- input samplesheet
- output directory
- reference genome files
- Bowtie2 index
- promoter annotation BED file
- motif database
- PE/SE read layout
- tool arguments
- container paths
- SLURM resources

Before running the workflow on a new system or dataset, update the relevant paths in the config file.

## Outputs

The main output directory is set by:

    params.outdir

The pipeline writes outputs into separate subdirectories, including:

    fastqc_raw/
    multiqc_raw/
    trimmomatic/
    fastqc_trimmed/
    multiqc_trimmed/
    bowtie2/
    bowtie2.filtered_q20/
    bowtie2.filtered_q20.markdup/
    macs3/
    frip/
    promoter_intersect/
    top_signal_peaks/
    meme_chip/
    sample_stats/
    pipeline_info/logs/

The most useful final outputs are usually:

- MACS3 peak files
- FRiP tables
- promoter overlap files
- MEME-ChIP output directories
- per-sample summary statistics
- MultiQC reports

## Downstream analysis

The `analysis/` directory contains R scripts and rendered HTML reports for comparing processed DAP-seq datasets and pipeline modes.

These scripts summarize metrics such as:

- input read counts
- read survival after trimming/filtering
- alignment rate
- duplicate rate
- FRiP
- number of peaks
- promoter overlaps
- motif detection

