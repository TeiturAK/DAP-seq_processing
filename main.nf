nextflow.enable.dsl=2

/*
  General DAP-seq pipeline

  Expected samplesheet columns:
  sample,role,fastq_1,fastq_2

  role must be either:
  - tf
  - input
*/

Channel
  .fromPath(params.samplesheet, checkIfExists: true)
  .splitCsv(header: true)
  .map { row ->
    tuple(
      row.sample as String,
      row.role.toLowerCase() as String,
      file(row.fastq_1),
      file(row.fastq_2)
    )
  }
  .set { samples }

process FASTQC_RAW {
  tag { sample }
  publishDir "${params.outdir}/fastqc_raw", mode: 'copy', pattern: '*.{html,zip}', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), val(role), path(read1), path(read2)

  output:
  path "*.html", emit: html
  path "*.zip", emit: zip
  path "${sample}.fastqc_raw.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.fastqc_sif} \
      fastqc \
        ${params.fastqc_args} \
        --outdir . \
        ${read1} ${read2}
  } 2> ${sample}.fastqc_raw.err
  """
}

process MULTIQC_RAW {
  tag 'raw_fastqc'
  publishDir "${params.outdir}/multiqc_raw", mode: 'copy', pattern: 'multiqc_*', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  path qc_files

  output:
  path "multiqc_report.html", emit: report
  path "multiqc_data", emit: data
  path "multiqc_raw.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.multiqc_sif} \
      multiqc \
        -o . \
        ${qc_files}
  } 2> multiqc_raw.err
  """
}

process TRIMMOMATIC_PE {
  tag { sample }
  publishDir "${params.outdir}/trimmomatic", mode: 'copy', pattern: '*.fq.gz', overwrite: true
  publishDir "${params.outdir}/trimmomatic", mode: 'copy', pattern: '*.trimmomatic.log', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), val(role), path(read1), path(read2)

  output:
  tuple val(sample), val(role), path("${sample}.trimmomatic_1.fq.gz"), path("${sample}.trimmomatic_2.fq.gz"), emit: reads
  path "${sample}.trimmomatic.log", emit: log
  path "${sample}.unpaired_1.fq.gz", emit: unpaired1
  path "${sample}.unpaired_2.fq.gz", emit: unpaired2
  path "${sample}.trimmomatic.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.trimmomatic_sif} \
      trimmomatic PE \
      -threads ${task.cpus} \
      -trimlog ${sample}.trimmomatic.log \
      ${read1} \
      ${read2} \
      ${sample}.trimmomatic_1.fq.gz \
      ${sample}.unpaired_1.fq.gz \
      ${sample}.trimmomatic_2.fq.gz \
      ${sample}.unpaired_2.fq.gz \
      ${params.trimmomatic_args}
  } 2> ${sample}.trimmomatic.err
  """
}

process FASTQC_TRIMMED {
  tag { sample }
  publishDir "${params.outdir}/fastqc_trimmed", mode: 'copy', pattern: '*.{html,zip}', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), val(role), path(read1), path(read2)

  output:
  path "*.html", emit: html
  path "*.zip", emit: zip
  path "${sample}.fastqc_trimmed.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.fastqc_sif} \
      fastqc \
        ${params.fastqc_args} \
        --outdir . \
        ${read1} ${read2}
  } 2> ${sample}.fastqc_trimmed.err
  """
}

process MULTIQC_TRIMMED {
  tag 'trimmed_fastqc'
  publishDir "${params.outdir}/multiqc_trimmed", mode: 'copy', pattern: 'multiqc_*', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  path qc_files

  output:
  path "multiqc_report.html", emit: report
  path "multiqc_data", emit: data
  path "multiqc_trimmed.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.multiqc_sif} \
      multiqc \
        -o . \
        ${qc_files}
  } 2> multiqc_trimmed.err
  """
}

process BOWTIE2_ALIGN_SORT {
  tag { sample }
  publishDir "${params.outdir}/bowtie2", mode: 'copy', pattern: '*.sorted.bam*', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), val(role), path(read1), path(read2)

  output:
  tuple val(sample), val(role), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.csi"), emit: bam
  path "${sample}.bowtie2_align_sort.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.bowtie2_sif} \
      bowtie2 \
      ${params.bowtie2_args} \
      -x ${params.bowtie2_index} \
      -1 ${read1} \
      -2 ${read2} \
      --threads ${task.cpus} | \
    singularity exec ${params.samtools_sif} \
      samtools sort \
        -@ ${task.cpus} \
        -O bam \
        -o ${sample}.sorted.bam \
        -

    singularity exec ${params.samtools_sif} \
      samtools index \
        -c \
        -@ ${task.cpus} \
        ${sample}.sorted.bam
  } 2> ${sample}.bowtie2_align_sort.err
  """
}

process SAMTOOLS_FILTER_Q20 {
  tag { sample }
  publishDir "${params.outdir}/bowtie2.filtered_q20", mode: 'copy', pattern: '*.filtered.sorted.bam*', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), val(role), path(bam), path(bai)

  output:
  tuple val(sample), val(role), path("${sample}.filtered.sorted.bam"), path("${sample}.filtered.sorted.bam.csi"), emit: bam
  path "${sample}.samtools_filter_q20.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.samtools_sif} \
      samtools view \
      -@ ${task.cpus} \
      ${params.samtools_filter_args} \
      ${bam} | \
    singularity exec ${params.samtools_sif} \
      samtools sort \
      -@ ${task.cpus} \
      -o ${sample}.filtered.sorted.bam \
      -

    singularity exec ${params.samtools_sif} \
      samtools index \
      -@ ${task.cpus} \
      -c \
      ${sample}.filtered.sorted.bam
  } 2> ${sample}.samtools_filter_q20.err
  """
}

process SAMTOOLS_MARKDUP_REMOVE {
  tag { sample }
  publishDir "${params.outdir}/bowtie2.filtered_q20.markdup", mode: 'copy', pattern: '*.filtered.markdup.sorted.bam*', overwrite: true
  publishDir "${params.outdir}/bowtie2.filtered_q20.markdup", mode: 'copy', pattern: '*.markdup.stats.txt', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), val(role), path(bam), path(bai)

  output:
  tuple val(sample), val(role), path("${sample}.filtered.markdup.sorted.bam"), path("${sample}.filtered.markdup.sorted.bam.csi"), emit: bam
  path "${sample}.markdup.stats.txt", emit: stats
  path "${sample}.samtools_markdup_remove.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.samtools_sif} \
      samtools collate \
      -@ ${task.cpus} \
      ${params.samtools_collate_args} \
      -o ${sample}.filtered.name-sorted.bam \
      ${bam}

    singularity exec ${params.samtools_sif} \
      samtools fixmate \
      -@ ${task.cpus} \
      ${params.samtools_fixmate_args} \
      ${sample}.filtered.name-sorted.bam \
      ${sample}.filtered.fixmate.bam

    singularity exec ${params.samtools_sif} \
      samtools sort \
      -@ ${task.cpus} \
      -o ${sample}.filtered.positionsorted.bam \
      ${sample}.filtered.fixmate.bam

    singularity exec ${params.samtools_sif} \
      samtools markdup \
      -@ ${task.cpus} \
      ${params.samtools_markdup_args} \
      ${sample}.filtered.positionsorted.bam \
      ${sample}.filtered.markdup.sorted.bam \
      2> >(tee ${sample}.markdup.stats.txt >&2)

    singularity exec ${params.samtools_sif} \
      samtools index \
      -@ ${task.cpus} \
      -c \
      ${sample}.filtered.markdup.sorted.bam
  } 2> ${sample}.samtools_markdup_remove.err
  """
}

process MACS3_CALLPEAK {
  tag { sample }
  publishDir "${params.outdir}/macs3", mode: 'copy', pattern: '*_peaks.*', overwrite: true
  publishDir "${params.outdir}/macs3", mode: 'copy', pattern: '*_summits.bed', overwrite: true
  publishDir "${params.outdir}/macs3", mode: 'copy', pattern: '*_model.r', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), path(treatment_bam), path(control_bam)

  output:
  tuple val(sample), path("${sample}_peaks.narrowPeak"), emit: narrowpeak
  path "${sample}_peaks.xls", emit: xls
  path "${sample}_summits.bed", emit: summits
  path "${sample}_model.r", optional: true, emit: model
  path "${sample}.macs3_callpeak.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.macs3_sif} \
      macs3 callpeak \
      -t ${treatment_bam} \
      -c ${control_bam} \
      ${params.macs3_args} \
      --outdir . \
      --name ${sample}
  } 2> ${sample}.macs3_callpeak.err
  """
}

process TOP_SIGNAL_PEAKS {
  tag { sample }
  publishDir "${params.outdir}/macs3.top${params.n_top_peaks}", mode: 'copy', pattern: '*.topSignal.*', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), path(narrowpeak)

  output:
  tuple val(sample), path("${sample}_peaks.topSignal.narrowPeak"), path("${sample}_peaks.topSignal.bed"), emit: top_peaks
  path "${sample}.top_signal_peaks.err", emit: err

  script:
  """
  set -euo pipefail

  {
    sort -t '	' -k7,7rn ${narrowpeak} | \
      head -n ${params.n_top_peaks} > ${sample}_peaks.tmp

    singularity exec ${params.bedtools_sif} \
      bedtools sort \
      -i ${sample}_peaks.tmp \
      > ${sample}_peaks.topSignal.narrowPeak

    awk -F '\t' -v OFS='\t' '{ midpoint = \$2 + int((\$3 - \$2) / 2); print \$1, midpoint, midpoint + 1 }' \
      ${sample}_peaks.topSignal.narrowPeak \
      > ${sample}_peaks.topSignal.bed
  } 2> ${sample}.top_signal_peaks.err
  """
}

process MEME_CHIP {
  tag { sample }
  publishDir "${params.outdir}/MEME", mode: 'copy', pattern: '*_peaks.topSignal*', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), path(narrowpeak), path(summit_bed)

  output:
  path "${sample}_peaks.topSignal", emit: meme_out
  path "${sample}_peaks.topSignal.summits.slop_${params.meme_slop_bp}bp.bed", emit: bed
  path "${sample}_peaks.topSignal.summits.slop_${params.meme_slop_bp}bp.fasta", emit: fasta
  path "${sample}.meme_chip.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.bedtools_sif} \
      bedtools slop \
        -b ${params.meme_slop_bp} \
        -i ${summit_bed} \
        -g ${params.genome_fai} \
        > ${sample}_peaks.topSignal.summits.slop_${params.meme_slop_bp}bp.bed

    singularity exec ${params.bedtools_sif} \
      bedtools getfasta \
        -fi ${params.genome_fasta} \
        -bed ${sample}_peaks.topSignal.summits.slop_${params.meme_slop_bp}bp.bed \
        > ${sample}_peaks.topSignal.summits.slop_${params.meme_slop_bp}bp.fasta

    singularity exec ${params.meme_sif} \
      meme-chip \
      -db ${params.motif_db} \
      -oc ${sample}_peaks.topSignal \
      ${params.meme_args} \
      ${sample}_peaks.topSignal.summits.slop_${params.meme_slop_bp}bp.fasta
  } 2> ${sample}.meme_chip.err
  """
}

workflow {
  FASTQC_RAW(samples)
  MULTIQC_RAW(FASTQC_RAW.out.zip.collect())

  TRIMMOMATIC_PE(samples)
  trimmed_reads = TRIMMOMATIC_PE.out.reads

  FASTQC_TRIMMED(trimmed_reads)
  MULTIQC_TRIMMED(FASTQC_TRIMMED.out.zip.collect())

  BOWTIE2_ALIGN_SORT(trimmed_reads)
  aligned_bams = BOWTIE2_ALIGN_SORT.out.bam

  SAMTOOLS_FILTER_Q20(aligned_bams)
  filtered_bams = SAMTOOLS_FILTER_Q20.out.bam

  SAMTOOLS_MARKDUP_REMOVE(filtered_bams)
  markdup_bams = SAMTOOLS_MARKDUP_REMOVE.out.bam

  control_bam = markdup_bams
    .filter { sample, role, bam, bai -> role == params.input_role }
    .map { sample, role, bam, bai -> bam }
    .first()

  treatment_bams = markdup_bams
    .filter { sample, role, bam, bai -> role == params.tf_role }
    .map { sample, role, bam, bai -> tuple(sample, bam) }

  macs3_input = treatment_bams.combine(control_bam)
    .map { sample, treatment_bam, control -> tuple(sample, treatment_bam, control) }

  MACS3_CALLPEAK(macs3_input)
  peaks = MACS3_CALLPEAK.out.narrowpeak

  TOP_SIGNAL_PEAKS(peaks)
  top_peaks = TOP_SIGNAL_PEAKS.out.top_peaks

  MEME_CHIP(top_peaks)
}
