nextflow.enable.dsl=2

/*
  General DAP-seq pipeline

  Expected samplesheet columns:
  sample,role,fastq_1,fastq_2

  Paired-end example:
  DAL1_1,tf,/path/R1.fastq.gz,/path/R2.fastq.gz
  DAL1_Input,input,/path/R1.fastq.gz,/path/R2.fastq.gz

  Single-end example:
  DAL1_1,tf,/path/read.fastq.gz,
  DAL1_Input,input,/path/input.fastq.gz,

  Select layout per run with:
  --read_layout PE
  or
  --read_layout SE

  role must be either:
  - tf
  - input
*/

def readLayout() {
  (params.read_layout ?: 'PE').toString().toUpperCase()
}

process FASTQC_RAW_PE {
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
    singularity exec ${params.singularity_bind} ${params.fastqc_sif} \
      fastqc \
        ${params.fastqc_args} \
        --outdir . \
        ${read1} ${read2}
  } 2> ${sample}.fastqc_raw.err
  """
}


process FASTQC_RAW_SE {
  tag { sample }
  publishDir "${params.outdir}/fastqc_raw", mode: 'copy', pattern: '*.{html,zip}', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), val(role), path(read1)

  output:
  path "*.html", emit: html
  path "*.zip", emit: zip
  path "${sample}.fastqc_raw.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.singularity_bind} ${params.fastqc_sif} \
      fastqc \
        ${params.fastqc_args} \
        --outdir . \
        ${read1}
  } 2> ${sample}.fastqc_raw.err
  """
}

process MULTIQC_RAW_PE {
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
    singularity exec ${params.singularity_bind} ${params.multiqc_sif} \
      multiqc \
        -o . \
        ${qc_files}
  } 2> multiqc_raw.err
  """
}


process MULTIQC_RAW_SE {
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
    singularity exec ${params.singularity_bind} ${params.multiqc_sif} \
      multiqc \
        -o . \
        ${qc_files}
  } 2> multiqc_raw.err
  """
}

process TRIMMOMATIC_PE {
  tag { sample }
  publishDir "${params.outdir}/trimmomatic", mode: 'symlink', pattern: '*.fq.gz', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), val(role), path(read1), path(read2)

  output:
  tuple val(sample), val(role), path("${sample}.trimmomatic_1.fq.gz"), path("${sample}.trimmomatic_2.fq.gz"), emit: reads
  path "${sample}.unpaired_1.fq.gz", emit: unpaired1
  path "${sample}.unpaired_2.fq.gz", emit: unpaired2
  path "${sample}.trimmomatic.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.singularity_bind} ${params.trimmomatic_sif} \
      trimmomatic PE \
      -threads ${task.cpus} \
      ${read1} \
      ${read2} \
      ${sample}.trimmomatic_1.fq.gz \
      ${sample}.unpaired_1.fq.gz \
      ${sample}.trimmomatic_2.fq.gz \
      ${sample}.unpaired_2.fq.gz \
      ${params.trimmomatic_pe_args}
  } 2> ${sample}.trimmomatic.err
  """
}


process TRIMMOMATIC_SE {
  tag { sample }
  publishDir "${params.outdir}/trimmomatic", mode: 'symlink', pattern: '*.fq.gz', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), val(role), path(read1)

  output:
  tuple val(sample), val(role), path("${sample}.trimmomatic.fq.gz"), emit: reads
  path "${sample}.trimmomatic.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.singularity_bind} ${params.trimmomatic_sif} \
      trimmomatic SE \
      -threads ${task.cpus} \
      ${read1} \
      ${sample}.trimmomatic.fq.gz \
      ${params.trimmomatic_se_args}
  } 2> ${sample}.trimmomatic.err
  """
}

process FASTQC_TRIMMED_PE {
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
    singularity exec ${params.singularity_bind} ${params.fastqc_sif} \
      fastqc \
        ${params.fastqc_args} \
        --outdir . \
        ${read1} ${read2}
  } 2> ${sample}.fastqc_trimmed.err
  """
}


process FASTQC_TRIMMED_SE {
  tag { sample }
  publishDir "${params.outdir}/fastqc_trimmed", mode: 'copy', pattern: '*.{html,zip}', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), val(role), path(read1)

  output:
  path "*.html", emit: html
  path "*.zip", emit: zip
  path "${sample}.fastqc_trimmed.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.singularity_bind} ${params.fastqc_sif} \
      fastqc \
        ${params.fastqc_args} \
        --outdir . \
        ${read1}
  } 2> ${sample}.fastqc_trimmed.err
  """
}

process MULTIQC_TRIMMED_PE {
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
    singularity exec ${params.singularity_bind} ${params.multiqc_sif} \
      multiqc \
        -o . \
        ${qc_files}
  } 2> multiqc_trimmed.err
  """
}


process MULTIQC_TRIMMED_SE {
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
    singularity exec ${params.singularity_bind} ${params.multiqc_sif} \
      multiqc \
        -o . \
        ${qc_files}
  } 2> multiqc_trimmed.err
  """
}

process BOWTIE2_ALIGN_SORT_PE {
  tag { sample }
  publishDir "${params.outdir}/bowtie2", mode: 'symlink', pattern: '*.sorted.bam*', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), val(role), path(read1), path(read2)

  output:
  tuple val(sample), val(role), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.*"), emit: bam
  path "${sample}.bowtie2_align_sort.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.singularity_bind} ${params.bowtie2_sif} \
      bowtie2 \
      ${params.bowtie2_pe_args} \
      -x ${params.bowtie2_index} \
      -1 ${read1} \
      -2 ${read2} \
      --threads ${task.cpus} | \
    singularity exec ${params.singularity_bind} ${params.samtools_sif} \
      samtools sort \
        -@ ${task.cpus} \
        -O bam \
        -o ${sample}.sorted.bam \
        -

    singularity exec ${params.singularity_bind} ${params.samtools_sif} \
      samtools index \
        ${params.samtools_index_args} \
        -@ ${task.cpus} \
        ${sample}.sorted.bam
  } 2> ${sample}.bowtie2_align_sort.err
  """
}


process BOWTIE2_ALIGN_SORT_SE {
  tag { sample }
  publishDir "${params.outdir}/bowtie2", mode: 'symlink', pattern: '*.sorted.bam*', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), val(role), path(read1)

  output:
  tuple val(sample), val(role), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.*"), emit: bam
  path "${sample}.bowtie2_align_sort.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.singularity_bind} ${params.bowtie2_sif} \
      bowtie2 \
      ${params.bowtie2_se_args} \
      -x ${params.bowtie2_index} \
      -U ${read1} \
      --threads ${task.cpus} | \
    singularity exec ${params.singularity_bind} ${params.samtools_sif} \
      samtools sort \
        -@ ${task.cpus} \
        -O bam \
        -o ${sample}.sorted.bam \
        -

    singularity exec ${params.singularity_bind} ${params.samtools_sif} \
      samtools index \
        ${params.samtools_index_args} \
        -@ ${task.cpus} \
        ${sample}.sorted.bam
  } 2> ${sample}.bowtie2_align_sort.err
  """
}

process SAMTOOLS_FILTER_Q20_PE {
  tag { sample }
  publishDir "${params.outdir}/bowtie2.filtered_q20", mode: 'symlink', pattern: '*.filtered.sorted.bam*', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), val(role), path(bam), path(bai)

  output:
  tuple val(sample), val(role), path("${sample}.filtered.sorted.bam"), path("${sample}.filtered.sorted.bam.*"), emit: bam
  path "${sample}.samtools_filter_q20.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.singularity_bind} ${params.samtools_sif} \
      samtools view \
      -@ ${task.cpus} \
      ${params.samtools_filter_pe_args} \
      ${bam} | \
    singularity exec ${params.singularity_bind} ${params.samtools_sif} \
      samtools sort \
      -@ ${task.cpus} \
      -o ${sample}.filtered.sorted.bam \
      -

    singularity exec ${params.singularity_bind} ${params.samtools_sif} \
      samtools index \
      -@ ${task.cpus} \
      ${params.samtools_index_args} \
      ${sample}.filtered.sorted.bam
  } 2> ${sample}.samtools_filter_q20.err
  """
}


process SAMTOOLS_FILTER_Q20_SE {
  tag { sample }
  publishDir "${params.outdir}/bowtie2.filtered_q20", mode: 'symlink', pattern: '*.filtered.sorted.bam*', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), val(role), path(bam), path(bai)

  output:
  tuple val(sample), val(role), path("${sample}.filtered.sorted.bam"), path("${sample}.filtered.sorted.bam.*"), emit: bam
  path "${sample}.samtools_filter_q20.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.singularity_bind} ${params.samtools_sif} \
      samtools view \
      -@ ${task.cpus} \
      ${params.samtools_filter_se_args} \
      ${bam} | \
    singularity exec ${params.singularity_bind} ${params.samtools_sif} \
      samtools sort \
      -@ ${task.cpus} \
      -o ${sample}.filtered.sorted.bam \
      -

    singularity exec ${params.singularity_bind} ${params.samtools_sif} \
      samtools index \
      -@ ${task.cpus} \
      ${params.samtools_index_args} \
      ${sample}.filtered.sorted.bam
  } 2> ${sample}.samtools_filter_q20.err
  """
}

process SAMTOOLS_MARKDUP_REMOVE_PE {
  tag { sample }
  publishDir "${params.outdir}/bowtie2.filtered_q20.markdup", mode: 'symlink', pattern: '*.filtered.markdup.sorted.bam*', overwrite: true
  publishDir "${params.outdir}/bowtie2.filtered_q20.markdup", mode: 'copy', pattern: '*.markdup.stats.txt', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), val(role), path(bam), path(bai)

  output:
  tuple val(sample), val(role), path("${sample}.filtered.markdup.sorted.bam"), path("${sample}.filtered.markdup.sorted.bam.*"), emit: bam
  path "${sample}.markdup.stats.txt", emit: stats
  path "${sample}.samtools_markdup_remove.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.singularity_bind} ${params.samtools_sif} \
      samtools collate \
      -@ ${task.cpus} \
      ${params.samtools_collate_args} \
      -o ${sample}.filtered.name-sorted.bam \
      ${bam}

    singularity exec ${params.singularity_bind} ${params.samtools_sif} \
      samtools fixmate \
      -@ ${task.cpus} \
      ${params.samtools_fixmate_args} \
      ${sample}.filtered.name-sorted.bam \
      ${sample}.filtered.fixmate.bam

    singularity exec ${params.singularity_bind} ${params.samtools_sif} \
      samtools sort \
      -@ ${task.cpus} \
      -o ${sample}.filtered.positionsorted.bam \
      ${sample}.filtered.fixmate.bam

    singularity exec ${params.singularity_bind} ${params.samtools_sif} \
      samtools markdup \
      -@ ${task.cpus} \
      ${params.samtools_markdup_args} \
      ${sample}.filtered.positionsorted.bam \
      ${sample}.filtered.markdup.sorted.bam \
      2> >(tee ${sample}.markdup.stats.txt >&2)

    singularity exec ${params.singularity_bind} ${params.samtools_sif} \
      samtools index \
      -@ ${task.cpus} \
      ${params.samtools_index_args} \
      ${sample}.filtered.markdup.sorted.bam
  } 2> ${sample}.samtools_markdup_remove.err
  """
}

process SAMTOOLS_MARKDUP_REMOVE_SE {
  tag { sample }
  publishDir "${params.outdir}/bowtie2.filtered_q20.markdup", mode: 'symlink', pattern: '*.filtered.markdup.sorted.bam*', overwrite: true
  publishDir "${params.outdir}/bowtie2.filtered_q20.markdup", mode: 'copy', pattern: '*.markdup.stats.txt', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), val(role), path(bam), path(bai)

  output:
  tuple val(sample), val(role), path("${sample}.filtered.markdup.sorted.bam"), path("${sample}.filtered.markdup.sorted.bam.*"), emit: bam
  path "${sample}.markdup.stats.txt", emit: stats
  path "${sample}.samtools_markdup_remove.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.singularity_bind} ${params.samtools_sif} \
      samtools markdup \
      -@ ${task.cpus} \
      ${params.samtools_markdup_args} \
      ${bam} \
      ${sample}.filtered.markdup.sorted.bam \
      2> >(tee ${sample}.markdup.stats.txt >&2)

    singularity exec ${params.singularity_bind} ${params.samtools_sif} \
      samtools index \
      -@ ${task.cpus} \
      ${params.samtools_index_args} \
      ${sample}.filtered.markdup.sorted.bam
  } 2> ${sample}.samtools_markdup_remove.err
  """
}

process MACS3_CALLPEAK {
  tag { sample }
  publishDir "${params.outdir}/macs3", mode: 'copy', pattern: '*_peaks.*', overwrite: true
  publishDir "${params.outdir}/macs3", mode: 'copy', pattern: '*_summits.bed', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), path(treatment_bam)
  path(control_bams)

  output:
  tuple val(sample), path("${sample}_peaks.narrowPeak"), emit: narrowpeak
  path "${sample}_summits.bed", emit: summits
  path "${sample}.macs3_callpeak.err", emit: err

  script:
  def control_args = control_bams.join(' ')

  """
  set -euo pipefail

  {
    singularity exec ${params.singularity_bind} ${params.macs3_sif} \
      macs3 callpeak \
        -t ${treatment_bam} \
        -c ${control_args} \
        ${readLayout() == 'SE' ? params.macs3_se_args : params.macs3_pe_args} \
        --outdir . \
        --name ${sample}
  } 2> ${sample}.macs3_callpeak.err
  """
}

process FRIP_CALCULATION {
  tag { sample }

  publishDir "${params.outdir}/frip", mode: 'copy', pattern: '*.FRiP.tsv', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), path(bam), path(narrowpeak)

  output:
  tuple val(sample), path("${sample}.FRiP.tsv"), emit: frip
  path "${sample}.frip.err", emit: err

  script:
  """
  set -euo pipefail

  {
    total_reads=\$(singularity exec ${params.singularity_bind} ${params.samtools_sif} \
      samtools view -c ${bam})

    reads_in_peaks=\$(singularity exec ${params.singularity_bind} ${params.bedtools_sif} \
      bedtools sort -i ${narrowpeak} \
      | singularity exec ${params.singularity_bind} ${params.bedtools_sif} \
        bedtools merge -i stdin \
      | singularity exec ${params.singularity_bind} ${params.bedtools_sif} \
        bedtools intersect \
          -u \
          -nonamecheck \
          -a ${bam} \
          -b stdin \
          -ubam \
      | singularity exec ${params.singularity_bind} ${params.samtools_sif} \
        samtools view -c)

    frip=\$(awk -v reads_in_peaks="\${reads_in_peaks}" -v total_reads="\${total_reads}" \
      'BEGIN { if (total_reads > 0) print reads_in_peaks / total_reads; else print "NA" }')

    printf "sample\\tFRiP\\treads_in_peaks\\ttotal_reads\\n" > ${sample}.FRiP.tsv
    printf "%s\\t%s\\t%s\\t%s\\n" "${sample}" "\${frip}" "\${reads_in_peaks}" "\${total_reads}" >> ${sample}.FRiP.tsv
  } 2> ${sample}.frip.err
  """
}

process PEAK_PROMOTER_INTERSECT {
  tag { sample }

  publishDir "${params.outdir}/peak_promoter_intersect", mode: 'copy', pattern: '*.bed', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), path(narrowpeak)

  output:
  tuple val(sample), path("${sample}_peaks.promoter_intersect.bed"), emit: intersect
  path "${sample}.peak_promoter_intersect.err", emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.singularity_bind} ${params.bedtools_sif} \
      bedtools intersect \
        ${params.bedtools_intersect_args} \
        -a ${narrowpeak} \
        -b ${params.promoter_bed} \
        > ${sample}_peaks.promoter_intersect.bed
  } 2> ${sample}.peak_promoter_intersect.err
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
    sort -k7,7rn ${narrowpeak} > ${sample}_peaks.sorted_by_signal.tmp
    head -n ${params.n_top_peaks} ${sample}_peaks.sorted_by_signal.tmp > ${sample}_peaks.tmp
    
    singularity exec ${params.singularity_bind} ${params.bedtools_sif} \
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
  tuple val(sample), path("${sample}_peaks.topSignal"), emit: meme_out
  tuple val(sample), path("${sample}_peaks.topSignal.summits.slop_${params.meme_slop_bp}bp.bed"), emit: bed
  tuple val(sample), path("${sample}_peaks.topSignal.summits.slop_${params.meme_slop_bp}bp.fasta"), emit: fasta
  tuple val(sample), path("${sample}.meme_chip.err"), emit: err

  script:
  """
  set -euo pipefail

  {
    singularity exec ${params.singularity_bind} ${params.bedtools_sif} \
      bedtools slop \
        -b ${params.meme_slop_bp} \
        -i ${summit_bed} \
        -g ${params.genome_fai} \
        > ${sample}_peaks.topSignal.summits.slop_${params.meme_slop_bp}bp.bed

    singularity exec ${params.singularity_bind} ${params.bedtools_sif} \
      bedtools getfasta \
        -fi ${params.genome_fasta} \
        -bed ${sample}_peaks.topSignal.summits.slop_${params.meme_slop_bp}bp.bed \
        > ${sample}_peaks.topSignal.summits.slop_${params.meme_slop_bp}bp.fasta

    singularity exec ${params.singularity_bind} ${params.meme_sif} \
      meme-chip \
      -db ${params.motif_db} \
      -oc ${sample}_peaks.topSignal \
      ${params.meme_args} \
      ${sample}_peaks.topSignal.summits.slop_${params.meme_slop_bp}bp.fasta
  } 2> ${sample}.meme_chip.err
  """
}

process COLLECT_SAMPLE_STATS {
  tag { sample }

  publishDir "${params.outdir}/sample_stats", mode: 'copy', pattern: '*.stats.tsv', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), path(frip_tsv), path(promoter_intersect), path(narrowpeak), path(meme_dir)

  output:
  tuple val(sample), path("${sample}.stats.tsv"), emit: stats
  path "${sample}.collect_sample_stats.err", emit: err

  script:
  """
  set -euo pipefail

  {
    layout="${params.read_layout}"

    trim_log="${params.outdir}/pipeline_info/logs/${sample}.trimmomatic.err"
    bowtie_log="${params.outdir}/pipeline_info/logs/${sample}.bowtie2_align_sort.err"
    macs3_log="${params.outdir}/pipeline_info/logs/${sample}.macs3_callpeak.err"
    markdup_stats="${params.outdir}/bowtie2.filtered_q20.markdup/${sample}.markdup.stats.txt"

    meme_summary_tsv="${meme_dir}/summary.tsv"
    centrimo_tsv="${meme_dir}/centrimo_out/centrimo.tsv"

    out="${sample}.stats.tsv"

    printf "sample\\tmetric\\tvalue\\n" > "\$out"

    add_stat() {
      metric="\$1"
      value="\$2"
      value="\${value:-NA}"
      printf "%s\\t%s\\t%s\\n" "${sample}" "\$metric" "\$value" >> "\$out"
    }

    # Trimmomatic
    if [[ "\$layout" == "PE" ]]; then
      trim_line=""
      if [[ -s "\$trim_log" ]]; then
        trim_line=\$(grep "Input Read Pairs:" "\$trim_log" || true)
      fi

      trim_input_pairs=\$(echo "\$trim_line" | sed -n 's/.*Input Read Pairs: \\([0-9]*\\).*/\\1/p')
      trim_both_surviving=\$(echo "\$trim_line" | sed -n 's/.*Both Surviving: \\([0-9]*\\).*/\\1/p')

      add_stat "trimmomatic_input_pairs" "\$trim_input_pairs"
      add_stat "trimmomatic_both_surviving" "\$trim_both_surviving"

    elif [[ "\$layout" == "SE" ]]; then
      trim_line=""
      if [[ -s "\$trim_log" ]]; then
        trim_line=\$(grep "Input Reads:" "\$trim_log" || true)
      fi

      trim_input_reads=\$(echo "\$trim_line" | sed -n 's/.*Input Reads: \\([0-9]*\\).*/\\1/p')
      trim_surviving_reads=\$(echo "\$trim_line" | sed -n 's/.*Surviving: \\([0-9]*\\).*/\\1/p')

      add_stat "trimmomatic_input_reads" "\$trim_input_reads"
      add_stat "trimmomatic_surviving_reads" "\$trim_surviving_reads"
    fi

    # Bowtie2
    bowtie_overall_alignment_rate="NA"

    if [[ "\$layout" == "PE" ]]; then
      bowtie_conc_0="NA"
      bowtie_conc_0_pct="NA"
      bowtie_conc_1="NA"
      bowtie_conc_1_pct="NA"
      bowtie_conc_multi="NA"
      bowtie_conc_multi_pct="NA"

      if [[ -s "\$bowtie_log" ]]; then
        bowtie_conc_0=\$(awk '/aligned concordantly 0 times/ && !seen {print \$1; seen=1}' "\$bowtie_log" || true)
        bowtie_conc_0_pct=\$(awk '/aligned concordantly 0 times/ && !seen {gsub(/[()%]/, "", \$2); print \$2; seen=1}' "\$bowtie_log" || true)

        bowtie_conc_1=\$(awk '/aligned concordantly exactly 1 time/ {print \$1}' "\$bowtie_log" || true)
        bowtie_conc_1_pct=\$(awk '/aligned concordantly exactly 1 time/ {gsub(/[()%]/, "", \$2); print \$2}' "\$bowtie_log" || true)

        bowtie_conc_multi=\$(awk '/aligned concordantly >1 times/ {print \$1}' "\$bowtie_log" || true)
        bowtie_conc_multi_pct=\$(awk '/aligned concordantly >1 times/ {gsub(/[()%]/, "", \$2); print \$2}' "\$bowtie_log" || true)

        bowtie_overall_alignment_rate=\$(awk '/overall alignment rate/ {gsub(/%/, "", \$1); print \$1}' "\$bowtie_log" || true)
      fi

      add_stat "bowtie2_concordant_0" "\$bowtie_conc_0"
      add_stat "bowtie2_concordant_0_pct" "\$bowtie_conc_0_pct"
      add_stat "bowtie2_concordant_1" "\$bowtie_conc_1"
      add_stat "bowtie2_concordant_1_pct" "\$bowtie_conc_1_pct"
      add_stat "bowtie2_concordant_multi" "\$bowtie_conc_multi"
      add_stat "bowtie2_concordant_multi_pct" "\$bowtie_conc_multi_pct"

    elif [[ "\$layout" == "SE" ]]; then
      bowtie_unaligned="NA"
      bowtie_unaligned_pct="NA"
      bowtie_aligned_1="NA"
      bowtie_aligned_1_pct="NA"
      bowtie_aligned_multi="NA"
      bowtie_aligned_multi_pct="NA"

      if [[ -s "\$bowtie_log" ]]; then
        bowtie_unaligned=\$(awk '/aligned 0 times/ {print \$1; exit}' "\$bowtie_log" || true)
        bowtie_unaligned_pct=\$(awk '/aligned 0 times/ {gsub(/[()%]/, "", \$2); print \$2; exit}' "\$bowtie_log" || true)

        bowtie_aligned_1=\$(awk '/aligned exactly 1 time/ {print \$1; exit}' "\$bowtie_log" || true)
        bowtie_aligned_1_pct=\$(awk '/aligned exactly 1 time/ {gsub(/[()%]/, "", \$2); print \$2; exit}' "\$bowtie_log" || true)

        bowtie_aligned_multi=\$(awk '/aligned >1 times/ {print \$1; exit}' "\$bowtie_log" || true)
        bowtie_aligned_multi_pct=\$(awk '/aligned >1 times/ {gsub(/[()%]/, "", \$2); print \$2; exit}' "\$bowtie_log" || true)

        bowtie_overall_alignment_rate=\$(awk '/overall alignment rate/ {gsub(/%/, "", \$1); print \$1}' "\$bowtie_log" || true)
      fi

      add_stat "bowtie2_aligned_0" "\$bowtie_unaligned"
      add_stat "bowtie2_aligned_0_pct" "\$bowtie_unaligned_pct"
      add_stat "bowtie2_aligned_1" "\$bowtie_aligned_1"
      add_stat "bowtie2_aligned_1_pct" "\$bowtie_aligned_1_pct"
      add_stat "bowtie2_aligned_multi" "\$bowtie_aligned_multi"
      add_stat "bowtie2_aligned_multi_pct" "\$bowtie_aligned_multi_pct"
    fi

    add_stat "bowtie2_overall_alignment_rate" "\$bowtie_overall_alignment_rate"

    # Samtools markdup
    markdup_examined="NA"
    markdup_duplicate_total="NA"

    if [[ -s "\$markdup_stats" ]]; then
      markdup_examined=\$(awk -F ': ' '\$1 == "EXAMINED" {print \$2}' "\$markdup_stats" || true)
      markdup_duplicate_total=\$(awk -F ': ' '\$1 == "DUPLICATE TOTAL" {print \$2}' "\$markdup_stats" || true)
    fi

    add_stat "samtools_markdup_examined" "\$markdup_examined"
    add_stat "samtools_markdup_duplicate_total" "\$markdup_duplicate_total"

    # MACS3
    macs3_mean_fragment_size_treatment="NA"

    if [[ -s "\$macs3_log" ]]; then
      macs3_mean_fragment_size_treatment=\$(grep "mean fragment size is determined" "\$macs3_log" | sed -n 's/.*as \\([0-9.]*\\) bp from treatment.*/\\1/p' || true)
    fi

    add_stat "macs3_mean_fragment_size_treatment" "\$macs3_mean_fragment_size_treatment"

    # Intersect count
    peak_count="NA"
    promoter_intersect_count="NA"

    if [[ -s "${narrowpeak}" ]]; then
      peak_count=\$(wc -l < "${narrowpeak}")
    fi

    if [[ -s "${promoter_intersect}" ]]; then
      promoter_intersect_count=\$(wc -l < "${promoter_intersect}")
    fi

    add_stat "macs3_peak_count" "\$peak_count"
    add_stat "promoter_intersect_count" "\$promoter_intersect_count"

    # FRiP
    frip="NA"
    frip_reads_in_peaks="NA"
    frip_total_reads="NA"

    if [[ -s "${frip_tsv}" ]]; then
      frip=\$(awk 'NR == 2 {print \$2}' "${frip_tsv}" || true)
      frip_reads_in_peaks=\$(awk 'NR == 2 {print \$3}' "${frip_tsv}" || true)
      frip_total_reads=\$(awk 'NR == 2 {print \$4}' "${frip_tsv}" || true)
    fi

    add_stat "FRiP" "\$frip"
    add_stat "frip_reads_in_peaks" "\$frip_reads_in_peaks"
    add_stat "frip_total_reads" "\$frip_total_reads"

    # MEME-ChIP summary top motif
    meme_has_motif="no"
    meme_top_evalue="NA"

    if [[ -s "\$meme_summary_tsv" ]]; then
      meme_top_evalue=\$(awk -F '\\t' 'NR == 2 {print \$8}' "\$meme_summary_tsv" || true)
      meme_top_evalue=\${meme_top_evalue:-NA}

      if [[ "\$meme_top_evalue" != "NA" ]]; then
        meme_has_motif="yes"
      fi
    fi

    add_stat "meme_has_motif" "\$meme_has_motif"
    add_stat "meme_top_evalue" "\$meme_top_evalue"

    # CentriMo central motif enrichment
    centrimo_has_central_motif="no"
    centrimo_top_evalue="NA"
    centrimo_top_adj_pvalue="NA"

    if [[ -s "\$centrimo_tsv" ]]; then
      centrimo_top_evalue=\$(awk -F '\\t' 'NR == 2 {print \$5}' "\$centrimo_tsv" || true)
      centrimo_top_adj_pvalue=\$(awk -F '\\t' 'NR == 2 {print \$6}' "\$centrimo_tsv" || true)

      centrimo_top_evalue=\${centrimo_top_evalue:-NA}
      centrimo_top_adj_pvalue=\${centrimo_top_adj_pvalue:-NA}

      centrimo_has_central_motif=\$(awk -v p="\$centrimo_top_adj_pvalue" '
        BEGIN {
          if (p != "NA" && p + 0 < 0.05) print "yes";
          else print "no";
        }
      ')
    fi

    add_stat "centrimo_has_central_motif" "\$centrimo_has_central_motif"
    add_stat "centrimo_top_evalue" "\$centrimo_top_evalue"
    add_stat "centrimo_top_adj_pvalue" "\$centrimo_top_adj_pvalue"

  } 2> ${sample}.collect_sample_stats.err
  """
}

workflow {
  layout = readLayout()

  if (!(layout in ['PE', 'SE'])) {
    error "params.read_layout must be PE or SE"
  }

  if (layout == 'PE') {
    samples = Channel
      .fromPath(params.input, checkIfExists: true)
      .splitCsv(header: true)
      .map { row ->
        if (!row.fastq_2) {
          error "PE mode requires fastq_2 for sample ${row.sample}"
        }
        tuple(
          row.sample as String,
          row.role.toLowerCase() as String,
          file(row.fastq_1),
          file(row.fastq_2)
        )
      }

    FASTQC_RAW_PE(samples)
    MULTIQC_RAW_PE(FASTQC_RAW_PE.out.zip.collect())

    TRIMMOMATIC_PE(samples)
    trimmed_reads = TRIMMOMATIC_PE.out.reads

    FASTQC_TRIMMED_PE(trimmed_reads)
    MULTIQC_TRIMMED_PE(FASTQC_TRIMMED_PE.out.zip.collect())

    BOWTIE2_ALIGN_SORT_PE(trimmed_reads)
    aligned_bams = BOWTIE2_ALIGN_SORT_PE.out.bam

    SAMTOOLS_FILTER_Q20_PE(aligned_bams)
    filtered_bams = SAMTOOLS_FILTER_Q20_PE.out.bam

    SAMTOOLS_MARKDUP_REMOVE_PE(filtered_bams)
    markdup_bams = SAMTOOLS_MARKDUP_REMOVE_PE.out.bam
  }

  if (layout == 'SE') {
    samples = Channel
      .fromPath(params.input, checkIfExists: true)
      .splitCsv(header: true)
      .map { row ->
        tuple(
          row.sample as String,
          row.role.toLowerCase() as String,
          file(row.fastq_1)
        )
      }

    FASTQC_RAW_SE(samples)
    MULTIQC_RAW_SE(FASTQC_RAW_SE.out.zip.collect())

    TRIMMOMATIC_SE(samples)
    trimmed_reads = TRIMMOMATIC_SE.out.reads

    FASTQC_TRIMMED_SE(trimmed_reads)
    MULTIQC_TRIMMED_SE(FASTQC_TRIMMED_SE.out.zip.collect())

    BOWTIE2_ALIGN_SORT_SE(trimmed_reads)
    aligned_bams = BOWTIE2_ALIGN_SORT_SE.out.bam

    SAMTOOLS_FILTER_Q20_SE(aligned_bams)
    filtered_bams = SAMTOOLS_FILTER_Q20_SE.out.bam

    SAMTOOLS_MARKDUP_REMOVE_SE(filtered_bams)
    markdup_bams = SAMTOOLS_MARKDUP_REMOVE_SE.out.bam
  }

  control_bams = markdup_bams
    .filter { sample, role, bam, bam_index -> role == params.input_role }
    .map { sample, role, bam, bam_index -> bam }
    .collect()

  treatment_bams = markdup_bams
    .filter { sample, role, bam, bam_index -> role == params.tf_role }
    .map { sample, role, bam, bam_index -> tuple(sample, bam) }

  MACS3_CALLPEAK(treatment_bams, control_bams)
  peaks = MACS3_CALLPEAK.out.narrowpeak

  frip_input = treatment_bams.join(peaks)

  FRIP_CALCULATION(frip_input)
  frip_scores = FRIP_CALCULATION.out.frip

  PEAK_PROMOTER_INTERSECT(peaks)
  peak_promoter_intersects = PEAK_PROMOTER_INTERSECT.out.intersect

  TOP_SIGNAL_PEAKS(peaks)
  top_peaks = TOP_SIGNAL_PEAKS.out.top_peaks

  MEME_CHIP(top_peaks)
  meme_out = MEME_CHIP.out.meme_out

  sample_stats_input = frip_scores
    .join(peak_promoter_intersects)
    .join(peaks)
    .join(meme_out)

  COLLECT_SAMPLE_STATS(sample_stats_input)
}
