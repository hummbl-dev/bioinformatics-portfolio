#!/usr/bin/env nextflow
// Germline short-variant pipeline: FASTQ -> VCF + annotation + QC report.
// Tested entry point:  nextflow run main.nf -profile test,conda
nextflow.enable.dsl = 2

params.reads     = "data/test/*_R{1,2}.fastq"
params.reference = "data/test/reference.fa"
params.snpeff_db = null            // e.g. "hg38"; null = skip annotation
params.min_qual  = 20
params.outdir    = "results"

process FASTQC {
    conda 'bioconda::fastqc=0.12.1'
    publishDir "${params.outdir}/qc", mode: 'copy'
    input:  tuple val(sample), path(reads)
    output: path "*_fastqc.*"
    script:
    """
    fastqc -t 2 ${reads.join(' ')}
    """
}

process TRIM {
    conda 'bioconda::fastp=0.23.4'
    input:  tuple val(sample), path(reads)
    output:
        tuple val(sample), path("${sample}_trim_R{1,2}.fastq.gz"), emit: reads
        path "${sample}.fastp.html", emit: report
    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
        -o ${sample}_trim_R1.fastq.gz -O ${sample}_trim_R2.fastq.gz \
        -h ${sample}.fastp.html
    """
}

process INDEX_REF {
    conda 'bioconda::bwa=0.7.18 bioconda::samtools=1.20'
    input:  path ref
    output: tuple path(ref), path("${ref}.*")
    script:
    """
    bwa index ${ref}
    samtools faidx ${ref}
    """
}

process ALIGN {
    conda 'bioconda::bwa=0.7.18 bioconda::samtools=1.20'
    input:
        tuple val(sample), path(reads)
        tuple path(ref), path(refidx)
    output: tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai")
    script:
    """
    bwa mem -t 2 ${ref} ${reads[0]} ${reads[1]} | \
        samtools sort -o ${sample}.sorted.bam -
    samtools markdup -r ${sample}.sorted.bam ${sample}.bam
    samtools index ${sample}.bam
    """
}

process CALL_VARIANTS {
    conda 'bioconda::bcftools=1.20'
    input:
        tuple val(sample), path(bam), path(bai)
        tuple path(ref), path(refidx)
    output: tuple val(sample), path("${sample}.raw.vcf")
    script:
    """
    bcftools mpileup -f ${ref} ${bam} | bcftools call -mv -Ov -o ${sample}.raw.vcf
    """
}

process FILTER_VARIANTS {
    conda 'bioconda::bcftools=1.20'
    input:  tuple val(sample), path(vcf)
    output: tuple val(sample), path("${sample}.filtered.vcf")
    script:
    """
    bcftools view -i 'QUAL>=${params.min_qual}' -Ov -o ${sample}.filtered.vcf ${vcf}
    """
}

process ANNOTATE {
    conda 'bioconda::snpeff=5.2'
    when: params.snpeff_db != null
    input:  tuple val(sample), path(vcf)
    output: tuple val(sample), path("${sample}.annotated.vcf")
    script:
    """
    snpEff ${params.snpeff_db} ${vcf} > ${sample}.annotated.vcf
    """
}

process MULTIQC {
    conda 'bioconda::multiqc=1.25'
    publishDir params.outdir, mode: 'copy'
    input:  path qc_files
    output: path "multiqc_report.html"
    script:
    """
    multiqc .
    """
}

workflow {
    reads_ch = Channel.fromFilePairs(params.reads, size: 2)
    ref_ch   = Channel.fromPath(params.reference)

    FASTQC(reads_ch)
    TRIM(reads_ch)
    INDEX_REF(ref_ch)
    ALIGN(TRIM.out.reads, INDEX_REF.out)
    CALL_VARIANTS(ALIGN.out, INDEX_REF.out)
    FILTER_VARIANTS(CALL_VARIANTS.out)
    ANNOTATE(FILTER_VARIANTS.out)
    MULTIQC(FASTQC.out.collect().concat(TRIM.out.report.collect()))
}
