#!/usr/bin/env nextflow
params.s3Path = "s3://dovetail-public/NF-test/fastqs/"
params.outDir = "./tmp"
params.genome = 'hg38'


Channel
    .fromFilePairs("s3://dovetail-public/NF-test/genome/nexflow/${params.genome}/*.{amb,sa,pac,ann,bwt,fa}", size: -1, checkIfExists: true)
    .set { bwa_index }

Channel
    .fromFilePairs("${params.s3Path}/*_{R1,R2}*.fastq.gz",
        flat: true)
    .set{s3_ch}


process bwa_mem {
    tag "${id}"
    label "cpu"
    container 'dovetailg/bwa-samtools'

    publishDir "${params.outDir}/fastqs",
         saveAs: {filename -> filename.endsWith('.fastq.gz') ? filename : null}
    
    input:
    tuple val(id), path(R1s), path(R2s) from s3_ch
	    .map { id, file1, file2 -> tuple(file1.name.toString().split('_')[0], file1, file2) }
	    .groupTuple()
    
    tuple val(index), path(index_files) from bwa_index.first()
    
    output:
    tuple val(id), path("*.sam") into sam4chrSize, sam_ch1
    tuple val(id), path(R1s), path(R2s) into fq_stat_ch
    
    script:
    """
    bwa mem -5SP -t ${task.cpus} \
        ${index} \
        ${R1s} \
        ${R2s} \
        > ${id}.sam
    """
}

process bam_sort {
    tag "${id}"
    label "cpu"
    container 'dovetailg/bwa-samtools'
    
    // publishDir "${params.outDir}/bam",
    //     saveAs: {filename -> filename.endsWith('.bam') ? filename : null}

    input:
    tuple val(id), path(bam) from sam_ch1
    
    output:
    tuple val(id), path("*.bam"), path("*.bai") into bam_ch, bam_preseq_ch, bam_bigwig_ch

    script:
    """
    samtools sort -m 2G \
	-@ ${task.cpus} \
	-o ${id}.bam \
	${bam} \
	&& samtools index -@${task.cpus} ${id}.bam
    """
}

process bam2bw {
    tag "_${id}"
    label "major"
    container 'dovetailg/r-cov'
    
    //publishDir "${params.outDir}/bigwigs"
    
    input:
    tuple val(id), path(bam), path(idx) from bam_bigwig_ch
    
    output:
    tuple val(id), path ("*.bw") into bigwig_ch
    
    script:
    """
    bam2bw ${bam} ${id}.bw ${task.cpus}
    """
}
