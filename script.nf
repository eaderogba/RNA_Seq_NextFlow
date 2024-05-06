nextflow.enable.dsl = 2

params.reads = "$baseDir/data/*_{1,2}.fastq"
params.transcriptome_file = "$baseDir/data/Homo_sapiens.GRCh38.cdna.all.fa"
params.multiqc = "$baseDir/multiqc"
params.outdir = "results"

log.info """\
        R N A S E Q - N F  P I P E L I N E
        ===================================
        Adapted by Adebowale
        ===================================
        transcriptome: ${params.transcriptome_file}
        reads        : ${params.reads}
        outdir       : ${params.outdir}
        """
        .stripIndent(true)

// Create transcriptome index file
/*
INDEX process defined
Binary index created using transcriptome file
*/
process INDEX {
    cpus 4

    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}

// Quality Control
process FASTQC {
    cpus 4
    tag "FASTQC on $sample_id"
    publishDir params.outdir

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

// MULTIQC
process MULTIQC {
    cpus 4

    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

// Expression Quantification
process QUANTIFICATION {
    cpus 4
    tag "Salmon on $sample_id"
    publishDir params.outdir, mode: 'copy'

    input:
    path salmon_index
    tuple val(sample_id), path(reads)

    output:
    path "$sample_id"

    script:
    """
    salmon quant --threads $task.cpus --libType=u -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
    """
}

workflow {
    // Channel declaration
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set {read_pairs_ch}
    //
    index_ch = INDEX(params.transcriptome_file)
    fastqc_ch = FASTQC(read_pairs_ch)
    multiqc_ch = MULTIQC(fastqc_ch.collect())
    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
}

workflow.onComplete {
    log.info (workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong")
}
