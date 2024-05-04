nextflow.enable.dsl = 2

params.reads = "$baseDir/data/*_{1,2}.fastq"
params.transcriptome_file = "$baseDir/data/Homo_sapiens.GRCh38.cdna.all.fa"
params.multiqc = "$baseDir/multiqc"
params.outdir = "results"

log.info """\
        R N A S E Q - N F  P I P E L I N E
        ===================================
        Created by Adebowale
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

workflow {
    // Channel declaration
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set {read_pairs_ch}
    //
    index_ch = INDEX(params.transcriptome_file)
    index_ch.view()
    read_pairs_ch.view()
}

// Expression Quantification