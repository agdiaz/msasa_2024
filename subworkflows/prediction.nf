nextflow.enable.dsl=2

workflow PREDICTION_WORKFLOW {
    take: directory_ch
    main:
        println "PREDICTION_WORKFLOW: $directory_ch"
        fasta_files = Channel.fromPath("${directory_ch}/*.tfa").view()

        // Channel.fromPath("${params.input_dir}/${directory}/${directory}.xml").view()
        // Channel.fromPath("${params.input_dir}/${directory}/${directory}.xml").view()
        // def reference_files = Channel.fromPath("${directoryectory}/*.reference.fasta").view()

        // alignCoincidences(fasta_files)
    emit:
        fasta_files = fasta_files.toList()
}
