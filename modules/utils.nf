params.sort_script = "$baseDir/src/sort_fasta.py"

process sortFasta {
    tag "${fasta_file.simpleName}"
    publishDir params.results_dir + "/references/${ref_file.simpleName}", overwrite: true, mode: 'copy'

    input:
        path fasta_file
        path ref_file

    output:
        path "${fasta_file}.sorted.fasta", emit: sorted_seqs

    script:
    """
    python ${params.sort_script} $fasta_file $ref_file ${fasta_file}.sorted.fasta
    """
}
