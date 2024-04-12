nextflow.enable.dsl=2

// Define parameters
params.input_dir = "$baseDir/examples"  // default input directory
params.msasa_script = "$baseDir/src/msasa_cli.py"
params.results_dir = "$baseDir/msasa_predictions_results"

workflow {
    // Find all *.tfa files in the specified directory and subdirectories
    fasta_files = Channel.fromPath("${params.input_dir}/**/*.tfa", checkIfExists: true)

    alignCoincidences(fasta_files)
    alignIdentity(fasta_files)
    alignSimilarityBlosum62(fasta_files)
    alignSimilarityPam250(fasta_files)
    alignGlobal(fasta_files)
    alignLocal(fasta_files)

    computeBaliScore()
}

process alignCoincidences {
    label 'msasa'
    tag "${fasta_file.simpleName} ✅"
    publishDir params.results_dir + "/${fasta_file.simpleName}", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.coincidences.aln"), emit: alnFasta
        path("${fasta_file}.msasa.coincidences.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.coincidences.log"), emit: log
        path("${fasta_file}.msasa.coincidences.log.png"), emit: plot

    script:
    """
    python ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.coincidences.aln \
        --log_file ${fasta_file}.msasa.coincidences.log \
        --extend \
        --temp 1 \
        --cooling_rate 0.999 \
        --min_temp 0.001 \
        --max_no_changes 10000 \
        --quality_function coincidences
    """
}

process alignIdentity {
    label 'msasa'
    tag "${fasta_file.simpleName} ✅"
    publishDir params.results_dir + "/${fasta_file.simpleName}", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.identity.aln"), emit: alnFasta
        path("${fasta_file}.msasa.identity.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.identity.log"), emit: log
        path("${fasta_file}.msasa.identity.log.png"), emit: plot

    script:
    """
    python ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.identity.aln \
        --log_file ${fasta_file}.msasa.identity.log \
        --extend \
        --temp 1 \
        --cooling_rate 0.999 \
        --min_temp 0.0001 \
        --max_no_changes 25000 \
        --quality_function identity
    """
}

process alignSimilarityBlosum62 {
    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir params.results_dir + "/${fasta_file.simpleName}", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.similarity_blosum62.aln"), emit: alnFasta
        path("${fasta_file}.msasa.similarity_blosum62.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.similarity_blosum62.log"), emit: log
        path("${fasta_file}.msasa.similarity_blosum62.log.png"), emit: plot

    script:
    """
    python ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.similarity_blosum62.aln \
        --log_file ${fasta_file}.msasa.similarity_blosum62.log \
        --extend \
        --temp 5 \
        --cooling_rate 0.99 \
        --min_temp 0.001 \
        --max_no_changes 5000 \
        --quality_function similarity_blosum62
    """
}

process alignSimilarityPam250 {
    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir params.results_dir + "/${fasta_file.simpleName}", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.similarity_pam250.aln"), emit: alnFasta
        path("${fasta_file}.msasa.similarity_pam250.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.similarity_pam250.log"), emit: log
        path("${fasta_file}.msasa.similarity_pam250.log.png"), emit: plot

    script:
    """
    python ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.similarity_pam250.aln \
        --log_file ${fasta_file}.msasa.similarity_pam250.log \
        --extend \
        --temp 5 \
        --cooling_rate 0.99 \
        --min_temp 0.001 \
        --max_no_changes 5000 \
        --quality_function similarity_pam250
    """
}

process alignGlobal {
    label 'msasa'
    tag "${fasta_file.simpleName} ✅"
    publishDir params.results_dir + "/${fasta_file.simpleName}", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.global.aln"), emit: alnFasta
        path("${fasta_file}.msasa.global.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.global.log"), emit: log
        path("${fasta_file}.msasa.global.log.png"), emit: plot

    script:
    """
    python ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.global.aln \
        --log_file ${fasta_file}.msasa.global.log \
        --extend \
        --temp 1.5 \
        --cooling_rate 0.999 \
        --min_temp 0.001 \
        --max_no_changes 5000 \
        --quality_function global
    """
}

process alignLocal {
    label 'msasa'
    tag "${fasta_file.simpleName} ✅"
    publishDir params.results_dir + "/${fasta_file.simpleName}", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.local.aln"), emit: alnFasta
        path("${fasta_file}.msasa.local.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.local.log"), emit: log
        path("${fasta_file}.msasa.local.log.png"), emit: plot

    script:
    """
    python ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.local.aln \
        --log_file ${fasta_file}.msasa.local.log \
        --extend \
        --temp 1 \
        --cooling_rate 0.999 \
        --min_temp 0.001 \
        --max_no_changes 10000 \
        --quality_function local
    """
}

process computeBaliScore {
    tag "${alignment.simpleName}"
    publishDir "$baseDir/tesis_results/${reference.simpleName}", mode: 'copy'
    errorStrategy 'ignore'

    input:
        val predictorName
        path alignment

    output:
        path "${reference.simpleName}_baliscore_${predictorName}.txt", emit: score

    shell:
    """
    /software/bali-score/target/release/bali-score -t ${alignment} \
        -r ${reference} \
        -o ${reference.simpleName}_baliscore_${predictorName}.txt
    """
}

process computeCoreIndex {
    tag "${alignment.simpleName}"
    publishDir "$baseDir/tesis_results/${reference.simpleName}", mode: 'copy'
    errorStrategy 'ignore'

    input:
        val predictorName
        path alignment

    output:
        path "${reference.simpleName}_coreindex_${predictorName}.txt"

    script:
    """
    t_coffee -infile='$alignment' -output=score_ascii -score -outfile='${reference.simpleName}_coreindex_${predictorName}.txt'
    """
}

process computeTransitiveConsistencyScore {
    tag "${alignment.simpleName}"
    publishDir "$baseDir/tesis_results/${reference.simpleName}", mode: 'copy'
    errorStrategy 'ignore'

    input:
        val predictorName
        path alignment

    output:
        path "${reference.simpleName}_tcs_${predictorName}.txt"

    script:
    """
    t_coffee -infile '$alignment' -evaluate -output=score_ascii -outfile='${reference.simpleName}_tcs_${predictorName}.txt'
    """
}