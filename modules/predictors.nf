process alignCoincidences {
    cpus 1
    memory '4 GB'
    time '5d'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.png", overwrite: true, mode: 'copy'
    // publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.pdf", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln.clustal", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/logs/${fasta_file.simpleName}", pattern: "*.log", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)
        val strict

    output:
        path("${fasta_file}.msasa.coincidences.${strict ? "strict" : "free" }.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.coincidences.${strict ? "strict" : "free" }.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.coincidences.${strict ? "strict" : "free" }.*.log"), emit: log
        path("${fasta_file}.msasa.coincidences.${strict ? "strict" : "free" }.*.log.png"), emit: plot
        path("${fasta_file}.msasa.coincidences.${strict ? "strict" : "free" }.*.log.sns.png"), emit: plot_sns
        // path("${fasta_file}.msasa.coincidences.${strict ? "strict" : "free" }.*.log.pdf"), emit: plot_pdf
        // path("${fasta_file}.msasa.coincidences.${strict ? "strict" : "free" }.*.log.sns.pdf"), emit: plot_sns_pdf

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.coincidences.${strict ? "strict" : "free" }.aln \
        --log_file ${fasta_file}.msasa.coincidences.${strict ? "strict" : "free" }.log \
        --experiments_log_file ${fasta_file}.msasa.coincidences.experiments.${strict ? "strict" : "free" }.log \
        --extend \
        --match_score 10 \
        --mismatch_score 0.75 \
        --gap_score 0.1 \
        --temp 0.90 \
        --cooling_rate 0.99 \
        --min_temp 0.00001 \
        --max_no_changes 2000 \
        --changes 10 \
        --iteration_neighbors 10 \
        --quality_function coincidences \
        --experiments ${params.experiments} \
        ${strict ? "--strict_mode" : "" }
    """
}

process alignIdentity {
    cpus 1
    memory '4 GB'
    time '5d'

    label 'msasa'
    tag "${fasta_file.simpleName}"

    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.png", overwrite: true, mode: 'copy'
    // publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.pdf", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln.clustal", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/logs/${fasta_file.simpleName}", pattern: "*.log", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)
        val strict

    output:
        path("${fasta_file}.msasa.identity.${strict ? "strict" : "free" }.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.identity.${strict ? "strict" : "free" }.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.identity.${strict ? "strict" : "free" }.*.log"), emit: log
        path("${fasta_file}.msasa.identity.${strict ? "strict" : "free" }.*.log.png"), emit: plot
        path("${fasta_file}.msasa.identity.${strict ? "strict" : "free" }.*.log.sns.png"), emit: plot_sns
        // path("${fasta_file}.msasa.identity.${strict ? "strict" : "free" }.*.log.pdf"), emit: plot_pdf
        // path("${fasta_file}.msasa.identity.${strict ? "strict" : "free" }.*.log.sns.pdf"), emit: plot_sns_pdf

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.identity.${strict ? "strict" : "free" }.aln \
        --log_file ${fasta_file}.msasa.identity.${strict ? "strict" : "free" }.log \
        --experiments_log_file ${fasta_file}.msasa.identity.experiments.${strict ? "strict" : "free" }.log \
        --extend \
        --match_score 8.0 \
        --mismatch_score 1.0 \
        --gap_score 2.0 \
        --temp 5 \
        --cooling_rate 0.99 \
        --min_temp 0.00009 \
        --max_no_changes 2000 \
        --changes 10 \
        --iteration_neighbors 5 \
        --quality_function identity \
        --experiments ${params.experiments} \
        ${strict ? "--strict_mode" : "" }
    """
}

process alignSimilarityBlosum62 {
    cpus 1
    memory '4 GB'
    time '5d'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.png", overwrite: true, mode: 'copy'
    // publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.pdf", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln.clustal", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/logs/${fasta_file.simpleName}", pattern: "*.log", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)
        val strict

    output:
        path("${fasta_file}.msasa.similarity_blosum62.${strict ? "strict" : "free" }.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.similarity_blosum62.${strict ? "strict" : "free" }.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.similarity_blosum62.${strict ? "strict" : "free" }.*.log"), emit: log
        path("${fasta_file}.msasa.similarity_blosum62.${strict ? "strict" : "free" }.*.log.png"), emit: plot
        path("${fasta_file}.msasa.similarity_blosum62.${strict ? "strict" : "free" }.*.log.sns.png"), emit: plot_sns
        // path("${fasta_file}.msasa.similarity_blosum62.${strict ? "strict" : "free" }.*.log.pdf"), emit: plot_pdf
        // path("${fasta_file}.msasa.similarity_blosum62.${strict ? "strict" : "free" }.*.log.sns.pdf"), emit: plot_sns_pdf

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.similarity_blosum62.${strict ? "strict" : "free" }.aln \
        --log_file ${fasta_file}.msasa.similarity_blosum62.${strict ? "strict" : "free" }.log \
        --experiments_log_file ${fasta_file}.msasa.similarity_blosum62.experiments.${strict ? "strict" : "free" }.log \
        --extend \
        --gap_score -4.0 \
        --mismatch_score -4.0 \
        --temp 1.0 \
        --cooling_rate 0.9925 \
        --min_temp 0.00001 \
        --max_no_changes 2000 \
        --changes 10 \
        --iteration_neighbors 5 \
        --quality_function similarity_blosum62 \
        --experiments ${params.experiments} \
        ${strict ? "--strict_mode" : "" }
    """
}

process alignSimilarityPam250 {
    cpus 1
    memory '4 GB'
    time '5d'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.png", overwrite: true, mode: 'copy'
    // publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.pdf", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln.clustal", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/logs/${fasta_file.simpleName}", pattern: "*.log", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)
        val strict

    output:
        path("${fasta_file}.msasa.similarity_pam250.${strict ? "strict" : "free" }.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.similarity_pam250.${strict ? "strict" : "free" }.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.similarity_pam250.${strict ? "strict" : "free" }.*.log"), emit: log
        path("${fasta_file}.msasa.similarity_pam250.${strict ? "strict" : "free" }.*.log.png"), emit: plot
        path("${fasta_file}.msasa.similarity_pam250.${strict ? "strict" : "free" }.*.log.sns.png"), emit: plot_sns
        // path("${fasta_file}.msasa.similarity_pam250.${strict ? "strict" : "free" }.*.log.pdf"), emit: plot_pdf
        // path("${fasta_file}.msasa.similarity_pam250.${strict ? "strict" : "free" }.*.log.sns.pdf"), emit: plot_sns_pdf

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.similarity_pam250.${strict ? "strict" : "free" }.aln \
        --log_file ${fasta_file}.msasa.similarity_pam250.${strict ? "strict" : "free" }.log \
        --experiments_log_file ${fasta_file}.msasa.similarity_pam250.experiments.${strict ? "strict" : "free" }.log \
        --extend \
        --gap_score -8.0 \
        --mismatch_score -8.0 \
        --temp 1.0 \
        --cooling_rate 0.9925 \
        --min_temp 0.00001 \
        --max_no_changes 2000 \
        --changes 10 \
        --iteration_neighbors 5 \
        --quality_function similarity_pam250 \
        --experiments ${params.experiments} \
        ${strict ? "--strict_mode" : "" }
    """
}

process alignGonnet92 {
    cpus 1
    memory '4 GB'
    time '5d'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.png", overwrite: true, mode: 'copy'
    // publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.pdf", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln.clustal", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/logs/${fasta_file.simpleName}", pattern: "*.log", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)
        val strict

    output:
        path("${fasta_file}.msasa.similarity_gonnet92.${strict ? "strict" : "free" }.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.similarity_gonnet92.${strict ? "strict" : "free" }.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.similarity_gonnet92.${strict ? "strict" : "free" }.*.log"), emit: log
        path("${fasta_file}.msasa.similarity_gonnet92.${strict ? "strict" : "free" }.*.log.png"), emit: plot
        path("${fasta_file}.msasa.similarity_gonnet92.${strict ? "strict" : "free" }.*.log.sns.png"), emit: plot_sns
        // path("${fasta_file}.msasa.similarity_gonnet92.${strict ? "strict" : "free" }.*.log.pdf"), emit: plot_pdf
        // path("${fasta_file}.msasa.similarity_gonnet92.${strict ? "strict" : "free" }.*.log.sns.pdf"), emit: plot_sns_pdf

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.similarity_gonnet92.${strict ? "strict" : "free" }.aln \
        --log_file ${fasta_file}.msasa.similarity_gonnet92.${strict ? "strict" : "free" }.log \
        --experiments_log_file ${fasta_file}.msasa.similarity_gonnet92.experiments.${strict ? "strict" : "free" }.log \
        --extend \
        --gap_score -6 \
        --mismatch_score -4.0 \
        --temp 1.0 \
        --cooling_rate 0.995 \
        --min_temp 0.00001 \
        --max_no_changes 2000 \
        --changes 10 \
        --iteration_neighbors 5 \
        --quality_function similarity_gonnet \
        --experiments ${params.experiments} \
        ${strict ? "--strict_mode" : "" }
    """
}

process alignGlobal {
    cpus 1
    memory '4 GB'
    time '5d'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.png", overwrite: true, mode: 'copy'
    // publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.pdf", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln.clustal", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/logs/${fasta_file.simpleName}", pattern: "*.log", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)
        val strict

    output:
        path("${fasta_file}.msasa.global.${strict ? "strict" : "free" }.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.global.${strict ? "strict" : "free" }.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.global.${strict ? "strict" : "free" }.*.log"), emit: log
        path("${fasta_file}.msasa.global.${strict ? "strict" : "free" }.*.log.png"), emit: plot
        path("${fasta_file}.msasa.global.${strict ? "strict" : "free" }.*.log.sns.png"), emit: plot_sns
        // path("${fasta_file}.msasa.global.${strict ? "strict" : "free" }.*.log.pdf"), emit: plot_pdf
        // path("${fasta_file}.msasa.global.${strict ? "strict" : "free" }.*.log.sns.pdf"), emit: plot_sns_pdf

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.global.${strict ? "strict" : "free" }.aln \
        --log_file ${fasta_file}.msasa.global.${strict ? "strict" : "free" }.log \
        --experiments_log_file ${fasta_file}.msasa.global.experiments.${strict ? "strict" : "free" }.log \
        --extend \
        --match_score 10.0 \
        --mismatch_score -1.0 \
        --gap_score -4.0 \
        --temp 0.80 \
        --cooling_rate 0.995 \
        --min_temp 0.00005 \
        --max_no_changes 2000 \
        --changes 10 \
        --iteration_neighbors 5 \
        --quality_function global \
        --experiments ${params.experiments} \
        ${strict ? "--strict_mode" : "" }
    """
}

process alignLocal {
    cpus 4
    memory '4 GB'
    time '5d'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.png", overwrite: true, mode: 'copy'
    // publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.pdf", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln.clustal", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/logs/${fasta_file.simpleName}", pattern: "*.log", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)
        val strict

    output:
        path("${fasta_file}.msasa.local.${strict ? "strict" : "free" }.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.local.${strict ? "strict" : "free" }.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.local.${strict ? "strict" : "free" }.*.log"), emit: log
        path("${fasta_file}.msasa.local.${strict ? "strict" : "free" }.*.log.png"), emit: plot
        path("${fasta_file}.msasa.local.${strict ? "strict" : "free" }.*.log.sns.png"), emit: plot_sns
        // path("${fasta_file}.msasa.local.${strict ? "strict" : "free" }.*.log.pdf"), emit: plot_pdf
        // path("${fasta_file}.msasa.local.${strict ? "strict" : "free" }.*.log.sns.pdf"), emit: plot_sns_pdf

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.local.${strict ? "strict" : "free" }.aln \
        --log_file ${fasta_file}.msasa.local.${strict ? "strict" : "free" }.log \
        --experiments_log_file ${fasta_file}.msasa.local.experiments.${strict ? "strict" : "free" }.log \
        --extend \
        --match_score 7.0 \
        --mismatch_score 1.0 \
        --gap_score -3.0 \
        --temp 0.80 \
        --cooling_rate 0.995 \
        --min_temp 0.00009 \
        --max_no_changes 2000 \
        --changes 10 \
        --iteration_neighbors 5 \
        --quality_function local \
        --experiments ${params.experiments} \
        ${strict ? "--strict_mode" : "" }
    """
}
