process alignCoincidences {
    cpus 4
    memory '6 GB'
    time '20h'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.png", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.pdf", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln.clustal", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/logs/${fasta_file.simpleName}", pattern: "*.log", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.coincidences.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.coincidences.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.coincidences.*.log"), emit: log
        path("${fasta_file}.msasa.coincidences.*.log.png"), emit: plot
        path("${fasta_file}.msasa.coincidences.*.log.sns.png"), emit: plot_sns
        path("${fasta_file}.msasa.coincidences.*.log.pdf"), emit: plot_pdf
        path("${fasta_file}.msasa.coincidences.*.log.sns.pdf"), emit: plot_sns_pdf

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.coincidences.aln \
        --log_file ${fasta_file}.msasa.coincidences.log \
        --experiments_log_file ${fasta_file}.msasa.coincidences.experiments.log \
        --extend \
        --match_score 10 \
        --mismatch_score 0.5 \
        --gap_score 0.1 \
        --temp 1.0 \
        --cooling_rate 0.99 \
        --min_temp 0.00009 \
        --max_no_changes 2000 \
        --changes 10 \
        --iteration_neighbors 10 \
        --quality_function coincidences \
        --experiments ${params.experiments}
    """
}

process alignIdentity {
    cpus 4
    memory '6 GB'
    time '20h'

    label 'msasa'
    tag "${fasta_file.simpleName}"

    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.png", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.pdf", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln.clustal", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/logs/${fasta_file.simpleName}", pattern: "*.log", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.identity.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.identity.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.identity.*.log"), emit: log
        path("${fasta_file}.msasa.identity.*.log.png"), emit: plot
        path("${fasta_file}.msasa.identity.*.log.sns.png"), emit: plot_sns
        path("${fasta_file}.msasa.identity.*.log.pdf"), emit: plot_pdf
        path("${fasta_file}.msasa.identity.*.log.sns.pdf"), emit: plot_sns_pdf

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.identity.aln \
        --log_file ${fasta_file}.msasa.identity.log \
        --experiments_log_file ${fasta_file}.msasa.identity.experiments.log \
        --extend \
        --match_score 8.0 \
        --mismatch_score 0.50 \
        --gap_score 0.25 \
        --temp 1.0 \
        --cooling_rate 0.99 \
        --min_temp 0.00009 \
        --max_no_changes 2000 \
        --changes 10 \
        --iteration_neighbors 10 \
        --quality_function identity \
        --experiments ${params.experiments}
    """
}

process alignSimilarityBlosum62 {
    cpus 4
    memory '6 GB'
    time '20h'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.png", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.pdf", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln.clustal", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/logs/${fasta_file.simpleName}", pattern: "*.log", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.similarity_blosum62.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.similarity_blosum62.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.similarity_blosum62.*.log"), emit: log
        path("${fasta_file}.msasa.similarity_blosum62.*.log.png"), emit: plot
        path("${fasta_file}.msasa.similarity_blosum62.*.log.sns.png"), emit: plot_sns
        path("${fasta_file}.msasa.similarity_blosum62.*.log.pdf"), emit: plot_pdf
        path("${fasta_file}.msasa.similarity_blosum62.*.log.sns.pdf"), emit: plot_sns_pdf

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.similarity_blosum62.aln \
        --log_file ${fasta_file}.msasa.similarity_blosum62.log \
        --experiments_log_file ${fasta_file}.msasa.similarity_blosum62.experiments.log \
        --extend \
        --gap_score -4.0 \
        --mismatch_score -4.0 \
        --temp 1.0 \
        --cooling_rate 0.995 \
        --min_temp 0.00001 \
        --max_no_changes 2000 \
        --changes 10 \
        --iteration_neighbors 10 \
        --quality_function similarity_blosum62 \
        --experiments ${params.experiments}
    """
}

process alignSimilarityPam250 {
    cpus 4
    memory '6 GB'
    time '20h'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.png", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.pdf", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln.clustal", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/logs/${fasta_file.simpleName}", pattern: "*.log", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.similarity_pam250.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.similarity_pam250.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.similarity_pam250.*.log"), emit: log
        path("${fasta_file}.msasa.similarity_pam250.*.log.png"), emit: plot
        path("${fasta_file}.msasa.similarity_pam250.*.log.sns.png"), emit: plot_sns
        path("${fasta_file}.msasa.similarity_pam250.*.log.pdf"), emit: plot_pdf
        path("${fasta_file}.msasa.similarity_pam250.*.log.sns.pdf"), emit: plot_sns_pdf

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.similarity_pam250.aln \
        --log_file ${fasta_file}.msasa.similarity_pam250.log \
        --experiments_log_file ${fasta_file}.msasa.similarity_pam250.experiments.log \
        --extend \
        --gap_score -8.0 \
        --mismatch_score -8.0 \
        --temp 1.0 \
        --cooling_rate 0.995 \
        --min_temp 0.00001 \
        --max_no_changes 2000 \
        --changes 10 \
        --iteration_neighbors 10 \
        --quality_function similarity_pam250 \
        --experiments ${params.experiments}
    """
}

process alignGonnet92 {
    cpus 4
    memory '6 GB'
    time '20h'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.png", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.pdf", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln.clustal", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/logs/${fasta_file.simpleName}", pattern: "*.log", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.similarity_gonnet92.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.similarity_gonnet92.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.similarity_gonnet92.*.log"), emit: log
        path("${fasta_file}.msasa.similarity_gonnet92.*.log.png"), emit: plot
        path("${fasta_file}.msasa.similarity_gonnet92.*.log.sns.png"), emit: plot_sns
        path("${fasta_file}.msasa.similarity_gonnet92.*.log.pdf"), emit: plot_pdf
        path("${fasta_file}.msasa.similarity_gonnet92.*.log.sns.pdf"), emit: plot_sns_pdf

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.similarity_gonnet92.aln \
        --log_file ${fasta_file}.msasa.similarity_gonnet92.log \
        --experiments_log_file ${fasta_file}.msasa.similarity_gonnet92.experiments.log \
        --extend \
        --gap_score -6 \
        --mismatch_score -4.0 \
        --temp 0.50 \
        --cooling_rate 0.995 \
        --min_temp 0.00001 \
        --max_no_changes 2000 \
        --changes 20 \
        --iteration_neighbors 5 \
        --quality_function similarity_gonnet \
        --experiments ${params.experiments}
    """
}

process alignGlobal {
    cpus 4
    memory '6 GB'
    time '20h'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.png", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.pdf", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln.clustal", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/logs/${fasta_file.simpleName}", pattern: "*.log", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.global.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.global.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.global.*.log"), emit: log
        path("${fasta_file}.msasa.global.*.log.png"), emit: plot
        path("${fasta_file}.msasa.global.*.log.sns.png"), emit: plot_sns
        path("${fasta_file}.msasa.global.*.log.pdf"), emit: plot_pdf
        path("${fasta_file}.msasa.global.*.log.sns.pdf"), emit: plot_sns_pdf

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.global.aln \
        --log_file ${fasta_file}.msasa.global.log \
        --experiments_log_file ${fasta_file}.msasa.global.experiments.log \
        --extend \
        --match_score 10.0 \
        --mismatch_score -1.0 \
        --gap_score -4.0 \
        --temp 1.0 \
        --cooling_rate 0.99 \
        --min_temp 0.00009 \
        --max_no_changes 2000 \
        --changes 20 \
        --iteration_neighbors 5 \
        --quality_function global \
        --experiments ${params.experiments}
    """
}

process alignLocal {
    cpus 4
    memory '6 GB'
    time '20h'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.png", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/plots/${fasta_file.simpleName}", pattern: "*.pdf", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/predictions/${fasta_file.simpleName}", pattern: "*.aln.clustal", overwrite: true, mode: 'copy'
    publishDir "${params.results_dir}/logs/${fasta_file.simpleName}", pattern: "*.log", overwrite: true, mode: 'copy'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.local.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.local.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.local.*.log"), emit: log
        path("${fasta_file}.msasa.local.*.log.png"), emit: plot
        path("${fasta_file}.msasa.local.*.log.sns.png"), emit: plot_sns
        path("${fasta_file}.msasa.local.*.log.pdf"), emit: plot_pdf
        path("${fasta_file}.msasa.local.*.log.sns.pdf"), emit: plot_sns_pdf

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.local.aln \
        --log_file ${fasta_file}.msasa.local.log \
        --experiments_log_file ${fasta_file}.msasa.local.experiments.log \
        --extend \
        --match_score 5.0 \
        --mismatch_score 0.5 \
        --gap_score -3.0 \
        --temp 1.0 \
        --cooling_rate 0.99 \
        --min_temp 0.00009 \
        --max_no_changes 2000 \
        --changes 20 \
        --iteration_neighbors 5 \
        --quality_function local \
        --experiments ${params.experiments}
    """
}