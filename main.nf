// nextflow.enable.dsl=2

// Define parameters
params.input_dir = "$baseDir/datasets/BB11004"  // default input directory
params.msasa_script = "$baseDir/src/msasa_cli.py"
params.mumsa_script = "$baseDir/src/mumsa_plot.py"
params.results_dir = "$baseDir/thesis_results"
params.experiments = 5

include {
    computeMumsaOverlapScore;

    computeCoreIndex as computeReferenceCoreIndex;
    computeCoreIndex as computeCoincidencesCoreIndex;
    computeCoreIndex as computeIdentityCoreIndex;
    computeCoreIndex as computeSimilarityBlosum62CoreIndex;
    computeCoreIndex as computeSimilarityPam250CoreIndex;
    computeCoreIndex as computeSimilarityGonnet92CoreIndex;
    computeCoreIndex as computeGlobalCoreIndex;
    computeCoreIndex as computeLocalCoreIndex;

    computeTransitiveConsistencyScore as computeReferenceTransitiveConsistencyScore;
    computeTransitiveConsistencyScore as computeCoincidencesTransitiveConsistencyScore;
    computeTransitiveConsistencyScore as computeIdentityTransitiveConsistencyScore;
    computeTransitiveConsistencyScore as computeSimilarityBlosum62TransitiveConsistencyScore;
    computeTransitiveConsistencyScore as computeSimilarityPam250TransitiveConsistencyScore;
    computeTransitiveConsistencyScore as computeSimilarityGonnet92TransitiveConsistencyScore;
    computeTransitiveConsistencyScore as computeGlobalTransitiveConsistencyScore;
    computeTransitiveConsistencyScore as computeLocalTransitiveConsistencyScore;

    computeBaliScore as computeCoincidencesBaliScore;
    computeBaliScore as computeIdentityBaliScore;
    computeBaliScore as computeSimilarityBlosum62BaliScore;
    computeBaliScore as computeSimilarityPam250BaliScore;
    computeBaliScore as computeSimilarityGonnet92BaliScore;
    computeBaliScore as computeGlobalBaliScore;
    computeBaliScore as computeLocalBaliScore;

} from "./modules/scores.nf"

workflow {
    fasta_files = Channel.fromPath("${params.input_dir}/*.tfa", checkIfExists: true)
    reference_fasta_files = Channel.fromPath("${params.input_dir}/*.reference.fasta", checkIfExists: true)
    bali_base_xml_files = Channel.fromPath("${params.input_dir}/*.xml", checkIfExists: true)

    clustal_fasta_files = Channel.fromPath("${params.input_dir}/*_clustal.fasta", checkIfExists: true)
    kaling_fasta_files = Channel.fromPath("${params.input_dir}/*_kalign.fasta", checkIfExists: true)
    mafft_fasta_files = Channel.fromPath("${params.input_dir}/*_mafft.fasta", checkIfExists: true)
    muscle_fasta_files = Channel.fromPath("${params.input_dir}/*_muscle.fasta", checkIfExists: true)
    t_coffee_fasta_files = Channel.fromPath("${params.input_dir}/*_tCoffee.fasta", checkIfExists: true)

    alignIdentity(fasta_files)
    alignCoincidences(fasta_files)
    alignSimilarityBlosum62(fasta_files)
    alignSimilarityPam250(fasta_files)
    alignGonnet92(fasta_files)
    alignGlobal(fasta_files)
    alignLocal(fasta_files)

    alignCoincidences.out.alnFasta.flatten().randomSample(1).set { coincidencesSample }
    alignIdentity.out.alnFasta.flatten().randomSample(1).set { identitySample }
    alignSimilarityBlosum62.out.alnFasta.flatten().randomSample(1).set { blosum62Sample }
    alignSimilarityPam250.out.alnFasta.flatten().randomSample(1).set { pam250Sample }
    alignGonnet92.out.alnFasta.flatten().randomSample(1).set { gonnet92Sample }
    alignGlobal.out.alnFasta.flatten().randomSample(1).set { globalSample }
    alignLocal.out.alnFasta.flatten().randomSample(1).set { localSample }

    computeMumsaOverlapScore(
        reference_fasta_files,
        muscle_fasta_files,
        clustal_fasta_files,
        mafft_fasta_files,
        t_coffee_fasta_files,
        kaling_fasta_files,
        coincidencesSample,
        identitySample,
        blosum62Sample,
        pam250Sample,
        gonnet92Sample,
        globalSample,
        localSample
    )

    plotMumsaOverlap(computeMumsaOverlapScore.out)

    computeReferenceCoreIndex(reference_fasta_files, "reference")
    computeCoincidencesCoreIndex(coincidencesSample, "coincidences")
    computeIdentityCoreIndex(identitySample, "identity")
    computeSimilarityBlosum62CoreIndex(blosum62Sample, "similarity_blosum62")
    computeSimilarityPam250CoreIndex(pam250Sample, "similarity_pam250")
    computeSimilarityGonnet92CoreIndex(gonnet92Sample, "similarity_pam250")
    computeGlobalCoreIndex(globalSample, "global")
    computeLocalCoreIndex(localSample, "local")

    computeReferenceTransitiveConsistencyScore(reference_fasta_files, "reference")
    computeCoincidencesTransitiveConsistencyScore(coincidencesSample, "coincidences")
    computeIdentityTransitiveConsistencyScore(identitySample, "identity")
    computeSimilarityBlosum62TransitiveConsistencyScore(blosum62Sample, "similarity_blosum62")
    computeSimilarityPam250TransitiveConsistencyScore(pam250Sample, "similarity_pam250")
    computeSimilarityGonnet92TransitiveConsistencyScore(gonnet92Sample, "similarity_gonnet92")
    computeGlobalTransitiveConsistencyScore(globalSample, "global")
    computeLocalTransitiveConsistencyScore(localSample, "local")

    computeCoincidencesBaliScore(bali_base_xml_files, coincidencesSample, "coincidences")
    computeIdentityBaliScore(bali_base_xml_files, identitySample, "identity")
    computeSimilarityBlosum62BaliScore(bali_base_xml_files, blosum62Sample, "similarity_blosum62")
    computeSimilarityPam250BaliScore(bali_base_xml_files, pam250Sample, "similarity_pam250")
    computeSimilarityGonnet92BaliScore(bali_base_xml_files, gonnet92Sample, "similarity_gonnet92")
    computeGlobalBaliScore(bali_base_xml_files, globalSample, "global")
    computeLocalBaliScore(bali_base_xml_files, localSample, "local")
}

process alignCoincidences {
    cpus 2
    memory '6 GB'
    time '20h'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir params.results_dir + "/predictions/${fasta_file.simpleName}", overwrite: true, mode: 'link'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.coincidences.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.coincidences.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.coincidences.*.log"), emit: log
        path("${fasta_file}.msasa.coincidences.*.log.png"), emit: plot

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.coincidences.aln \
        --log_file ${fasta_file}.msasa.coincidences.log \
        --experiments_log_file ${fasta_file}.msasa.coincidences.experiments.log \
        --extend \
        --match_score 10 \
        --mismatch_score 0.5 \
        --gap_score 0.1 \
        --temp 0.10 \
        --cooling_rate 0.999 \
        --min_temp 0.000001 \
        --max_no_changes 10000 \
        --changes 25 \
        --iteration_neighbors 100 \
        --quality_function coincidences \
        --experiments ${params.experiments}
    """
}

process alignIdentity {
    cpus 2
    memory '6 GB'
    time '20h'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir params.results_dir + "/predictions/${fasta_file.simpleName}", overwrite: true, mode: 'link'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.identity.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.identity.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.identity.*.log"), emit: log
        path("${fasta_file}.msasa.identity.*.log.png"), emit: plot

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.identity.aln \
        --log_file ${fasta_file}.msasa.identity.log \
        --experiments_log_file ${fasta_file}.msasa.identity.experiments.log \
        --extend \
        --match_score 10 \
        --mismatch_score 1.0 \
        --gap_score 0.0 \
        --temp 1.0 \
        --cooling_rate 0.995 \
        --min_temp 0.000001 \
        --max_no_changes 10000 \
        --changes 50 \
        --iteration_neighbors 10 \
        --quality_function identity \
        --experiments ${params.experiments}
    """
}

process alignSimilarityBlosum62 {
    cpus 2
    memory '6 GB'
    time '20h'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir params.results_dir + "/predictions/${fasta_file.simpleName}", overwrite: true, mode: 'link'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.similarity_blosum62.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.similarity_blosum62.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.similarity_blosum62.*.log"), emit: log
        path("${fasta_file}.msasa.similarity_blosum62.*.log.png"), emit: plot

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.similarity_blosum62.aln \
        --log_file ${fasta_file}.msasa.similarity_blosum62.log \
        --experiments_log_file ${fasta_file}.msasa.similarity_blosum62.experiments.log \
        --extend \
        --gap_score -8.0 \
        --mismatch_score -5.0 \
        --temp 1.0 \
        --cooling_rate 0.999 \
        --min_temp 0.00001 \
        --max_no_changes 5000 \
        --changes 50 \
        --iteration_neighbors 10 \
        --quality_function similarity_blosum62 \
        --experiments ${params.experiments}
    """
}

process alignSimilarityPam250 {
    cpus 2
    memory '6 GB'
    time '20h'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir params.results_dir + "/predictions/${fasta_file.simpleName}", overwrite: true, mode: 'link'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.similarity_pam250.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.similarity_pam250.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.similarity_pam250.*.log"), emit: log
        path("${fasta_file}.msasa.similarity_pam250.*.log.png"), emit: plot

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.similarity_pam250.aln \
        --log_file ${fasta_file}.msasa.similarity_pam250.log \
        --experiments_log_file ${fasta_file}.msasa.similarity_pam250.experiments.log \
        --extend \
        --gap_score -10.0 \
        --mismatch_score -1.0 \
        --temp 1.0 \
        --cooling_rate 0.999 \
        --min_temp 0.00001 \
        --max_no_changes 5000 \
        --changes 500 \
        --iteration_neighbors 10 \
        --quality_function similarity_pam250 \
        --experiments ${params.experiments}
    """
}

process alignGonnet92 {
    cpus 2
    memory '6 GB'
    time '20h'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir params.results_dir + "/predictions/${fasta_file.simpleName}", overwrite: true, mode: 'link'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.similarity_gonnet92.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.similarity_gonnet92.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.similarity_gonnet92.*.log"), emit: log
        path("${fasta_file}.msasa.similarity_gonnet92.*.log.png"), emit: plot

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.similarity_gonnet92.aln \
        --log_file ${fasta_file}.msasa.similarity_gonnet92.log \
        --experiments_log_file ${fasta_file}.msasa.similarity_gonnet92.experiments.log \
        --extend \
        --gap_score -10.0 \
        --mismatch_score -2.0 \
        --temp 1.0 \
        --cooling_rate 0.999 \
        --min_temp 0.00001 \
        --max_no_changes 5000 \
        --changes 25 \
        --iteration_neighbors 10 \
        --quality_function similarity_gonnet \
        --experiments ${params.experiments}
    """
}

process alignGlobal {
    cpus 2
    memory '6 GB'
    time '20h'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir params.results_dir + "/predictions/${fasta_file.simpleName}", overwrite: true, mode: 'link'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.global.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.global.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.global.*.log"), emit: log
        path("${fasta_file}.msasa.global.*.log.png"), emit: plot

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.global.aln \
        --log_file ${fasta_file}.msasa.global.log \
        --experiments_log_file ${fasta_file}.msasa.global.experiments.log \
        --extend \
        --match_score 10.0 \
        --mismatch_score -1.0 \
        --gap_score -2.5 \
        --temp 0.1 \
        --cooling_rate 0.992 \
        --min_temp 0.00001 \
        --max_no_changes 7500 \
        --changes 50 \
        --iteration_neighbors 10 \
        --quality_function global \
        --experiments ${params.experiments}
    """
}

process alignLocal {
    cpus 2
    memory '6 GB'
    time '20h'

    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir params.results_dir + "/predictions/${fasta_file.simpleName}", overwrite: true, mode: 'link'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.local.*.aln"), emit: alnFasta
        path("${fasta_file}.msasa.local.*.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.local.*.log"), emit: log
        path("${fasta_file}.msasa.local.*.log.png"), emit: plot

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.local.aln \
        --log_file ${fasta_file}.msasa.local.log \
        --experiments_log_file ${fasta_file}.msasa.local.experiments.log \
        --extend \
        --match_score 5.0 \
        --mismatch_score 1.0 \
        --gap_score -2.0 \
        --temp 0.1 \
        --cooling_rate 0.992 \
        --min_temp 0.00001 \
        --max_no_changes 7500 \
        --changes 50 \
        --iteration_neighbors 10 \
        --quality_function local \
        --experiments ${params.experiments}
    """
}

process plotMumsaOverlap {
    label 'msasa'
    cache false
    tag "${file_path.simpleName}"
    publishDir params.results_dir + "/plots/${file_path.simpleName}", overwrite: true, mode: 'link'

    input:
        path(file_path)
    output:
        path("${file_path.simpleName}_plot.png")

    script:
    """
    python ${params.mumsa_script} ${file_path} ${file_path.simpleName}_plot.png
    """
}