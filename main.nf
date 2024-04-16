// nextflow.enable.dsl=2

// Define parameters
params.input_dir = "$baseDir/datasets/BB11004"  // default input directory
params.msasa_script = "$baseDir/src/msasa_cli.py"
params.results_dir = "$baseDir/msasa_predictions_results"

params.match_score = 1
params.mismatch_score = 0
params.gap_score = -1

include {
    computeMumsaOverlapScore;

    computeCoreIndex as computeReferenceCoreIndex;
    computeCoreIndex as computeCoincidencesCoreIndex;
    computeCoreIndex as computeIdentityCoreIndex;
    computeCoreIndex as computeSimilarityBlosum62CoreIndex;
    computeCoreIndex as computeSimilarityPam250CoreIndex;
    computeCoreIndex as computeGlobalCoreIndex;
    computeCoreIndex as computeLocalCoreIndex;

    computeTransitiveConsistencyScore as computeReferenceTransitiveConsistencyScore;
    computeTransitiveConsistencyScore as computeCoincidencesTransitiveConsistencyScore;
    computeTransitiveConsistencyScore as computeIdentityTransitiveConsistencyScore;
    computeTransitiveConsistencyScore as computeSimilarityBlosum62TransitiveConsistencyScore;
    computeTransitiveConsistencyScore as computeSimilarityPam250TransitiveConsistencyScore;
    computeTransitiveConsistencyScore as computeGlobalTransitiveConsistencyScore;
    computeTransitiveConsistencyScore as computeLocalTransitiveConsistencyScore;

    computeBaliScore as computeCoincidencesBaliScore;
    computeBaliScore as computeIdentityBaliScore;
    computeBaliScore as computeSimilarityBlosum62BaliScore;
    computeBaliScore as computeSimilarityPam250BaliScore;
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

    alignCoincidences(fasta_files)
    alignIdentity(fasta_files)
    alignSimilarityBlosum62(fasta_files)
    alignSimilarityPam250(fasta_files)
    alignGlobal(fasta_files)
    alignLocal(fasta_files)

    computeMumsaOverlapScore(
        reference_fasta_files,
        muscle_fasta_files,
        clustal_fasta_files,
        mafft_fasta_files,
        t_coffee_fasta_files,
        kaling_fasta_files,
        alignCoincidences.out.alnFasta,
        alignIdentity.out.alnFasta,
        alignSimilarityBlosum62.out.alnFasta,
        alignSimilarityPam250.out.alnFasta,
        alignGlobal.out.alnFasta,
        alignLocal.out.alnFasta
    )

    plotMumsaOverlap(computeMumsaOverlapScore.out)

    // computeReferenceCoreIndex(reference_fasta_files, "reference")
    // computeCoincidencesCoreIndex(alignCoincidences.out.alnFasta, "coincidences")
    // computeIdentityCoreIndex(alignIdentity.out.alnFasta, "identity")
    // computeSimilarityBlosum62CoreIndex(alignSimilarityBlosum62.out.alnFasta, "similarity_blosum62")
    // computeSimilarityPam250CoreIndex(alignSimilarityPam250.out.alnFasta, "similarity_pam250")
    // computeGlobalCoreIndex(alignGlobal.out.alnFasta, "global")
    // computeLocalCoreIndex(alignLocal.out.alnFasta, "local")

    // computeReferenceTransitiveConsistencyScore(reference_fasta_files, "reference")
    // computeCoincidencesTransitiveConsistencyScore(alignCoincidences.out.alnFasta, "coincidences")
    // computeIdentityTransitiveConsistencyScore(alignIdentity.out.alnFasta, "identity")
    // computeSimilarityBlosum62TransitiveConsistencyScore(alignSimilarityBlosum62.out.alnFasta, "similarity_blosum62")
    // computeSimilarityPam250TransitiveConsistencyScore(alignSimilarityPam250.out.alnFasta, "similarity_pam250")
    // computeGlobalTransitiveConsistencyScore(alignGlobal.out.alnFasta, "global")
    // computeLocalTransitiveConsistencyScore(alignLocal.out.alnFasta, "local")

    // computeCoincidencesBaliScore(bali_base_xml_files, alignCoincidences.out.alnFasta, "coincidences")
    // computeIdentityBaliScore(bali_base_xml_files, alignIdentity.out.alnFasta, "identity")
    // computeSimilarityBlosum62BaliScore(bali_base_xml_files, alignSimilarityBlosum62.out.alnFasta, "similarity_blosum62")
    // computeSimilarityPam250BaliScore(bali_base_xml_files, alignSimilarityPam250.out.alnFasta, "similarity_pam250")
    // computeGlobalBaliScore(bali_base_xml_files, alignGlobal.out.alnFasta, "global")
    // computeLocalBaliScore(bali_base_xml_files, alignLocal.out.alnFasta, "local")
}

process alignCoincidences {
    // cache false
    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir params.results_dir + "/predictions/${fasta_file.simpleName}", overwrite: true, mode: 'link'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.coincidences.aln"), emit: alnFasta
        path("${fasta_file}.msasa.coincidences.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.coincidences.log"), emit: log
        path("${fasta_file}.msasa.coincidences.log.png"), emit: plot

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.coincidences.aln \
        --log_file ${fasta_file}.msasa.coincidences.log \
        --extend \
        --match_score 10 \
        --mismatch_score 0.0 \
        --gap_score 0 \
        --temp 0.1 \
        --cooling_rate 0.995 \
        --min_temp 0.000001 \
        --max_no_changes 7500 \
        --changes 5 \
        --iteration_neighbors 50 \
        --quality_function coincidences
    """
}

process alignIdentity {
    // cache false
    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir params.results_dir + "/predictions/${fasta_file.simpleName}", overwrite: true, mode: 'link'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.identity.aln"), emit: alnFasta
        path("${fasta_file}.msasa.identity.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.identity.log"), emit: log
        path("${fasta_file}.msasa.identity.log.png"), emit: plot

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.identity.aln \
        --log_file ${fasta_file}.msasa.identity.log \
        --extend \
        --match_score 10 \
        --mismatch_score 0 \
        --gap_score 0\
        --temp 0.1 \
        --cooling_rate 0.995 \
        --min_temp 0.000001 \
        --max_no_changes 7500 \
        --changes 5 \
        --iteration_neighbors 50 \
        --quality_function identity
    """
}

process alignSimilarityBlosum62 {
    // cache false
    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir params.results_dir + "/predictions/${fasta_file.simpleName}", overwrite: true, mode: 'link'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.similarity_blosum62.aln"), emit: alnFasta
        path("${fasta_file}.msasa.similarity_blosum62.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.similarity_blosum62.log"), emit: log
        path("${fasta_file}.msasa.similarity_blosum62.log.png"), emit: plot

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.similarity_blosum62.aln \
        --log_file ${fasta_file}.msasa.similarity_blosum62.log \
        --extend \
        --gap_score -8.0 \
        --temp 0.50 \
        --cooling_rate 0.995 \
        --min_temp 0.00001 \
        --max_no_changes 2500 \
        --changes 20 \
        --iteration_neighbors 25 \
        --quality_function similarity_blosum62
    """
}

process alignSimilarityPam250 {
    // cache false
    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir params.results_dir + "/predictions/${fasta_file.simpleName}", overwrite: true, mode: 'link'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.similarity_pam250.aln"), emit: alnFasta
        path("${fasta_file}.msasa.similarity_pam250.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.similarity_pam250.log"), emit: log
        path("${fasta_file}.msasa.similarity_pam250.log.png"), emit: plot

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.similarity_pam250.aln \
        --log_file ${fasta_file}.msasa.similarity_pam250.log \
        --extend \
        --gap_score -2.0 \
        --temp 1.0 \
        --cooling_rate 0.999 \
        --min_temp 0.00001 \
        --max_no_changes 2500 \
        --changes 10 \
        --iteration_neighbors 25 \
        --quality_function similarity_pam250
    """
}

process alignGlobal {
    cache false
    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir params.results_dir + "/predictions/${fasta_file.simpleName}", overwrite: true, mode: 'link'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.global.aln"), emit: alnFasta
        path("${fasta_file}.msasa.global.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.global.log"), emit: log
        path("${fasta_file}.msasa.global.log.png"), emit: plot

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.global.aln \
        --log_file ${fasta_file}.msasa.global.log \
        --extend \
        --match_score 5.0 \
        --mismatch_score 0.0 \
        --gap_score 0.0 \
        --temp 0.1 \
        --cooling_rate 0.992 \
        --min_temp 0.00001 \
        --max_no_changes 5000 \
        --changes 25 \
        --iteration_neighbors 50 \
        --quality_function global
    """
}

process alignLocal {
    cache false
    label 'msasa'
    tag "${fasta_file.simpleName}"
    publishDir params.results_dir + "/predictions/${fasta_file.simpleName}", overwrite: true, mode: 'link'

    input:
        path(fasta_file)

    output:
        path("${fasta_file}.msasa.local.aln"), emit: alnFasta
        path("${fasta_file}.msasa.local.aln.clustal"), emit: alnClustal
        path("${fasta_file}.msasa.local.log"), emit: log
        path("${fasta_file}.msasa.local.log.png"), emit: plot

    script:
    """
    python -W ignore ${params.msasa_script} ${fasta_file} ${fasta_file}.msasa.local.aln \
        --log_file ${fasta_file}.msasa.local.log \
        --extend \
        --match_score 10.0 \
        --mismatch_score 1.0 \
        --gap_score 0.0 \
        --temp 0.1 \
        --cooling_rate 0.995 \
        --min_temp 0.00001 \
        --max_no_changes 5000 \
        --changes 25 \
        --iteration_neighbors 50 \
        --quality_function local
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
    #!/Users/adrian/miniconda3/envs/msasa_2024/bin/python

    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt

    with open('${file_path}', 'r') as file:
        lines = file.readlines()[1:]  # Skip the header line

    # Parsing the data into lists
    alignments = []
    scores = []

    for line in lines:
        parts = line.rsplit(maxsplit=1)
        alignments.append(parts[0].strip())
        scores.append(float(parts[1]))

    # Plotting the data
    plt.figure(figsize=(10, 8))
    bars = plt.barh(alignments, scores, color='skyblue')
    plt.xlabel('Overlap Score')
    plt.title('Overlap Scores to the reference alignment')
    plt.gca().invert_yaxis()
    plt.grid(axis='x', linestyle='--', alpha=0.6)

    for bar, alignment in zip(bars, alignments):
        if '.msasa.' in alignment:
            # Highlight label background if it contains '.msasa.'
            plt.gca().text(bar.get_width(), bar.get_y() + bar.get_height()/2, f'{bar.get_width():.2f}',
                        va='center', ha='left', color='black', fontsize=9, backgroundcolor='yellow')
        else:
            # Regular label without highlighting
            plt.gca().text(bar.get_width(), bar.get_y() + bar.get_height()/2, f'{bar.get_width():.2f}',
                        va='center', ha='left', color='black', fontsize=9)

    plt.tight_layout()
    plt.savefig('${file_path.simpleName}_plot.png', dpi=300)
    """
}