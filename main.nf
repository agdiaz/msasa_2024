// nextflow.enable.dsl=2

// Define parameters
params.input_dir = "$baseDir/datasets/BB11004"  // default input directory
params.msasa_script = "$baseDir/src/msasa_cli.py"
params.mumsa_script = "$baseDir/src/mumsa_plot.py"
params.results_dir = "$baseDir/thesis_results"
params.experiments = 1

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

include {
    alignIdentity;
    alignCoincidences;
    alignSimilarityBlosum62;
    alignSimilarityPam250;
    alignGonnet92;
    alignGlobal;
    alignLocal;
} from "./modules/predictors.nf"

include {
    sortFasta as sortReferenceFasta;
    sortFasta as sortClustalFasta;
    sortFasta as sortKAlignFasta;
    sortFasta as sortMAFFTFasta;
    sortFasta as sortMuscleFasta;
    sortFasta as sortTCoffeeFasta;
} from "./modules/utils.nf"

workflow predict {
    take:
        fasta_files

    main:
        alignIdentity(fasta_files)
        alignCoincidences(fasta_files)
        alignSimilarityBlosum62(fasta_files)
        alignSimilarityPam250(fasta_files)
        alignGonnet92(fasta_files)
        alignGlobal(fasta_files)
        alignLocal(fasta_files)

    emit:
        coincidencesSample = alignCoincidences.out.alnFasta.flatten().randomSample(1)
        identitySample = alignIdentity.out.alnFasta.flatten().randomSample(1)
        blosum62Sample = alignSimilarityBlosum62.out.alnFasta.flatten().randomSample(1)
        pam250Sample = alignSimilarityPam250.out.alnFasta.flatten().randomSample(1)
        gonnet92Sample = alignGonnet92.out.alnFasta.flatten().randomSample(1)
        globalSample = alignGlobal.out.alnFasta.flatten().randomSample(1)
        localSample = alignLocal.out.alnFasta.flatten().randomSample(1)
}

workflow evaluateAlignments {
    take:
        bali_base_xml_files
        reference_fasta_files
        muscle_fasta_files
        clustal_fasta_files
        mafft_fasta_files
        t_coffee_fasta_files
        kaling_fasta_files
        coincidencesSample
        identitySample
        blosum62Sample
        pam250Sample
        gonnet92Sample
        globalSample
        localSample
    main:
        sortMuscleFasta(muscle_fasta_files, reference_fasta_files).set { sorted_muscle_fasta_files }
        sortClustalFasta(clustal_fasta_files, reference_fasta_files).set { sorted_clustal_fasta_files }
        sortMAFFTFasta(mafft_fasta_files, reference_fasta_files).set { sorted_mafft_fasta_files }
        sortTCoffeeFasta(t_coffee_fasta_files, reference_fasta_files).set { sorted_t_coffee_fasta_files }
        sortKAlignFasta(kaling_fasta_files, reference_fasta_files).set { sorted_kaling_fasta_files }

        computeCoincidencesCoreIndex(coincidencesSample, "coincidences")
        computeIdentityCoreIndex(identitySample, "identity")
        computeSimilarityBlosum62CoreIndex(blosum62Sample, "similarity_blosum62")
        computeSimilarityPam250CoreIndex(pam250Sample, "similarity_pam250")
        computeSimilarityGonnet92CoreIndex(gonnet92Sample, "similarity_pam250")
        computeGlobalCoreIndex(globalSample, "global")
        computeLocalCoreIndex(localSample, "local")

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

        computeMumsaOverlapScore(
            reference_fasta_files,
            sorted_muscle_fasta_files,
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
}

workflow {
    fasta_files = Channel.fromPath("${params.input_dir}/*.tfa", checkIfExists: true)
    reference_fasta_files = Channel.fromPath("${params.input_dir}/*.reference.fasta", checkIfExists: true)
    bali_base_xml_files = Channel.fromPath("${params.input_dir}/*.xml", checkIfExists: true)

    clustal_fasta_files = Channel.fromPath("${params.input_dir}/*_clustal.fasta", checkIfExists: true)
    kaling_fasta_files = Channel.fromPath("${params.input_dir}/*_kalign.fasta", checkIfExists: true)
    mafft_fasta_files = Channel.fromPath("${params.input_dir}/*_mafft.fasta", checkIfExists: true)
    muscle_fasta_files = Channel.fromPath("${params.input_dir}/*_muscle.fasta", checkIfExists: true)
    t_coffee_fasta_files = Channel.fromPath("${params.input_dir}/*_tCoffee.fasta", checkIfExists: true)

    predict(fasta_files)
    evaluateAlignments(
        bali_base_xml_files,
        reference_fasta_files,
        muscle_fasta_files,
        clustal_fasta_files,
        mafft_fasta_files,
        t_coffee_fasta_files,
        kaling_fasta_files,
        predict.out.coincidencesSample,
        predict.out.identitySample,
        predict.out.blosum62Sample,
        predict.out.pam250Sample,
        predict.out.gonnet92Sample,
        predict.out.globalSample,
        predict.out.localSample
    )

    computeReferenceCoreIndex(reference_fasta_files, "reference")
    computeReferenceTransitiveConsistencyScore(reference_fasta_files, "reference")
}



process plotMumsaOverlap {
    label 'msasa'
    cache false
    tag "${file_path.simpleName}"
    publishDir params.results_dir + "/plots/${file_path.simpleName}", overwrite: true, mode: 'copy'

    input:
        path(file_path)
    output:
        path("${file_path.simpleName}_mumsa_plot.png")

    script:
    """
    python ${params.mumsa_script} ${file_path} ${file_path.simpleName}_mumsa_plot.png
    """
}

workflow.onComplete {
    println "Pipeline for {params.input_dir} completed at: $workflow.complete"
    println "Command line: $workflow.commandLine"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}