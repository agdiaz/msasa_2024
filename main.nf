// nextflow.enable.dsl=2

// Define parameters
params.input_dir = "$baseDir/datasets/BB11004"  // default input directory
params.msasa_script = "$baseDir/src/msasa_cli.py"
params.mumsa_script = "$baseDir/src/mumsa_plot.py"
params.results_dir = "$baseDir/thesis_results"
params.experiments = 1

include {
    computeMumsaOverlapScore;

    // computeCoreIndex as computeReferenceCoreIndex;
    // computeCoreIndex as computeCoincidencesCoreIndex;
    // computeCoreIndex as computeIdentityCoreIndex;
    // computeCoreIndex as computeSimilarityBlosum62CoreIndex;
    // computeCoreIndex as computeSimilarityPam250CoreIndex;
    // computeCoreIndex as computeSimilarityGonnet92CoreIndex;
    // computeCoreIndex as computeGlobalCoreIndex;
    // computeCoreIndex as computeLocalCoreIndex;

    // computeTransitiveConsistencyScore as computeReferenceTransitiveConsistencyScore;
    // computeTransitiveConsistencyScore as computeCoincidencesTransitiveConsistencyScore;
    // computeTransitiveConsistencyScore as computeIdentityTransitiveConsistencyScore;
    // computeTransitiveConsistencyScore as computeSimilarityBlosum62TransitiveConsistencyScore;
    // computeTransitiveConsistencyScore as computeSimilarityPam250TransitiveConsistencyScore;
    // computeTransitiveConsistencyScore as computeSimilarityGonnet92TransitiveConsistencyScore;
    // computeTransitiveConsistencyScore as computeGlobalTransitiveConsistencyScore;
    // computeTransitiveConsistencyScore as computeLocalTransitiveConsistencyScore;

    computeBaliScore as computeCoincidencesBaliScoreFree;
    computeBaliScore as computeIdentityBaliScoreFree;
    computeBaliScore as computeSimilarityBlosum62BaliScoreFree;
    computeBaliScore as computeSimilarityPam250BaliScoreFree;
    computeBaliScore as computeSimilarityGonnet92BaliScoreFree;
    computeBaliScore as computeGlobalBaliScoreFree;
    computeBaliScore as computeLocalBaliScoreFree;

    computeBaliScore as computeCoincidencesBaliScoreStrict;
    computeBaliScore as computeIdentityBaliScoreStrict;
    computeBaliScore as computeSimilarityBlosum62BaliScoreStrict;
    computeBaliScore as computeSimilarityPam250BaliScoreStrict;
    computeBaliScore as computeSimilarityGonnet92BaliScoreStrict;
    computeBaliScore as computeGlobalBaliScoreStrict;
    computeBaliScore as computeLocalBaliScoreStrict;
} from "./modules/scores.nf"

include {
    alignIdentity as alignIdentityFree;
    alignIdentity as alignIdentityStrict;
    alignCoincidences as alignCoincidencesFree;
    alignCoincidences as alignCoincidencesStrict;
    alignSimilarityBlosum62 as alignSimilarityBlosum62Free;
    alignSimilarityBlosum62 as alignSimilarityBlosum62Strict;
    alignSimilarityPam250 as alignSimilarityPam250Free;
    alignSimilarityPam250 as alignSimilarityPam250Strict;
    alignGonnet92 as alignGonnet92Free;
    alignGonnet92 as alignGonnet92Strict;
    alignGlobal as alignGlobalFree;
    alignGlobal as alignGlobalStrict;
    alignLocal as alignLocalFree;
    alignLocal as alignLocalStrict;
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
        alignIdentityFree(fasta_files, false)
        alignIdentityStrict(fasta_files, true)
        alignCoincidencesFree(fasta_files, false)
        alignCoincidencesStrict(fasta_files, true)
        alignSimilarityBlosum62Free(fasta_files, false)
        alignSimilarityBlosum62Strict(fasta_files, true)
        alignSimilarityPam250Free(fasta_files, false)
        alignSimilarityPam250Strict(fasta_files, true)
        alignGonnet92Free(fasta_files, false)
        alignGonnet92Strict(fasta_files, true)
        alignGlobalFree(fasta_files, false)
        alignGlobalStrict(fasta_files, true)
        alignLocalFree(fasta_files, false)
        alignLocalStrict(fasta_files, true)

    emit:
        coincidencesSampleFree = alignCoincidencesFree.out.alnFasta.flatten().randomSample(1)
        coincidencesSampleStrict = alignCoincidencesStrict.out.alnFasta.flatten().randomSample(1)
        identitySampleFree = alignIdentityFree.out.alnFasta.flatten().randomSample(1)
        identitySampleStrict = alignIdentityStrict.out.alnFasta.flatten().randomSample(1)
        blosum62SampleFree = alignSimilarityBlosum62Free.out.alnFasta.flatten().randomSample(1)
        blosum62SampleStrict = alignSimilarityBlosum62Strict.out.alnFasta.flatten().randomSample(1)
        pam250SampleFree = alignSimilarityPam250Free.out.alnFasta.flatten().randomSample(1)
        pam250SampleStrict = alignSimilarityPam250Strict.out.alnFasta.flatten().randomSample(1)
        gonnet92SampleFree = alignGonnet92Free.out.alnFasta.flatten().randomSample(1)
        gonnet92SampleStrict = alignGonnet92Strict.out.alnFasta.flatten().randomSample(1)
        globalSampleFree = alignGlobalFree.out.alnFasta.flatten().randomSample(1)
        globalSampleStrict = alignGlobalStrict.out.alnFasta.flatten().randomSample(1)
        localSampleFree = alignLocalFree.out.alnFasta.flatten().randomSample(1)
        localSampleStrict = alignLocalStrict.out.alnFasta.flatten().randomSample(1)
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
        coincidencesSampleFree
        coincidencesSampleStrict
        identitySampleFree
        identitySampleStrict
        blosum62SampleFree
        blosum62SampleStrict
        pam250SampleFree
        pam250SampleStrict
        gonnet92SampleFree
        gonnet92SampleStrict
        globalSampleFree
        globalSampleStrict
        localSampleFree
        localSampleStrict

    main:
        sortMuscleFasta(muscle_fasta_files, reference_fasta_files).set { sorted_muscle_fasta_files }
        sortClustalFasta(clustal_fasta_files, reference_fasta_files).set { sorted_clustal_fasta_files }
        sortMAFFTFasta(mafft_fasta_files, reference_fasta_files).set { sorted_mafft_fasta_files }
        sortTCoffeeFasta(t_coffee_fasta_files, reference_fasta_files).set { sorted_t_coffee_fasta_files }
        sortKAlignFasta(kaling_fasta_files, reference_fasta_files).set { sorted_kaling_fasta_files }

        // computeCoincidencesCoreIndex(coincidencesSample, "coincidences")
        // computeIdentityCoreIndex(identitySample, "identity")
        // computeSimilarityBlosum62CoreIndex(blosum62Sample, "similarity_blosum62")
        // computeSimilarityPam250CoreIndex(pam250Sample, "similarity_pam250")
        // computeSimilarityGonnet92CoreIndex(gonnet92Sample, "similarity_pam250")
        // computeGlobalCoreIndex(globalSample, "global")
        // computeLocalCoreIndex(localSample, "local")

        // computeCoincidencesTransitiveConsistencyScore(coincidencesSample, "coincidences")
        // computeIdentityTransitiveConsistencyScore(identitySample, "identity")
        // computeSimilarityBlosum62TransitiveConsistencyScore(blosum62Sample, "similarity_blosum62")
        // computeSimilarityPam250TransitiveConsistencyScore(pam250Sample, "similarity_pam250")
        // computeSimilarityGonnet92TransitiveConsistencyScore(gonnet92Sample, "similarity_gonnet92")
        // computeGlobalTransitiveConsistencyScore(globalSample, "global")
        // computeLocalTransitiveConsistencyScore(localSample, "local")

        computeCoincidencesBaliScoreFree(bali_base_xml_files, coincidencesSampleFree, "coincidences_free")
        computeIdentityBaliScoreFree(bali_base_xml_files, identitySampleFree, "identity_free")
        computeSimilarityBlosum62BaliScoreFree(bali_base_xml_files, blosum62SampleFree, "similarity_blosum62_free")
        computeSimilarityPam250BaliScoreFree(bali_base_xml_files, pam250SampleFree, "similarity_pam250_free")
        computeSimilarityGonnet92BaliScoreFree(bali_base_xml_files, gonnet92SampleFree, "similarity_gonnet92_free")
        computeGlobalBaliScoreFree(bali_base_xml_files, globalSampleFree, "global_free")
        computeLocalBaliScoreFree(bali_base_xml_files, localSampleFree, "local_free")

        computeCoincidencesBaliScoreStrict(bali_base_xml_files, coincidencesSampleStrict, "coincidences_strict")
        computeIdentityBaliScoreStrict(bali_base_xml_files, identitySampleStrict, "identity_strict")
        computeSimilarityBlosum62BaliScoreStrict(bali_base_xml_files, blosum62SampleStrict, "similarity_blosum62_strict")
        computeSimilarityPam250BaliScoreStrict(bali_base_xml_files, pam250SampleStrict, "similarity_pam250_strict")
        computeSimilarityGonnet92BaliScoreStrict(bali_base_xml_files, gonnet92SampleStrict, "similarity_gonnet92_strict")
        computeGlobalBaliScoreStrict(bali_base_xml_files, globalSampleStrict, "global_strict")
        computeLocalBaliScoreStrict(bali_base_xml_files, localSampleStrict, "local_strict")

        computeMumsaOverlapScore(
            reference_fasta_files,

            sorted_muscle_fasta_files,
            clustal_fasta_files,
            mafft_fasta_files,
            t_coffee_fasta_files,
            kaling_fasta_files,

            coincidencesSampleFree,
            coincidencesSampleStrict,
            identitySampleFree,
            identitySampleStrict,
            blosum62SampleFree,
            blosum62SampleStrict,
            pam250SampleFree,
            pam250SampleStrict,
            gonnet92SampleFree,
            gonnet92SampleStrict,
            globalSampleFree,
            globalSampleStrict,
            localSampleFree,
            localSampleStrict
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
        predict.out.coincidencesSampleFree,
        predict.out.coincidencesSampleStrict,
        predict.out.identitySampleFree,
        predict.out.identitySampleStrict,
        predict.out.blosum62SampleFree,
        predict.out.blosum62SampleStrict,
        predict.out.pam250SampleFree,
        predict.out.pam250SampleStrict,
        predict.out.gonnet92SampleFree,
        predict.out.gonnet92SampleStrict,
        predict.out.globalSampleFree,
        predict.out.globalSampleStrict,
        predict.out.localSampleFree,
        predict.out.localSampleStrict
    )

    // computeReferenceCoreIndex(reference_fasta_files, "reference")
    // computeReferenceTransitiveConsistencyScore(reference_fasta_files, "reference")
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
    println "Pipeline for {params.input_dir} completed at: $workflow.complete ($workflow.duration)"
    println "Command line: $workflow.commandLine"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}