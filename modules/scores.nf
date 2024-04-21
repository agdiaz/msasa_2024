// https://github.com/gcollet/MstatX
// https://www.mdpi.com/2079-3197/11/11/212
// https://tcoffee.readthedocs.io/en/latest/tcoffee_main_documentation.html#evaluating-your-alignment
// https://teacheng.illinois.edu/LectureNotes/SequenceAlignment2_SubstitutionMatrix.pdf

process computeMumsaOverlapScore {
    tag "${referenceAlignment.simpleName}"

    publishDir params.results_dir + "/scores/${referenceAlignment.simpleName}", overwrite: true, mode: 'copy'

    input:
        path referenceAlignment

        path muscleAlignment
        path clustalAlignment
        path mafftAlignment
        path tCoffeeAlignment
        path kalignAlignment

        path msasaCoincidencesAlignmentFree
        path msasaCoincidencesAlignmentStrict
        path msasaIdentityAlignmentFree
        path msasaIdentityAlignmentStrict
        path msasaSimilarityBlosum62AlignmentFree
        path msasaSimilarityBlosum62AlignmentStrict
        path msasaSimilarityPam250AlignmentFree
        path msasaSimilarityPam250AlignmentStrict
        path msasaSimilarityGonnet92AlignmentFree
        path msasaSimilarityGonnet92AlignmentStrict
        path msasaGlobalAlignmentFree
        path msasaGlobalAlignmentStrict
        path msasaLocalAlignmentFree
        path msasaLocalAlignmentStrict

    output:
        path "${referenceAlignment.simpleName}.mumsa.txt"

    script:
    """
    /software/mumsa-1.0/mumsa -r ${referenceAlignment} \
        ${muscleAlignment} \
        ${clustalAlignment} \
        ${mafftAlignment} \
        ${tCoffeeAlignment} \
        ${kalignAlignment} \
        ${msasaCoincidencesAlignmentFree} \
        ${msasaCoincidencesAlignmentStrict} \
        ${msasaIdentityAlignmentFree} \
        ${msasaIdentityAlignmentStrict} \
        ${msasaSimilarityBlosum62AlignmentFree} \
        ${msasaSimilarityBlosum62AlignmentStrict} \
        ${msasaSimilarityPam250AlignmentFree} \
        ${msasaSimilarityPam250AlignmentStrict} \
        ${msasaSimilarityGonnet92AlignmentFree} \
        ${msasaSimilarityGonnet92AlignmentStrict} \
        ${msasaGlobalAlignmentFree} \
        ${msasaGlobalAlignmentStrict} \
        ${msasaLocalAlignmentFree} \
        ${msasaLocalAlignmentStrict} > ${referenceAlignment.simpleName}.mumsa.txt
    """
}

// t_coffee -infile='$fasta_file' -output=score_ascii -score -outfile='${fasta_file.simpleName}_coreindex_${predictorName}.txt'
process computeCoreIndex {
    label 'msasa'
    errorStrategy 'ignore'

    tag "${fasta_file.simpleName}"
    publishDir params.results_dir + "/scores/${fasta_file.simpleName}", overwrite: true, mode: 'copy'

    input:
        path fasta_file
        val predictorName

    output:
        path "${fasta_file.simpleName}_coreindex_${predictorName}.txt"

    script:
    """
    t_coffee -infile=${fasta_file} -output=score_ascii -score -outfile=${fasta_file.simpleName}_coreindex_${predictorName}.txt
    """
}
// t_coffee -infile '$fasta_file' -evaluate -output=score_ascii -outfile='${fasta_file.simpleName}_tcs_${predictorName}.txt'
process computeTransitiveConsistencyScore {
    label 'msasa'
    errorStrategy 'ignore'
    tag "${fasta_file.simpleName}"
    publishDir params.results_dir + "/scores/${fasta_file.simpleName}", overwrite: true, mode: 'copy'

    input:
        path fasta_file
        val predictorName

    output:
        path "${fasta_file.simpleName}_tcs_${predictorName}.txt"

    script:
    """
    t_coffee -infile ${fasta_file} -evaluate -output=score_ascii -outfile=${fasta_file.simpleName}_tcs_${predictorName}.txt
    """
}

process computeBaliScore {
    tag "${fasta_file.simpleName} ${predictorName}"
    errorStrategy 'ignore'
    cache false

    publishDir params.results_dir + "/scores/${fasta_file.simpleName}", overwrite: true, mode: 'copy'

    input:
        path reference
        path fasta_file
        val predictorName

    output:
        path "${fasta_file.name}_baliscore_${predictorName}.txt", emit: score

    shell:
    """
    /software/bali-score/target/release/bali-score -t ${fasta_file} \
        -r ${reference} \
        -o ${fasta_file.name}_baliscore_${predictorName}.txt
    """
}
