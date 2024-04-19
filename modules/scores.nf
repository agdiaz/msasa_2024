// https://github.com/gcollet/MstatX
// https://www.mdpi.com/2079-3197/11/11/212
// https://tcoffee.readthedocs.io/en/latest/tcoffee_main_documentation.html#evaluating-your-alignment
// https://teacheng.illinois.edu/LectureNotes/SequenceAlignment2_SubstitutionMatrix.pdf

process computeMumsaOverlapScore {
    tag "${referenceAlignment.simpleName}"

    publishDir params.results_dir + "/scores/${msasaCoincidencesAlignment.simpleName}", overwrite: true, mode: 'copy'

    input:
    path referenceAlignment
    path muscleAlignment
    path clustalAlignment
    path mafftAlignment
    path tCoffeeAlignment
    path kalignAlignment
    path msasaCoincidencesAlignment
    path msasaIdentityAlignment
    path msasaSimilarityBlosum62Alignment
    path msasaSimilarityPam250Alignment
    path msasaSimilarityGonnet92
    path msasaGlobalAlignment
    path msasaLocalAlignment

    output:
        path "${referenceAlignment.name}.mumsa.txt"

    script:
    """
    /software/mumsa-1.0/mumsa -r ${referenceAlignment} \
        ${muscleAlignment} \
        ${clustalAlignment} \
        ${mafftAlignment} \
        ${tCoffeeAlignment} \
        ${tCoffeeAlignment} \
        ${kalignAlignment} \
        ${msasaCoincidencesAlignment} \
        ${msasaIdentityAlignment} \
        ${msasaSimilarityBlosum62Alignment} \
        ${msasaSimilarityPam250Alignment} \
        ${msasaSimilarityGonnet92} \
        ${msasaGlobalAlignment} \
        ${msasaLocalAlignment} > ${referenceAlignment.name}.mumsa.txt
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
        path "${fasta_file.simpleName}_coreindex_${predictorName}.*"

    script:
    """
    t_coffee -infile=${fasta_file} -output=html -score -outfile=${fasta_file.simpleName}_coreindex_${predictorName}.html
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
        path "${fasta_file.simpleName}_tcs_${predictorName}.html"

    script:
    """
    t_coffee -infile ${fasta_file} -evaluate -output=score_html -outfile=${fasta_file.simpleName}_tcs_${predictorName}.html
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
