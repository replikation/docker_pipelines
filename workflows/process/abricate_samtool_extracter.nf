process abricate_samtool_extracter {
      label 'samtools'
    input:
      tuple path(tabfile), path(fasta) 
    output:
      path("*.fa")
    script:
      """
        samtools_clinker_parser.sh ${tabfile} ${params.searchterm} ${params.range}

        samtools faidx ${fasta}

        while read fastaheader; do    
                samtools faidx ${fasta} "\$fastaheader" >> "\${fastaheader}".fa
        done < seq_subsections.txt
      """
}