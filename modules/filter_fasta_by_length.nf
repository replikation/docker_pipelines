process filter_fasta_by_length {
    label 'ubuntu'
  input:
    tuple val(name), path(fasta)
  output:
	  tuple val(name), path("${name}_filtered.fasta")
  script:
    """
  	# make fasta files to one liner
    sed ':a;N;/^>/M!s/\n//;ta;P;D' ${fasta} |\
    awk '/^>/ { getline seq } length(seq) > 1500 { print \$0 "\\n" seq }' |\
    awk '/^>/ { getline seq } length(seq) < 1000000 { print \$0 "\\n" seq }' > ${name}_filtered.fasta 
    """
}
