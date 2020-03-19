process makeblastdatabase {
    label 'blast'   
  input:
    tuple val(name), path(fasta) 
  output:
	  tuple val(name), path("blast_database", type: 'dir') 
  script:
    """
  	makeblastdb -in ${fasta} -dbtype nucl -parse_seqids -out blast_database -title ${name}_db
    """
}

process blastn {
    label 'blast'   
  input:
    tuple val(name), path(fasta), path(database)
    each method
  output:
	  tuple val(name), path("*.tab") 
  script:
    """
  	blastn -query ${fasta} -db ${database} -out ${name}.blast -outfmt "6 sseqid qlen length sstart send" -num_threads ${task.cpus} -evalue 10E-80
    """
}
