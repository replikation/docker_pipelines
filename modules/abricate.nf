process abricate {
    publishDir "${params.output}/${name}/ABR-Screening", mode: 'copy', pattern: "*.tab"
    label 'abricate'
    //if (workflow.profile == 'standard' && params.fastq) { maxForks 1 }   
  input:
    tuple val(name), path(fasta) 
    each method
  output:
	  tuple val(name), val("${method}"), path("*.tab")
  script:
    """
  	abricate ${fasta} --nopath --quiet --mincov 85 --db ${method} > ${method}.tab
    """
}

process abricate_compare {
    label 'abricate'
    //def random = (Math.random() + Math.random()).toString().md5().toString()
  input:
    tuple val(name), val(splitname), val(type), path(fasta) 
    each method
  output:
	  tuple val(name), path("${method}_${type}_${name}_*.abricate")
  script:
    """
  	abricate ${fasta} --nopath --quiet --mincov 85 --db ${method} > ${method}_${type}_${name}_\${PWD##*/}.abricate
    """
}

process abricate_transposon {
    label 'abricate'
    publishDir "${params.output}/${name}", mode: 'copy', pattern: "*.tab"
  input:
    tuple val(name), path(fasta) 
    path(database)
  output:
	  tuple val(name), path("*.tab")
  script:
    """
    mkdir -p db/mobile_elements
    cp ${database} db/mobile_elements/sequences
    abricate --setupdb --datadir \$PWD/db

  	abricate ${fasta} --datadir \$PWD/db --nopath --quiet --mincov 85 --db mobile_elements > mobile_elements_\${PWD##*/}.tab
    """
}

/*

    % cd /path/to/abricate/db     # this is the --datadir default option
    % mkdir tinyamr
    % cp /path/to/tinyamr.fa sequences
    % head -n 1 sequences
    >tinyamr~~~GENE_ID~~~GENE_ACC~~RESISTANCES some description here
    % abricate --setupdb
    % # or just do this: makeblastdb -in sequences -title tinyamr -dbtype nucl -hash_index


    % abricate --db tinyamr screen_this.fasta
  */