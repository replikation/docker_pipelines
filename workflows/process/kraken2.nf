process kraken2_illumina_pe {
        label 'kraken2'
        publishDir "${params.output}/${name}/Read_classification", mode: 'copy'
    input:
        tuple val(name), path(reads)
        path(database)
  	output:
    	tuple val(name), path("${name}.kraken.out"), path("${name}.kreport")
  	script:
    """
    mkdir -p kraken_db && tar xzf ${database} -C kraken_db 


    kraken2 --db kraken_db --threads ${task.cpus} --paired --output ${name}.kraken.out --report ${name}.kreport ${reads}

    #  kraken has the opertunity to emit also unclassified reads!
    #  kraken2 --paired --classified-out cseqs#.fq seqs_1.fq seqs_2.fq

    # cleanup to reduce footprint
    rm -rf kraken_db/
    """
    stub:
    """
    touch ${name}.kraken.out ${name}.kreport
    """
  }