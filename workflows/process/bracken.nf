process bracken {
        label 'bracken'
        publishDir "${params.output}/${name}/Read_classification", mode: 'copy'
    input:
        tuple val(name), path(krakenout), path(kreport)
        path(database)
  	output:
    	tuple val(name), path("${name}.bracken"), path("${name}.breport")
  	script:
    """
    mkdir -p kraken_db && tar xzf ${database} -C kraken_db

    bracken -d kraken_db -i ${name}.kreport -r ${params.readlength} -l S -t ${task.cpus} \
     -o ${name}.bracken -w ${name}.breport


    # cleanup to reduce footprint
    rm -rf kraken_db/
    """
    stub:
    """
    touch ${name}.bracken ${name}.breport
    """
  }