process krakentools {
        label 'krakentools'
        publishDir "${params.output}/${name}/Read_classification/alpha_diversity", mode: 'copy', pattern: "${name}_alpha-diversity.txt"
        publishDir "${params.output}/${name}/Read_classification", mode: 'copy', pattern: "${name}.b.krona.txt"
    input:
        tuple val(name), path(brackenout), path(breport)
  	output:
    	tuple val(name), path("${name}_alpha-diversity.txt"), path("${name}.b.krona.txt")
  	script:
    """
    alpha_diversity.py -f ${brackenout} -a BP > ${name}_alpha-diversity.txt
    alpha_diversity.py -f ${brackenout} -a Sh >> ${name}_alpha-diversity.txt
    alpha_diversity.py -f ${brackenout} -a F  >> ${name}_alpha-diversity.txt
    alpha_diversity.py -f ${brackenout} -a Si  >> ${name}_alpha-diversity.txt
    alpha_diversity.py -f ${brackenout} -a ISi >> ${name}_alpha-diversity.txt

    # krona report
    kreport2krona.py -r ${breport} -o ${name}.b.krona.txt --no-intermediate-ranks 

    """
    stub:
    """
    touch ${name}_alpha-diversity.txt ${name}.b.krona.txt
    """
  }