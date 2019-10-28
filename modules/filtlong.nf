process filtlong {
    //publishDir "${params.output}/${name}/ABR-Screening", mode: 'copy', pattern: "*.tab"
    label 'filtlong'
  input:
    tuple val(name), file(reads) 
  output:
	  tuple val(name), file("${name}_filtered.fastq") 
  script:
    """
  	filtlong --min_length 2000 ${reads} > ${name}_filtered.fastq
    """
}