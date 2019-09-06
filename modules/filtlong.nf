process filtlong {
    //publishDir "${params.output}/${name}/ABR-Screening", mode: 'copy', pattern: "*.tab"
    label 'filtlong'
  input:
    set val(name), file(reads) 
  output:
	  set val(name), file("${name}_filtered.fastq") 
  script:
    """
  	filtlong --min_length 2000 ${reads} > ${name}_filtered.fastq
    """
}