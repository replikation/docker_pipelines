process prokka {
      publishDir "${params.output}/${name}/${params.assemblydir}", mode: 'copy', pattern: "${name}.gff"
      publishDir "${params.output}/${name}/${params.assemblydir}", mode: 'copy', pattern: "${name}.gbk"
      label 'prokka'
    input:
      tuple val(name), file(assembly) 
    output:
      tuple val(name), file("${name}.gff"), file("${name}.gbk")
    script:
      """
    	prokka --cpus ${task.cpus} --outdir output --prefix annotation ${assembly}
      mv output/annotation.gff ${name}.gff
      mv output/annotation.gbk ${name}.gbk
      """
}