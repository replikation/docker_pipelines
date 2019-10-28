// For a new release you have to change  the "dirname" within the gtdbtkgetdatabase.nf
// cpus are 1 due to RAM error of pplacer 

process gtdbtk {
      publishDir "${params.output}/${name}/gtdbtk", mode: 'copy', pattern: "${name}-results"
      label 'gtdbtk'
    input:
      tuple val(name), file(dir) 
      file(database) 
    output:
      tuple val(name), file("${name}-results")
    script:
      """
      tar xzf ${database}
      export GTDBTK_DATA_PATH=./release89
      gtdbtk classify_wf  --genome_dir ${dir} --cpus 1 --out_dir ${name}-results -x fa
      """
}