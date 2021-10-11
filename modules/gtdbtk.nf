// For a new release you have to change  the "dirname" within the gtdbtkgetdatabase.nf
// cpus are 1 due to RAM error of pplacer 

process gtdbtk {
      publishDir "${params.output}/${name}/gtdbtk", mode: 'copy', pattern: "${name}-results"
      label 'gtdbtk'
    input:
      tuple val(name), path(dir) 
      file(database) 
    output:
      tuple val(name), file("${name}-results")
    script:
      """
      tar xzf ${database}
      rm ${database}

      DBNAME=\$(echo release*/)
      export GTDBTK_DATA_PATH=./\${DBNAME}
      mkdir -p scratch_dir

      # file suffix get
      SUFFIXNAME=\$(ls ${dir} | head -1 | rev | cut -f1 -d "." | rev)

      gtdbtk classify_wf  --genome_dir ${dir} --cpus 1 --out_dir ${name}-results -x \${SUFFIXNAME} --scratch_dir scratch_dir
      """
}