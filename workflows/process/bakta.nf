process bakta {
      publishDir "${params.output}/${name}/bakta", mode: 'copy', pattern: "${name}-results"
      label 'bakta'
    input:
      tuple val(name), path(fasta) 
      file(database) 
    output:
      tuple val(name), file("${name}-results")
    script:
      """
      tar xzf ${database}
      rm ${database}
      export BAKTA_DB=./db
      amrfinder_update --force_update --database db/amrfinderplus-db/

      bakta --output ${name}-results --threads ${task.cpus} ${fasta}

      # reduce fingerprint on local systems
      rm -rf db
      """
}