process dev {
      publishDir "DEV_WORKFLOW_OUT/${repeater}", mode: 'copy', pattern: "*.cf"
      label 'dev'
      //errorStrategy 'ignore'
    input:
      file(database)
      //each repeater
    output:
      file("*.cf")
    shell:
      """
      tar xzf ${database}
      ls
      mv */* .
      ls
      centrifuge-build -p 7 --seed 42 --conversion-table ex.conv --taxonomy-tree ex.tree --name-table ex.name ex.fa ex
      """
}

/*
cpus = 40 ; memory = '320 GB'
with -p 8
builds up a Database with 30k genomes
*/


