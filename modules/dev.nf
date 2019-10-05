process dev {
      publishDir "DEV_WORKFLOW_OUT/${repeater}", mode: 'copy', pattern: "*.cf"
      label 'dev'
      errorStrategy 'ignore'
    input:
      file(database)
      each repeater
    output:
      set val(name), file("*.cf")
    shell:
      """
      tar xzf ${database}
      ls
      mv *centr*/* .
      ls
      centrifuge-build -p ${repeater} --seed 42 --conversion-table ex.conv --taxonomy-tree ex.tree --name-table ex.name ex.fa ex
      """
}




