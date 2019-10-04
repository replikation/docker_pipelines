process dev {
      publishDir "DEV_WORKFLOW_OUT/", mode: 'copy', pattern: "*.cf"
      label 'dev'
    input:
      file(database)
    output:
      set val(name), file("*.cf")
    shell:
      """
      tar xzf ${database}
      ls
      mv *centr*/* .
      ls
      centrifuge-build -p 4 --seed 42 --conversion-table ex.conv --taxonomy-tree ex.tree --name-table ex.name ex.fa ex
      """
}




