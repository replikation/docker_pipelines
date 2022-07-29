process prokka_gbk {
      label 'prokka'
    input:
      path(fastas) 
    output:
      path("*.gbk")
    script:
      """
        FILENAME=\$(basename -s .fa "${fastas}")
        prokka --cpus ${task.cpus} --outdir output --prefix annotation ${fastas} --centre X --compliant
        mv output/annotation.gbk "\${FILENAME}.gbk"
      """
}