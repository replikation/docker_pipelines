
process bedtools {
  label 'bedtools'
    input:
      tuple val(name) , file(bamfile)
    output:
      tuple val(name) , file("myseq.bedGraph")
    shell:
    """
    bedtools genomecov -bg -ibam ${bamfile} > myseq.bedGraph
    """
}


