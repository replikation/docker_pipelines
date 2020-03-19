
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



/*
process minimap2_to_bam introduces a Math.random() function to add a random long number to the output.
This is introduced as a quick fix, because otherweise we create to much files of the same name which will be a issue in
the process later on when we use all of them together (differential binning)
*/