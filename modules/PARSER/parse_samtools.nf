process parse_samtools {
    publishDir "${params.output}/", mode: 'copy'
    label 'samtools'   
  input:
    tuple val(name), file(annotation), file(fasta)
  output:
	  tuple val(name), file("chromosome_file.txt"), file(annotation)
  shell:
    """
    contiglist=\$(< ${annotation} cut -f 2 | sort | uniq )


    # parse contig stats into chromosome_file
    samtools faidx ${fasta}
    
    ## only using annotated contigs
    while IFS= read -r contig ; do
      grep "\${contig}" ${fasta}.fai >> index.list
    done < <(echo "\${contiglist}")

    cut -f1,2 index.list |  sed \$'s/\\\\t/\\\\t1\\\\t/' | sort -h -r -k 3  > chromosome_file.txt
    """
}

