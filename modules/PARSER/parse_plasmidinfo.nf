process parse_plasmidinfo {
    publishDir "${params.output}/", mode: 'copy'
    label 'ubuntu'   
  input:
    tuple val(name), path(results_abricate), path(results_IS)
  output:
	  tuple val(name), file("${name}_annotation_file.txt") 
  script:
    """
    # parse results into annotation_file
    grep -v "#FILE" ${results_abricate} | grep -w "ncbi" | cut -f2,6,3,4 | sed \$'s/\$/\\\tresistance_genes/' > annotation.tmp
    grep -v "#FILE" ${results_abricate} | grep -w "plasmidfinder" | cut -f2,6,3,4 | sed \$'s/\$/\\\tplasmid_genes/' >> annotation.tmp
    grep -v "#FILE" ${results_IS} | cut -f2,6,3,4 | sed \$'s/\$/\\\tInsertion_sequences/' >> annotation.tmp

    awk '{print \$4,\$1,\$2,\$3,\$5}' OFS='\\t' annotation.tmp > ${name}_annotation_file.txt
    """
}


    /* 
        //parser
        cat *_plasmid_*.abricate | grep -v "#FILE"
            grep "ncbi" -> get gene name (6) and start (3) and stop (4)
            grep "plasmidfinder" -> get gene name (6) and start (3) and stop (4)
        cat *_plasmid_*.tab | grep -v "#FILE"
            get gene name (6) and start (3) and stop (4)
        // samtools parser to get contig names and range
            parse away the "only positive" contigs from the other parser

        // chromo-map -> you should adjust this based on the "n of input" to modify die color vector
            // if it fails do the one with only one input lable
    */