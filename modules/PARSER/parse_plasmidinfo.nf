process parse_plasmidinfo {
    publishDir "${params.output}/", mode: 'copy'
    label 'ubuntu'   
  input:
    tuple val(name), path(results_abricate), path(results_IS), path(results_fargene)
  output:
	  tuple val(name), file("${name}_annotation_file.txt") 
  script:
    """
    # parse results into annotation_file
    grep -v "#FILE" ${results_abricate} | grep -w "ncbi" | cut -f2,6,3,4 | sed \$'s/\$/\\\tresistance_genes/' > annotation.tmp
    grep -v "#FILE" ${results_abricate} | grep -w "plasmidfinder" | cut -f2,6,3,4 | sed \$'s/\$/\\\tplasmid_genes/' >> annotation.tmp
    grep -v "#FILE" ${results_IS} | cut -f2,6,3,4 | sed \$'s/\$/\\\tInsertion_sequences/' >> annotation.tmp


    # format to tab delimited and add uniq numbers behind each gene
    awk '{print \$4,\$1,\$2,\$3,\$5}' OFS='\\t' annotation.tmp |\
      awk '{printf "%d\\t%s\\n", NR, \$0}' |\
      awk '{print \$2"_"\$1"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6}' |\
      tr -d "'" \
      > ${name}_annotation_file.txt

    # add fargene

    cat ${results_fargene} | awk '{print \$4"\\t"\$1"\\t"\$16"\\t"\$17"\\tfargene"}' | sed 's|_seq._.||g' |\
      awk '{printf "%d\\t%s\\n", NR, \$0}' |\
      awk '{print \$2"_"\$1"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6}' |\
      tr -d "'" \
      >> ${name}_annotation_file.txt
    """
}


/* 
fargene infos

cat hmmsearchresults/retrieved-genes-class_A-hmmsearched.out | grep -v "^#"



*/



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