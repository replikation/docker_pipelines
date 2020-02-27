process overview_parser {
    publishDir "${params.output}/", mode: 'copy'
    label 'ubuntu'   
  input:
    path(results_abri_chromosomes)
    path(results_abri_plasmids)
    path(results_abri_unclassified)
	  path(results_far_chromosomes)
    path(results_far_plasmids)
    path(results_far_unclassified)
    val(sampleIDs)
  output:
	  file("input.csv") 
  shell:
    """


    all_sample_IDs=\$(echo "${sampleIDs}" | tr -d " []" | tr "," "\\n")

    header=\$(echo "${sampleIDs}" | tr -d " []")
    printf "type,\${header}\\n" > input.csv

    ABRIMETHODS=\$(cat *.abricate | grep -v "^#FILE" | cut -f15 | sort | uniq )
    FARMETHODS=\$(head -n 2 -q *.fargene | grep "The used HMM-model was:" | awk '{print \$5}' | sort | uniq | sort )

    ##############
    # ABRICATE
    ##############

    while IFS= read -r method ; do
      # chromosome
      printf "\${method}-genome" >> input.csv

      while IFS= read -r sampleID ; do
        amount=\$(cat *_chromosome_\${sampleID}.abricate *_unclassified_\${sampleID}.abricate | grep "\${method}" | wc -l)
  
        if test \$amount -gt 0 
        then 
          printf ",\$amount" >> input.csv 
        else  
          printf ",NA" >> input.csv 
        fi

      done < <(echo "\${all_sample_IDs}")

      printf "\\n" >> input.csv

      # plasmid
      printf "\${method}-plasmid" >> input.csv
      while IFS= read -r sampleID ; do
        amount=\$(grep -c "\${method}" *_plasmid_\${sampleID}.abricate | wc -l)
  
        if test \$amount -gt 0  
        then 
          printf ",\$amount" >> input.csv 
        else  
          printf ",NA" >> input.csv 
        fi

      done < <(echo "\${all_sample_IDs}")
      printf "\\n" >> input.csv




    done < <(echo "\${ABRIMETHODS}")  

    ##############
    # FARGENE
    ##############

    """
}

/*


      # unclassified
      printf "\${method}-uncl" >> input.csv
      while IFS= read -r sampleID ; do
        amount=\$(grep -c "\${method}" *_unclassified_\${sampleID}.abricate | wc -l)
        if test \$amount -gt 0 
        then 
          printf ",\$amount" >> input.csv 
        else  
          printf ",NA" >> input.csv 
        fi

      done < <(echo "\${all_sample_IDs}")
      printf "\\n" >> input.csv

*/