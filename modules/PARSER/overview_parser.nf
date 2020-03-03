process overview_parser {
    publishDir "${params.output}/", mode: 'copy'
    label 'ubuntu'   
  input:
    path(results)
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
        amount=\$(cat *_chromosome_\${sampleID}_*.abricate *_unclassified_\${sampleID}_*.abricate | grep "\${method}" | wc -l)
  
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
        amount=\$(grep "\${method}" *_plasmid_\${sampleID}_*.abricate | wc -l)
  
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
    # FARGENE                 ERROR FAR METHODE IS NOT INCLUDED via HEAD
    ##############

    while IFS= read -r method ; do
      # chromosome
      printf "\${method%.hmm}-genome" >> input.csv

      while IFS= read -r sampleID ; do
        amount=\$(cat *_chromosome_\${sampleID}_*.fargene *_unclassified_\${sampleID}_*.fargene | grep -A 5 "\${method}" |\
                grep "Number of predicted genes" |\
                awk '{printf "%s\\n",\$5}' | awk '{s+=\$1} END {print s}')

        if test \$amount -gt 0 
        then 
          printf ",\$amount" >> input.csv 
        else  
          printf ",NA" >> input.csv 
        fi

      done < <(echo "\${all_sample_IDs}")

      printf "\\n" >> input.csv

      # plasmid
      printf "\${method%.hmm}-plasmid" >> input.csv

      while IFS= read -r sampleID ; do
        amount=\$(cat *_plasmid_\${sampleID}_*.fargene | grep -A 5 "\${method}" |\
                grep "Number of predicted genes" |\
                awk '{printf "%s\\n",\$5}' | awk '{s+=\$1} END {print s}')
  
        if test \$amount -gt 0  
        then 
          printf ",\$amount" >> input.csv 
        else  
          printf ",NA" >> input.csv 
        fi

      done < <(echo "\${all_sample_IDs}")

      printf "\\n" >> input.csv

    done < <(echo "\${FARMETHODS}")  
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