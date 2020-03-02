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

    ABRIMETHODS=\$(cat *.abricate_* | grep -v "^#FILE" | cut -f15 | sort | uniq )
    FARMETHODS=\$(head -n 2 -q *.fargene_* | grep "The used HMM-model was:" | awk '{print \$5}' | sort | uniq | sort )

    ##############
    # ABRICATE
    ##############

    while IFS= read -r method ; do
      # chromosome
      printf "\${method}-genome" >> input.csv

      while IFS= read -r sampleID ; do
        amount=\$(cat *_chromosome_\${sampleID}.abricate_* *_unclassified_\${sampleID}.abricate_* | grep "\${method}" | wc -l)
  
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
        amount=\$(grep "\${method}" *_plasmid_\${sampleID}.abricate_* | wc -l)
  
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

    while IFS= read -r method ; do
      # chromosome
      printf "\${method}-genome" >> input.csv

      while IFS= read -r sampleID ; do
        amount=\$(cat *_chromosome_\${sampleID}.fargene_* *_unclassified_\${sampleID}.fargene_* |  grep "Number of predicted genes" |\
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
      printf "\${method}-plasmid" >> input.csv

      while IFS= read -r sampleID ; do
        amount=\$(cat *_plasmid_\${sampleID}.fargene_* |  grep "Number of predicted genes" |\
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