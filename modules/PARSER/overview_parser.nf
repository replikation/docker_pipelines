process overview_parser {
    publishDir "${params.output}/", mode: 'copy'
    label 'ubuntu'   
  input:
    path(results)
    val(sampleIDs)
  output:
	  file("input*.csv") 
  shell:
  if (params.coverage)
    """
    all_sample_IDs=\$(echo "${sampleIDs}" | tr -d " []" | tr "," "\\n")

    header=\$(echo "${sampleIDs}" | tr -d " []")
    printf "ID,type,type,method,\${header}\\n" > input.csv
    printf "ID,type,type,method,\${header}\\n" > input_coverage.csv

    ABRIMETHODS=\$(cat *.abricate | grep -v "^#FILE" | cut -f15 | sort | uniq )
    FARMETHODS=\$(head -n 2 -q *gencount.fargene | grep "The used HMM-model was:" | awk '{print \$5}' | sort | uniq | sort )

    ##############
    # ABRICATE
    ##############

    while IFS= read -r method ; do
      # chromosome
      printf "\${method}-genome,\${method},Genome,Antibiotic resistance" >> input.csv
      printf "\${method}-genome,\${method},Genome,Antibiotic resistance" >> input_coverage.csv

      while IFS= read -r sampleID ; do
        amount=\$(cat *_chromosome_\${sampleID}_*.abricate *_unclassified_\${sampleID}_*.abricate | grep "\${method}" | wc -l)
        coverage=\$(cat *_chromosome_\${sampleID}_*.abricate *_unclassified_\${sampleID}_*.abricate | grep "\${method}" | cut -f2 | rev | cut -f1 -d "_" | rev | awk '{s+=\$1} END {print s}')

        if test \$amount -gt 0 
        then 
          printf ",\$amount" >> input.csv 
        else  
          printf ",NA" >> input.csv 
        fi

          printf ",\$coverage" >> input_coverage.csv 

      done < <(echo "\${all_sample_IDs}")

      printf "\\n" >> input.csv
      printf "\\n" >> input_coverage.csv

      # plasmid
      printf "\${method}-plasmid,\${method},Plasmid,Antibiotic resistance" >> input.csv
      printf "\${method}-plasmid,\${method},Plasmid,Antibiotic resistance" >> input_coverage.csv

      while IFS= read -r sampleID ; do
        amount=\$(grep "\${method}" *_plasmid_\${sampleID}_*.abricate | wc -l)
        coverage=\$(cat *_plasmid_\${sampleID}_*.abricate | grep "\${method}" | cut -f2 | rev | cut -f1 -d "_" | rev | awk '{s+=\$1} END {print s}')
  
        if test \$amount -gt 0  
        then 
          printf ",\$amount" >> input.csv 
        else  
          printf ",NA" >> input.csv 
        fi

          printf ",\$coverage" >> input_coverage.csv 


      done < <(echo "\${all_sample_IDs}")

      printf "\\n" >> input.csv
      printf "\\n" >> input_coverage.csv

    done < <(echo "\${ABRIMETHODS}")  

    ##############
    # FARGENE 
    ##############

    while IFS= read -r method ; do
      # chromosome
      printf "\${method%.hmm}-genome,\${method%.hmm},Genome,beta-lactamase class" >> input.csv
      printf "\${method%.hmm}-genome,\${method%.hmm},Genome,beta-lactamase class" >> input_coverage.csv

      while IFS= read -r sampleID ; do
        amount=\$(cat *_chromosome_\${sampleID}_*.fargene *_unclassified_\${sampleID}_*.fargene | grep -A 5 "\${method}" |\
                grep "Number of predicted genes" |\
                awk '{printf "%s\\n",\$5}' | awk '{s+=\$1} END {print s}')

        coverage=\$(ls *_chromosome_\${sampleID}_*coverage.fargene *_unclassified_\${sampleID}_*coverage.fargene |\
            grep -i "\${method%.hmm}" | xargs cat | cut -f1 -d " " | sed 's|_seq._.||g' | rev | cut -f1 -d"_" | rev | awk '{s+=\$1} END {print s}')

        if test \$amount -gt 0 
        then 
          printf ",\$amount" >> input.csv 
        else  
          printf ",NA" >> input.csv 
        fi

        printf ",\$coverage" >> input_coverage.csv 

      done < <(echo "\${all_sample_IDs}")

      printf "\\n" >> input.csv
      printf "\\n" >> input_coverage.csv

      # plasmid
      printf "\${method%.hmm}-plasmid,\${method%.hmm},Plasmid,beta-lactamase class" >> input.csv
      printf "\${method%.hmm}-plasmid,\${method%.hmm},Plasmid,beta-lactamase class" >> input_coverage.csv

      while IFS= read -r sampleID ; do
        amount=\$(cat *_plasmid_\${sampleID}_*.fargene | grep -A 5 "\${method}" |\
                grep "Number of predicted genes" |\
                awk '{printf "%s\\n",\$5}' | awk '{s+=\$1} END {print s}')
  
      coverage=\$(ls *_plasmid_\${sampleID}_*coverage.fargene |\
          grep -i "\${method%.hmm}" | xargs cat | cut -f1 -d " " | sed 's|_seq._.||g' | rev | cut -f1 -d"_" | rev | awk '{s+=\$1} END {print s}')

       if test \$amount -gt 0  
        then 
          printf ",\$amount" >> input.csv 
        else  
          printf ",NA" >> input.csv 
        fi

        printf ",\$coverage" >> input_coverage.csv 

      done < <(echo "\${all_sample_IDs}")

      printf "\\n" >> input.csv
      printf "\\n" >> input_coverage.csv

    done < <(echo "\${FARMETHODS}")  
    """

else if (!params.coverage)
    """
    all_sample_IDs=\$(echo "${sampleIDs}" | tr -d " []" | tr "," "\\n")

    header=\$(echo "${sampleIDs}" | tr -d " []")
    printf "ID,type,type,method,\${header}\\n" > input.csv

    ABRIMETHODS=\$(cat *.abricate | grep -v "^#FILE" | cut -f15 | sort | uniq )
    FARMETHODS=\$(head -n 2 -q *.fargene | grep "The used HMM-model was:" | awk '{print \$5}' | sort | uniq | sort )

    ##############
    # ABRICATE
    ##############

    while IFS= read -r method ; do
      # chromosome
      printf "\${method}-genome,\${method},Genome,Antibiotic resistance" >> input.csv

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
      printf "\${method}-plasmid,\${method},Plasmid,Antibiotic resistance" >> input.csv

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
    # FARGENE 
    ##############

    while IFS= read -r method ; do
      # chromosome
      printf "\${method%.hmm}-genome,\${method%.hmm},Genome,beta-lactamase class" >> input.csv

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
      printf "\${method%.hmm}-plasmid,\${method%.hmm},Plasmid,beta-lactamase class" >> input.csv

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


ls *_chromosome_Ger-01_*coverage.fargene *_plasmid_Ger-01_*coverage.fargene | grep -iT "class_B_3" | xargs cat | cut -f1 -d " " | sed 's|_seq._.||g' | rev | cut -f1 -d"_" | rev | awk '{s+=$1} END {print s}'



*/