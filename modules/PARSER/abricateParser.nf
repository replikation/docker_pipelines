process abricateParser {
    //publishDir "${params.output}/${name}/ABR-Screening", mode: 'copy', pattern: "*.csv"
    label 'ubuntu'   
  input:
    tuple val(name), file(results) 
  output:
	  tuple val(name), file("*.csv") 
  shell:
    """
    printf "amount;type;group\\n" > !{name}_ncbi.csv
  	# beta-lactamase
      blaData=\$(grep -w "ncbi" !{results} | cut -f5 | grep "bla" | cut -f1 -d "-"| sort | uniq -c |\
        sed -e 's/^[ \\t]*//' | tr " " ";" | sed -e 's/\$/;beta-lactamase/') 
      printf "\${blaData}\\n" >> !{name}_ncbi.csv
    # tetracycline
      tetData=\$(grep -w "ncbi" !{results} | grep "tetracycline" | cut -f5 | grep -v "bla" | sort |\
        uniq -c | sed -e 's/^[ \\t]*//'| tr -d "'" | tr " " ";" | sed -e 's/\$/;tetracycline-resistance/') 
      printf "\${tetData}\\n" >> !{name}_ncbi.csv
    # aminoglycoside
      aminoData=\$(grep -w "ncbi" !{results} | grep -v "efflux" | grep "aminoglycoside" | cut -f5 | grep -v "bla" | sort | awk '{print substr(\$0,0,4)}' |\
        uniq -c | sed -e 's/^[ \\t]*//'| tr -d "'" | tr " " ";" | sed -e 's/\$/;aminoglycoside-resistance/') 
      printf "\${aminoData}\\n" >> !{name}_ncbi.csv
    # efflux
      effluxData=\$(grep -w "ncbi" !{results} | grep -v "tetracycline" | grep "efflux" | cut -f5 | grep -v "bla" | sort | awk '{print substr(\$0,0,5)}' |\
        uniq -c | sed -e 's/^[ \\t]*//'| tr -d "'" | tr " " ";" | sed -e 's/\$/;efflux-system/') 
      printf "\${effluxData}\\n" >> !{name}_ncbi.csv
    # quinolones
      quinoData=\$(grep -w "ncbi" !{results} |  grep -v "efflux" | grep "quinolone" | cut -f5 | grep -v "bla" | sort | awk '{print substr(\$0,0,4)}' |\
        uniq -c | sed -e 's/^[ \\t]*//'| tr -d "'" | tr " " ";" | sed -e 's/\$/;quinolone-resistance/') 
      printf "\${quinoData}\\n" >> !{name}_ncbi.csv
    # other
      otherData=\$(grep -w "ncbi" !{results} |  grep -v "efflux" | grep -v "tetracycline" | grep -v "aminoglycoside" | grep -v "quinolone" |\
        cut -f5 | grep -v "bla" |  sort |\
        uniq -c | sed -e 's/^[ \\t]*//'| tr -d "'" | tr " " ";" | sed -e 's/\$/;other-resistance-genes/') 
      printf "\${otherData}\\n" >> !{name}_ncbi.csv
    """
}

// tempDir may be neccessary