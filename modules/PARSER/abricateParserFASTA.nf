process abricateParserFASTA {
    //publishDir "${params.output}/${name}/ABR-Screening", mode: 'copy', pattern: "*.csv"
    label 'ubuntu'   
  input:
    set val(name), val(method), file(results) 
  output:
	  set val(name), val(method), file("*.csv") 
  shell:
    """
    printf "amount;type;group\\n" > !{name}_!{method}.csv
  	# beta-lactamase
      blaData=\$(tail -n+2 !{results} | grep -E "bla[A-Z]|lactamase" | cut -f5 | sort | uniq -c |\
        sed -e 's/^[ \\t]*//' | tr " " ";" | sed -e 's/\$/;beta-lactamase/') 
      printf "\${blaData}\\n" >> !{name}_!{method}.csv
    # tetracycline
      tetData=\$(tail -n+2 !{results} | grep "tetracycline" | cut -f5 | grep -v "bla" | sort |\
        uniq -c | sed -e 's/^[ \\t]*//'| tr -d "'" | tr " " ";" | sed -e 's/\$/;tetracycline-resistance/') 
      printf "\${tetData}\\n" >> !{name}_!{method}.csv
    # aminoglycoside
      aminoData=\$(tail -n+2 !{results} | grep -v "efflux" | grep "aminoglycoside" | cut -f5 | grep -v "bla" | sort |\
        uniq -c | sed -e 's/^[ \\t]*//'| tr -d "'" | tr " " ";" | sed -e 's/\$/;aminoglycoside-resistance/') 
      printf "\${aminoData}\\n" >> !{name}_!{method}.csv
    # efflux
      effluxData=\$(tail -n+2 !{results} | grep -v "tetracycline" | grep "efflux" | cut -f5 | grep -v "bla" | sort |\
        uniq -c | sed -e 's/^[ \\t]*//'| tr -d "'" | tr " " ";" | sed -e 's/\$/;efflux-system/') 
      printf "\${effluxData}\\n" >> !{name}_!{method}.csv
    # quinolones
      quinoData=\$(tail -n+2 !{results} |  grep -v "efflux" | grep "quinolone" | cut -f5 | grep -v "bla" | sort |\
        uniq -c | sed -e 's/^[ \\t]*//'| tr -d "'" | tr " " ";" | sed -e 's/\$/;quinolone-resistance/') 
      printf "\${quinoData}\\n" >> !{name}_!{method}.csv
    # other
      otherData=\$(tail -n+2 !{results} |  grep -v "efflux" | grep -v "tetracycline" | grep -v "aminoglycoside" | grep -v "quinolone" |\
        cut -f5 | grep -vE "bla[A-Z]|lactamase" |  sort |\
        uniq -c | sed -e 's/^[ \\t]*//'| tr -d "'" | tr " " ";" | sed -e 's/\$/;other-resistance-genes/') 
      printf "\${otherData}\\n" >> !{name}_!{method}.csv
    """
}

// tempDir may be neccessary