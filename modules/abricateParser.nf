process abricateParser {
    publishDir "${params.output}/${name}/ABR-Screening", mode: 'copy', pattern: "*.csv"
    label 'ubuntu'   
  input:
    set val(name), file(results) 
  output:
	  set val(name), file("*.csv") 
  shell:
    """
    printf "amount;type\\n" > !{name}_ncbi.csv
  	blaData=\$(grep -w "ncbi" !{results} | cut -f5 | grep "bla" | cut -f1 -d "-"| sort | uniq -c | sed -e 's/^[ \\t]*//' | tr " " ";") 
    printf "\${blaData}\\n" >> !{name}_ncbi.csv
    """
}

// tempDir may be neccessary