process krona {
        label 'krona'
        publishDir "${params.output}/${params.readqcdir}/${name}/", mode: 'copy'
    input:
        tuple val(name), path(kraken2), path(kreport)
  	output:
    	tuple val(name), file("${name}_krona.html")
  	script:
    """
    cat ${kreport} | cut -f 3,5 > file.krona
    ktImportTaxonomy file.krona -m 1
    mv *.html ${name}_krona.html
    """
    stub:
    """
    touch ${name}_krona.html
    """
}

process krona_from_bracken {
        label 'krona'
        publishDir "${params.output}/${name}/Read_classification/", mode: 'copy'
    input:
        tuple val(name), path(alphadiversity), path(kronatextinput)
  	output:
    	tuple val(name), file("${name}.krona.html")
  	script:
    """
    ktImportText ${kronatextinput} -o ${name}.krona.html
    """
    stub:
    """
    touch ${name}_krona.html
    """
}


/*

python KrakenTools/kreport2krona.py -r breports/SRR14143424.breport -o b_krona_txt/SRR14143424.b.krona.txt --no-intermediate-ranks 
KronaScripts/ktImportText b_krona_txt/SRR14143424.b.krona.txt \ -o krona_html/SRR14143424.krona.html


*/