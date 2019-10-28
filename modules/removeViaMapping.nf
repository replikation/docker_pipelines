process removeViaMapping {
        publishDir "${params.output}/${name}", mode: 'copy', pattern: "${name}_*_non_mappers.fastq.gz"
        label 'seqtk'
    input:
        tuple val(name), file(readlist), file(shortreads)
    output:
        tuple val(name), file("${name}_forward_non_mappers.fastq.gz"), file("${name}_reverse_non_mappers.fastq.gz")
    script:
        """    
        seqtk subseq ${shortreads[0]} ${readlist} | gzip > ${name}_forward_non_mappers.fastq.gz
        seqtk subseq ${shortreads[1]} ${readlist} | gzip > ${name}_reverse_non_mappers.fastq.gz
        """
}