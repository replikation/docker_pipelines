process spades {
    label 'spades'  
    publishDir "${params.output}/${name}/", mode: 'copy', pattern: "${name}.*fa*"

    input:
        tuple val(name), path(reads)
    output:
        tuple val(name), path(reads), path("${name}.fasta"), emit: spade_fasta_ch
        tuple val(name), path("${name}.gfa"), emit: spade_assembly_graph_ch
    script:
        """
        spades.py --careful --cov-cutoff auto -1 ${reads[0]} -2 ${reads[1]} -t ${task.cpus} -o output
        mv output/scaffolds.fasta ${name}.fasta
        mv output/assembly_graph_after_simplification.gfa ${name}.gfa
        """
}