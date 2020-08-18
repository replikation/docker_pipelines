process cluster {
    label 'cd_hit'   
    publishDir "${params.output}/${name}", mode: 'copy', pattern: "${name}_representatives.fasta"
  input:
    tuple val(name), file(fasta)
  output:
    tuple val(name), file("${name}_representatives.fasta")
  script:
    """
    psi-cd-hit.pl -i ${fasta} -o output -aL 0.6 -prog blastn -para ${task.cpus} -c 0.6 -g 1 \
    -s "-evalue 10E-100 -max_target_seqs 100000 -qcov_hsp_perc 10 -max_hsps 10"
    grep ">" output > representatives.txt
    mv output ${name}_representatives.fasta
    """
}