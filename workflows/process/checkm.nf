process checkm {
        label 'checkm'
        publishDir "${params.output}/${name}/quality_of_bins/", mode: 'copy', pattern: "summary.txt"
        publishDir "${params.output}/${name}/quality_of_bins/", mode: 'copy', pattern: "taxonomy.txt"
        publishDir "${params.output}/${name}/quality_of_bins/", mode: 'copy', pattern: "bin_qa_plot.png"
    input:
        tuple val(name), path(bins)
    output:
        tuple val(name), path("summary.txt")
        tuple val(name), path("taxonomy.txt"), file("bin_qa_plot.png")
    script:
        """
        SUFFIXNAME=\$(ls ${bins} | head -1 | rev | cut -f1 -d "." | rev)

        mkdir temporary
        mkdir ${name}_bin
        mv ${bins}/*.fa* ${name}_bin/

        checkm lineage_wf --tmpdir temporary --pplacer_threads 4 -t ${task.cpus} --reduced_tree -x \$SUFFIXNAME ${name}_bin ${name}_checkm > summary.txt
        checkm bin_qa_plot --image_type png -x \$SUFFIXNAME ${name}_checkm ${name}_bin ${name}_checkm_plot
        checkm tree_qa ${name}_checkm > taxonomy.txt

        mv ${name}_checkm_plot/bin_qa_plot.png bin_qa_plot.png
        """
    stub:
    """
    touch summary.txt taxonomy.txt bin_qa_plot.png
    """
}
