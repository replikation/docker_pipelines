process fasttree {
    label "fasttree"
    input:
        tuple val(name), path(alignment)
    output:
        tuple val(name), path("tree.nwk")
    script:
        """
        export OMP_NUM_THREADS=${task.cpus}
        FastTree -gtr  ${alignment} > tree.nwk
        """
}