process toytree {
    publishDir "${params.output}/${name}/tree/", mode: 'copy'
    label 'toytree'
  input:
    tuple val(name), path(tree)
  output:
	  tuple path("tree.svg"), path("tree.pdf"), path("tree_radial.svg"), path("tree_radial.pdf")
  script:
    """
    cp ${tree} tree_INPUT.newick
    vis_tree_features.py
    vis_tree_features_radial.py
    """
}



/*


*/
