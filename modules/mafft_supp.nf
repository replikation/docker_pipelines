process mafft_supp {
    label 'mafft'
  input:
    tuple val(name), path(dir)
    path(supplement_file)
  output:
	  tuple val(name), path("protein_alignments.msa")
  script:
    if (!params.filenames)
    """
    cat ${dir}/* > all_proteins.aa
    cat ${supplement_file} >> all_proteins.aa
  	mafft --thread ${task.cpus} --auto all_proteins.aa > protein_alignments.msa
    """
    else if (params.filenames)
    """
    for file in ${dir}/* ; do
      filename=\$(basename \${file})
      cat \${file} | sed "1s/.*/>\${filename%.*}/"   >> all_proteins.aa
    done
    cat ${supplement_file} >> all_proteins.aa
    
  	mafft --thread ${task.cpus} --auto all_proteins.aa > protein_alignments.msa
    """
}

