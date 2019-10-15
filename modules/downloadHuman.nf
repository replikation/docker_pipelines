process downloadHuman {
    label 'ubuntu'
    storeDir 'stored-human-genome/'
  output:
    file("*.fna")
  shell:
    """
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
    gzip -d GCF_000001405.39_GRCh38.p13_genomic.fna.gz
    """
}