process downloadHuman {
    label 'ubuntu'
    storeDir '/tmp/stored-human-genome'
  output:
    file("GCF_000001405.39_GRCh38.p13_genomic.fna")
  shell:
    """
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
    zcat GCF_000001405.39_GRCh38.p13_genomic.fna.gz > GCF_000001405.39_GRCh38.p13_genomic.fna
    rm GCF_000001405.39_GRCh38.p13_genomic.fna.gz
    """
}