process sourmash_download_db {
        publishDir "db_auto-build"
      output:
        file("genbank-k31.lca.json.gz")
      script:
        """
        wget https://osf.io/4f8n3/download -O genbank-k31.lca.json.gz
        """
    }