process sourmash_download_db {
        storeDir 'nextflow-autodownload-databases/sourmash'
      output:
        file("genbank-k31.lca.json")
      script:
        """
        wget https://osf.io/4f8n3/download -O genbank-k31.lca.json.gz
        gunzip genbank-k31.lca.json.gz
        """
    }