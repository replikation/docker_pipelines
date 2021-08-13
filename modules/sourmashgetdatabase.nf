process sourmash_download_db {
        storeDir "${params.databases}/sourmash"
        label 'sourmash' 
      output:
        path("gtdb-release202-k31.lca.json.gz")
      script:
        """
        wget https://osf.io/9xdg2/download -O gtdb-release202-k31.lca.json.gz
        """
    }
    