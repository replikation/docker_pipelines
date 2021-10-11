process bakta_database {
        storeDir "${params.databases}/bakta" 
        label 'ubuntu'    
      output:
        path("db.tar.gz")
      script:
        """
        wget --no-check-certificate https://zenodo.org/record/5215743/files/db.tar.gz
        """
    }