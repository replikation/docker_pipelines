process bakta_database {
        storeDir "${params.databases}/bakta_V4" 
        label 'ubuntu'    
      output:
        path("db.tar.gz")
      script:
        """
        wget --no-check-certificate https://zenodo.org/record/7025248/files/db.tar.gz
        """
    }