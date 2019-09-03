// check storeDit

process gtdbtk_download_db {
        storeDir 'nextflow-autodownload-databases/gtdbtk'
      output:
        file("release89")
      script:
        """
        wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz
        tar xvzf gtdbtk_r89_data.tar.gz
        rm gtdbtk_r89_data.tar.gz
        """
    }