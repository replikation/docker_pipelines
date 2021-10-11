process gtdbtk_download_db {
        storeDir "${params.databases}/gtdbtk" 
        label 'ubuntu'    
      output:
        path("gtdbtk_data.tar.gz")
      script:
        """
        wget --no-check-certificate https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz
        """
    }

// from here:
// https://www.biorxiv.org/content/biorxiv/early/2019/07/23/712166.full.pdf