process gtdbtk_download_db {
        if (workflow.profile == 'gcloud') { publishDir 'gs://databases-nextflow/databases/gtdbtk', mode: 'copy', pattern: "gtdbtk_r89_data.tar.gz"}
        else { storeDir 'nextflow-autodownload-databases/gtdbtk' }  
        label 'ubuntu'    
      output:
        file("gtdbtk_r89_data.tar.gz")
      script:
        """
        wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz
        """
    }