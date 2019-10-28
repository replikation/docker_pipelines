process centrifuge_download_db {
        if (workflow.profile == 'gcloud') { publishDir 'gs://databases-nextflow/databases/centrifuge', mode: 'copy', pattern: "gtdb_r89_54k_centrifuge.tar"}
        else { storeDir 'nextflow-autodownload-databases/centrifuge' } 
        label 'ubuntu'    
      output:
        file("gtdb_r89_54k_centrifuge.tar")
      script:
        """
        wget https://monash.figshare.com/ndownloader/files/16378439 -O gtdb_r89_54k_centrifuge.tar
        """
    }