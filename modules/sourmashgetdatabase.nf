process sourmash_download_db {
        if (workflow.profile == 'gcloud') { publishDir 'gs://databases-nextflow/databases/sourmash', mode: 'copy', pattern: "GTDB-k31.lca.json.gz" }
        else { storeDir 'nextflow-autodownload-databases/sourmash' }  
        label 'ubuntu' 
      output:
        file("GTDB-k31.lca.json.gz")
      script:
        """
        wget https://osf.io/gs29b/download -O GTDB-k31.lca.json.gz
        """
    }

    