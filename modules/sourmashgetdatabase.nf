process sourmash_download_db {
        if (workflow.profile == 'gcloud') { publishDir "${params.database}/sourmash", mode: 'copy', pattern: "GTDB-k31.lca.json.gz" }
        else { storeDir "${params.database}/sourmash" }  
        label 'ubuntu' 
      output:
        file("GTDB-k31.lca.json.gz")
      script:
        """
        wget https://osf.io/gs29b/download -O GTDB-k31.lca.json.gz
        """
    }

    