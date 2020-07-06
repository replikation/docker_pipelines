process gtdbtk_download_db {
        if (workflow.profile == 'gcloud') { publishDir "${params.database}/gtdbtk", mode: 'copy', pattern: "gtdb.lca.json"}
        else { storeDir "${params.database}/gtdbtk" }  
        label 'ubuntu'    
      output:
        file("gtdb.lca.json")
      script:
        """
        wget https://ndownloader.figshare.com/files/18809423?private_link=ed98a281ef089c033352 -O gtdb.lca.json
        """
    }

// from here:
// https://www.biorxiv.org/content/biorxiv/early/2019/07/23/712166.full.pdf