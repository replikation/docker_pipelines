process gtdbtk_download_db {
        storeDir "${params.databases}/gtdbtk" 
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