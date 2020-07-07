process centrifuge_download_db {
        if (params.cloudProcess) { publishDir "${params.database}/centrifuge", mode: 'copy', pattern: "gtdb_r89_54k_centrifuge.tar"}
        else { storeDir "${params.database}/centrifuge" } 
        label 'ubuntu'    
      output:
        file("gtdb_r89_54k_centrifuge.tar")
      script:
        """
        wget https://monash.figshare.com/ndownloader/files/16378439 -O gtdb_r89_54k_centrifuge.tar
        """
    }