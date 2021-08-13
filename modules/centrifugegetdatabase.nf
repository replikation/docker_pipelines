process centrifuge_download_db {
        storeDir "${params.databases}/centrifuge"
        label 'ubuntu'    
      output:
        file("gtdb_r89_54k_centrifuge.tar")
      script:
        """
        wget https://monash.figshare.com/ndownloader/files/16378439 -O gtdb_r89_54k_centrifuge.tar
        """
    }