process download_database_kraken2 {
        label "ubuntu"
        storeDir "${params.databases}/kraken2_k2standard_20220926"
        errorStrategy 'retry'
        maxRetries 1
    output:
        path("kraken.tar.gz")
    script:
    if (task.attempt.toInteger() == 1)
        """
        echo ${task.attempt}
        wget --no-check-certificate https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20220926.tar.gz -O kraken.tar.gz
        """
    else if (task.attempt.toInteger() > 1)
        """
        echo ${task.attempt}
        wget --no-check-certificate https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20220926.tar.gz -O kraken.tar.gz
        """
    stub:
        """
        touch kraken.tar.gz
        """
}


/* 

DATABASES

https://benlangmead.github.io/aws-indexes/k2

*/