# Building a database bacteria viral human

## 1. Start docker

````bash
docker run -it -v $PWD:/input replikation/centrifuge_small /bin/bash
````

## 2. Within the docker

* ``centrifuge-download`` takes a while
* to much cores in ``centrifuge-build`` can do a "memory leak" and signal 9 kill
  * 4 cores are tested and worked
  * building took 20 h

````bash
cd /input && \
centrifuge-download -o taxonomy taxonomy && \
centrifuge-download -o library -m -d "archaea,bacteria,viral" refseq > seqid2taxid.map && \
cat library/*/*.fna > input-sequences.fna && \
centrifuge-build -p 4 --conversion-table seqid2taxid.map \
                 --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
                 input-sequences.fna abv_4t --noauto
````

# Database sizes for bacterial genomes

* all bacteria, doesnt matter wich assembly level:

## genBank complete versus incomplete

* genBank, including incomplete:
  * ``centrifuge-download -o library -a Any -m -d "bacteria" genbank > seqid2taxid.map``
  * 216860 bacteria genomes at assembly level Any (genBank)
* genBank, including complete:
  * ``centrifuge-download -o library -a Any -m -d "bacteria" genbank > seqid2taxid.map``
  *  13779 bacteria genomes at assembly level Complete Genome

## refseq complete versus incomplete

* refseq, including complete:
  * ``centrifuge-download -o library -m -d "bacteria" refseq > seqid2taxid.map``
  * 13193 bacteria genomes at assembly level Complete Genome (refseq)
* refseq, including incomplete:
  * ``centrifuge-download -o library -a Any -m -d "bacteria" refseq > seqid2taxid.map``
  * 150900 bacteria genomes at assembly level Any

## good but big database:

````bash
centrifuge-download -o library -a Any -m -d "archaea,bacteria" refseq > seqid2taxid.map
````

HELP for centrifuge-download:


````bash
centrifuge-download [<options>] <database>

ARGUMENT
 <database>        One of refseq, genbank, contaminants or taxonomy:
                     - use refseq or genbank for genomic sequences,
                     - contaminants gets contaminant sequences from UniVec and EmVec,
                     - taxonomy for taxonomy mappings.

COMMON OPTIONS
 -o <directory>         Folder to which the files are downloaded. Default: '.'.
 -P <# of threads>      Number of processes when downloading (uses xargs). Default: '1'

WHEN USING database refseq OR genbank:
 -d <domain>            What domain to download. One or more of bacteria, viral, archaea, fungi, protozoa, invertebrate, plant, vertebrate_mammalian, vertebrate_other (comma separated).
 -a <assembly level>    Only download genomes with the specified assembly level. Default: 'Complete Genome'. Use 'Any' for any assembly level.
 -c <refseq category>   Only download genomes in the specified refseq category. Default: any.
 -t <taxids>            Only download the specified taxonomy IDs, comma separated. Default: any.
 -g <program>           Download using program. Options: rsync, curl, wget. Default wget (auto-detected).
 -r                     Download RNA sequences, too.
 -u                     Filter unplaced sequences.
 -m                     Mask low-complexity regions using dustmasker. Default: off.
 -l                     Modify header to include taxonomy ID. Default: off.
 -g                     Download GI map.
 -v                     Verbose mode
 ````






## 3. Need more RAM?
### create a big temporay swap until next reboot
````bash
# check the name of the actual swap file!!
# on ubuntu 18 its swap.img
  sudo fallocate -l 50G /swapfile
  sudo chmod 600 /swapfile
  sudo mkswap /swapfile
  sudo swapon /swapfile
# remove /swap after reboot
  sudo rm -f /swapfile
````

* more here: https://linuxize.com/post/create-a-linux-swap-file/

# sourmash

* include DB into docker (its only a few MB)
docker run -v $PWD:/input -it replikation/sourmash /bin/bash

sourmash compute --scaled 1000 -k 31  polished.assembly.fa -o /input/signature.sig

sourmash gather -k 31 signature.sig genbank-k31.lca.json.gz -o results

**for taxonomic classification it seems i have to extract each fasta into a file an compute them???**
