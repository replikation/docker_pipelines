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
