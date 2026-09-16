#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node 24
#SBATCH -t 4:00:00
#SBATCH --mail-user=haltomj@chop.edu
#SBATCH --mail-type=ALL

```
########################################################
# Download
########################################################
mkdir reference

# GRCh38.d1.vd1.fa.tar.gz:
wget -O reference/GRCh38.d1.vd1.fa.tar.gz https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834

# gencode.v36.annotation.gtf.gz:
wget -O reference/gencode.v36.annotation.gtf.gz https://api.gdc.cancer.gov/data/be002a2c-3b27-43f3-9e0f-fd47db92a6b5


#  gencode.v36.transcripts.fa.gz       
wget -O reference/gencode.v36.transcripts.fa.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.transcripts.fa.gz

########################################################
# DECOMPRESS
########################################################

# Extract .tar.gz genome archive
tar -xzf reference/GRCh38.d1.vd1.fa.tar.gz \
    -C reference/

# Decompress ordinary .gz files
gunzip reference/gencode.v36.annotation.gtf.gz
gunzip reference/gencode.v36.transcripts.fa.gz


########################################################
# generate tx2gene
########################################################

awk '$3=="transcript" {
    match($0, /gene_id "([^"]+)"/, g);
    match($0, /transcript_id "([^"]+)"/, t);
    if (g[1] != "" && t[1] != "")
        print t[1] "\t" g[1]
}' OFS='\t' reference/gencode.v36.annotation.gtf \
>  reference/tx2gene.tsv
sed -i '1i transcript_id\tgene_id' reference/tx2gene.tsv


rm reference/GRCh38.d1.vd1.fa.tar.gz
```
