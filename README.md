# Salmon Snakemake RNA-seq pipeline

A reusable workflow that accepts **SRA accessions, FASTQ files, or BAM files** and converges them to a common Salmon quantification workflow.

## Workflow

SRA -> fasterq-dump -> FASTQ
BAM -> samtools fastq -> FASTQ
FASTQ -----------------> FASTQ
                         |
                         +-> FastQC
                         +-> Salmon quantification -> quant.sf
                                               |
                                               +-> tximport gene counts / TPM

The workflow supports single-end (SE) and paired-end (PE) samples.

## 1. Install Snakemake

A minimal option is:

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
```

The workflow itself uses per-rule Conda environments.

## 2. Prepare references

You need:

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



rm reference/GRCh38.d1.vd1.fa.tar.gz
```




- genome FASTA
- transcriptome FASTA
- transcript-to-gene table (`tx2gene.tsv`)

`tx2gene.tsv` must contain two tab-separated columns:

```text
transcript_id    gene_id
ENST000...       ENSG000...
```

The workflow builds a **decoy-aware Salmon index** by concatenating the transcriptome and genome and using genome contig names as decoys.

If you already have a Salmon index, set `salmon_index` in `config.yaml` and set `build_index: false`.

## 3. Edit samples.tsv

Required columns:

```text
sample    type    layout    read1    read2    sra    bam
```

Examples are included.

Rules:

- `type=fastq`: provide `read1`; provide `read2` for PE.
- `type=sra`: provide the accession in `sra`.
- `type=bam`: provide the BAM path in `bam`.
- `layout` must be `SE` or `PE`.
- Leave unused fields blank.
- Sample names must be unique.

```
sample      type    layout  read1                 read2                 sra          bam
Patient1    fastq   PE      Patient1_R1.fastq.gz  Patient1_R2.fastq.gz
Patient2    sra     PE                                                    SRR123456
Patient3    bam     PE                                                                 Patient3.bam
Patient4    fastq   SE      Patient4.fastq.gz
```
For BAM input, reads are extracted with `samtools fastq`. This is appropriate for ordinary genome-aligned BAMs because Salmon is then run from the reconstructed FASTQs rather than treating a genomic BAM as a transcriptome alignment.

## 4. Configure

Edit `config.yaml`:

```yaml
samples: "samples.tsv"

reference:
  genome: "/path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
  transcripts: "/path/to/Homo_sapiens.GRCh38.cdna.all.fa"
  tx2gene: "/path/to/tx2gene.tsv"

build_index: true
salmon_index: "resources/salmon_index"

salmon:
  libtype: "A"
  validateMappings: true
  seqBias: true
  gcBias: true
```

`libtype: A` lets Salmon infer library type.

## 5. Dry run

```bash
snakemake -n -p
```

## 6. Run

```bash
snakemake --use-conda --cores 16
```

On a cluster, use your site's Snakemake executor/profile.

## Outputs

```text
results/
├── fastqc/
├── multiqc/
│   └── multiqc_report.html
├── salmon/
│   └── SAMPLE/
│       └── quant.sf
└── matrices/
    ├── gene_counts.tsv
    ├── gene_tpm.tsv
    └── gene_length.tsv
```

Intermediate normalized FASTQs are placed under `work/fastq/`.

## Important BAM caveat

Converting a BAM back to FASTQ can lose some information from the original raw FASTQs and depends on the BAM retaining the reads you need. If the BAM was heavily filtered, only retained alignments/reads can be recovered. Raw FASTQ is preferable when available.

## Recommended first test

Before launching hundreds/thousands of samples:

```bash
snakemake --use-conda --cores 8 results/salmon/YOUR_SAMPLE/quant.sf
```

Then inspect `results/salmon/YOUR_SAMPLE/logs/salmon_quant.log`.
