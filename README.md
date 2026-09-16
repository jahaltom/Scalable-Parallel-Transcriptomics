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
conda create -n snakemake -c conda-forge -c bioconda snakemake=7.32.4
conda activate snakemake
```

The workflow itself uses per-rule Conda environments.

## 2. Prepare references

You need:


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

## 3. Edit metadata.tsv

Required columns:

```text
sample    type    layout    read1    read2    sra    bam    Treatment
```

Examples are included.

Rules:

- `type=fastq`: provide `read1`; provide `read2` for PE.
- `type=sra`: provide the accession in `sra`.
- `type=bam`: provide the BAM path in `bam`.
- `layout` must be `SE` or `PE`.
- Leave unused fields blank.
- Sample names must be unique.

| sample   | type  | layout | read1                 | read2                 | sra       | bam          | Treatment    |
|----------|-------|--------|-----------------------|-----------------------|-----------|--------------|--------------|
| Patient1 | fastq | PE     | Patient1_R1.fastq.gz  | Patient1_R2.fastq.gz  |           |              |   Tumor      |
| Patient2 | sra   | PE     |                       |                       | SRR123456 |              |   Control    |
| Patient3 | bam   | PE     |                       |                       |           | Patient3.bam |              |
| Patient4 | fastq | SE     | Patient4.fastq.gz     |                       |           |              |              |


For BAM input, reads are extracted with `samtools fastq`. This is appropriate for ordinary genome-aligned BAMs because Salmon is then run from the reconstructed FASTQs rather than treating a genomic BAM as a transcriptome alignment.

## 4. Configure

Edit `config.yaml`:

```yaml
samples: "metadata.tsv"

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
snakemake -s Snakefile -n
```

## 6. Run

```bash
snakemake -j 24 -s Snakefile --use-conda --rerun-incomplete --latency-wait 60 --cluster "sbatch -t 05:00:00 -c {threads} -N 1"
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

Gene_Level_Raw_Counts.txt

DESeq2_results/
├── QC/
├── PCA/
├── DGE/
├── Volcano/
├── Heatmaps/
├── GSEA/
├── DESeq2_Normalized_Counts.txt
├── VST_Counts.txt
├── sessionInfo.txt
└── .complete
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

## Stand alone DGE.r

Needs Gene_Level_Raw_Counts.txt and metadata.tsv

Gene_Level_Raw_Counts.txt:


| Gene            | SRR..  |
|-----------------|--------|
| ENSG00000223972 | fastq  |

```
sbatch scripts/DGE.sh
```


