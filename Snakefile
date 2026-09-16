import os
import pandas as pd
from pathlib import Path

configfile: "config.yaml"

SAMPLES_FILE = config["samples"]
samples = pd.read_csv(SAMPLES_FILE, sep="\t", dtype=str).fillna("")

required_cols = {"sample", "type", "layout", "read1", "read2", "sra", "bam"}
missing = required_cols - set(samples.columns)
if missing:
    raise ValueError(f"metadata.tsv is missing columns: {sorted(missing)}")

if samples["sample"].duplicated().any():
    dup = samples.loc[samples["sample"].duplicated(), "sample"].tolist()
    raise ValueError(f"Duplicate sample names: {dup}")

samples = samples.set_index("sample", drop=False)
SAMPLE_NAMES = samples.index.tolist()

for s, row in samples.iterrows():
    if row["type"] not in {"fastq", "sra", "bam"}:
        raise ValueError(f"{s}: type must be fastq, sra, or bam")
    if row["layout"] not in {"SE", "PE"}:
        raise ValueError(f"{s}: layout must be SE or PE")
    if row["type"] == "fastq":
        if not row["read1"]:
            raise ValueError(f"{s}: FASTQ sample requires read1")
        if row["layout"] == "PE" and not row["read2"]:
            raise ValueError(f"{s}: paired FASTQ sample requires read2")
    elif row["type"] == "sra" and not row["sra"]:
        raise ValueError(f"{s}: SRA sample requires sra accession")
    elif row["type"] == "bam" and not row["bam"]:
        raise ValueError(f"{s}: BAM sample requires bam path")

BUILD_INDEX = bool(config.get("build_index", True))
SALMON_INDEX = config["salmon_index"]
GENOME = config["reference"]["genome"]
TRANSCRIPTS = config["reference"]["transcripts"]
TX2GENE = config["reference"]["tx2gene"]

def fq1(sample):
    return f"work/fastq/{sample}_R1.fastq.gz"

def fq2(sample):
    return f"work/fastq/{sample}_R2.fastq.gz"

def normalized_fastqs(wc):
    row = samples.loc[wc.sample]
    if row["layout"] == "PE":
        return [fq1(wc.sample), fq2(wc.sample)]
    return [fq1(wc.sample)]

def salmon_reads(wc):
    row = samples.loc[wc.sample]
    if row["layout"] == "PE":
        return [fq1(wc.sample), fq2(wc.sample)]
    return [fq1(wc.sample)]

def quant_cmd(wc, input):
    row = samples.loc[wc.sample]
    if row["layout"] == "PE":
        return f"-1 {input[0]} -2 {input[1]}"
    return f"-r {input[0]}"

rule all:
    input:
        expand("results/salmon/{sample}/quant.sf", sample=SAMPLE_NAMES),
        "results/matrices/gene_counts.tsv",
        "results/matrices/gene_tpm.tsv",
        "results/matrices/gene_length.tsv",
        "results/multiqc/multiqc_report.html",
        "DESeq2_results/.complete"

# -------------------------------------------------------------------------
# Reference / Salmon index
# -------------------------------------------------------------------------

if BUILD_INDEX:
    rule make_decoys:
        input:
            genome=GENOME
        output:
            "resources/decoys.txt"
        conda:
            "envs/salmon.yaml"
        shell:
            r"""
            mkdir -p resources
            grep '^>' {input.genome:q} | cut -d ' ' -f 1 | sed 's/^>//' > {output:q}
            """

    rule make_gentrome:
        input:
            transcripts=TRANSCRIPTS,
            genome=GENOME
        output:
            "resources/gentrome.fa"
        shell:
            r"""
            mkdir -p resources
            cat {input.transcripts:q} {input.genome:q} > {output:q}
            """

    rule salmon_index:
        input:
            gentrome="resources/gentrome.fa",
            decoys="resources/decoys.txt"
        output:
            directory(SALMON_INDEX)
        threads:
            config["threads"].get("index", 12)
        conda:
            "envs/salmon.yaml"
        shell:
            r"""
            salmon index \
              -t {input.gentrome:q} \
              -d {input.decoys:q} \
              -i {output:q} \
              -p {threads}
            """
else:
    # Validate existence at parse time for a pre-built index.
    if not os.path.isdir(SALMON_INDEX):
        raise ValueError(
            f"build_index is false but Salmon index directory does not exist: {SALMON_INDEX}"
        )

# -------------------------------------------------------------------------
# Input normalization: FASTQ
# -------------------------------------------------------------------------

FASTQ_SE = [s for s in SAMPLE_NAMES if samples.loc[s, "type"] == "fastq" and samples.loc[s, "layout"] == "SE"]
FASTQ_PE = [s for s in SAMPLE_NAMES if samples.loc[s, "type"] == "fastq" and samples.loc[s, "layout"] == "PE"]
SRA_SE   = [s for s in SAMPLE_NAMES if samples.loc[s, "type"] == "sra"   and samples.loc[s, "layout"] == "SE"]
SRA_PE   = [s for s in SAMPLE_NAMES if samples.loc[s, "type"] == "sra"   and samples.loc[s, "layout"] == "PE"]
BAM_SE   = [s for s in SAMPLE_NAMES if samples.loc[s, "type"] == "bam"   and samples.loc[s, "layout"] == "SE"]
BAM_PE   = [s for s in SAMPLE_NAMES if samples.loc[s, "type"] == "bam"   and samples.loc[s, "layout"] == "PE"]

rule stage_fastq_se:
    input:
        lambda wc: samples.loc[wc.sample, "read1"]
    output:
        "work/fastq/{sample}_R1.fastq.gz"
    wildcard_constraints:
        sample="|".join(map(str, FASTQ_SE)) if FASTQ_SE else r"a^"
    conda:
        "envs/core.yaml"
    shell:
        r"""
        mkdir -p work/fastq
        scripts/stage_fastq.sh {input:q} {output:q}
        """

rule stage_fastq_pe:
    input:
        r1=lambda wc: samples.loc[wc.sample, "read1"],
        r2=lambda wc: samples.loc[wc.sample, "read2"]
    output:
        r1="work/fastq/{sample}_R1.fastq.gz",
        r2="work/fastq/{sample}_R2.fastq.gz"
    wildcard_constraints:
        sample="|".join(map(str, FASTQ_PE)) if FASTQ_PE else r"a^"
    conda:
        "envs/core.yaml"
    shell:
        r"""
        mkdir -p work/fastq
        scripts/stage_fastq.sh {input.r1:q} {output.r1:q}
        scripts/stage_fastq.sh {input.r2:q} {output.r2:q}
        """

rule sra_to_fastq_se:
    output:
        "work/fastq/{sample}_R1.fastq.gz"
    params:
        acc=lambda wc: samples.loc[wc.sample, "sra"]
    threads:
        config["threads"].get("sra_to_fastq", 6)
    wildcard_constraints:
        sample="|".join(map(str, SRA_SE)) if SRA_SE else r"a^"
    conda:
        "envs/sra.yaml"
    shell:
        r"""
        mkdir -p work/fastq work/sra_tmp/{wildcards.sample}
        fasterq-dump {params.acc:q} \
          --threads {threads} \
          --outdir work/sra_tmp/{wildcards.sample}

        f="work/sra_tmp/{wildcards.sample}/{params.acc}.fastq"
        f1="work/sra_tmp/{wildcards.sample}/{params.acc}_1.fastq"

        if [ -s "$f" ]; then
            gzip -c "$f" > {output:q}
        elif [ -s "$f1" ]; then
            gzip -c "$f1" > {output:q}
        else
            echo "Could not find single-end FASTQ generated by fasterq-dump" >&2
            exit 1
        fi
        rm -rf work/sra_tmp/{wildcards.sample}
        """

rule sra_to_fastq_pe:
    output:
        r1="work/fastq/{sample}_R1.fastq.gz",
        r2="work/fastq/{sample}_R2.fastq.gz"
    params:
        acc=lambda wc: samples.loc[wc.sample, "sra"]
    threads:
        config["threads"].get("sra_to_fastq", 6)
    wildcard_constraints:
        sample="|".join(map(str, SRA_PE)) if SRA_PE else r"a^"
    conda:
        "envs/sra.yaml"
    shell:
        r"""
        mkdir -p work/fastq work/sra_tmp/{wildcards.sample}
        fasterq-dump {params.acc:q} \
          --split-files \
          --threads {threads} \
          --outdir work/sra_tmp/{wildcards.sample}

        r1="work/sra_tmp/{wildcards.sample}/{params.acc}_1.fastq"
        r2="work/sra_tmp/{wildcards.sample}/{params.acc}_2.fastq"

        test -s "$r1"
        test -s "$r2"

        gzip -c "$r1" > {output.r1:q}
        gzip -c "$r2" > {output.r2:q}
        rm -rf work/sra_tmp/{wildcards.sample}
        """

rule bam_to_fastq_se:
    input:
        lambda wc: samples.loc[wc.sample, "bam"]
    output:
        "work/fastq/{sample}_R1.fastq.gz"
    threads:
        config["threads"].get("bam_to_fastq", 4)
    wildcard_constraints:
        sample="|".join(map(str, BAM_SE)) if BAM_SE else r"a^"
    conda:
        "envs/samtools.yaml"
    shell:
        r"""
        mkdir -p work/fastq
        samtools collate -u -O {input:q} \
          | samtools fastq -@ {threads} - \
          | gzip -c > {output:q}
        """

rule bam_to_fastq_pe:
    input:
        lambda wc: samples.loc[wc.sample, "bam"]
    output:
        r1="work/fastq/{sample}_R1.fastq.gz",
        r2="work/fastq/{sample}_R2.fastq.gz"
    threads:
        config["threads"].get("bam_to_fastq", 4)
    wildcard_constraints:
        sample="|".join(map(str, BAM_PE)) if BAM_PE else r"a^"
    conda:
        "envs/samtools.yaml"
    shell:
        r"""
        mkdir -p work/fastq
        samtools collate -u -O {input:q} \
          | samtools fastq -@ {threads} \
              -1 >(gzip -c > {output.r1:q}) \
              -2 >(gzip -c > {output.r2:q}) \
              -0 /dev/null \
              -s /dev/null \
              -n -
        """

# -------------------------------------------------------------------------
# QC
# -------------------------------------------------------------------------

rule fastqc:
    input:
        normalized_fastqs
    output:
        html="results/fastqc/{sample}.done"
    threads:
        config["threads"].get("fastqc", 2)
    conda:
        "envs/qc.yaml"
    shell:
        r"""
        mkdir -p results/fastqc
        fastqc --threads {threads} --outdir results/fastqc {input:q}
        touch {output.html:q}
        """

rule multiqc:
    input:
        expand("results/fastqc/{sample}.done", sample=SAMPLE_NAMES),
        expand("results/salmon/{sample}/quant.sf", sample=SAMPLE_NAMES)
    output:
        "results/multiqc/multiqc_report.html"
    conda:
        "envs/qc.yaml"
    shell:
        r"""
        mkdir -p results/multiqc
        multiqc results/fastqc results/salmon \
          --outdir results/multiqc \
          --filename multiqc_report.html \
          --force
        """

# -------------------------------------------------------------------------
# Salmon
# -------------------------------------------------------------------------

rule salmon_quant:
    input:
        reads=salmon_reads,
        index=lambda wc: SALMON_INDEX
    output:
        quant="results/salmon/{sample}/quant.sf"
    params:
        reads=quant_cmd,
        libtype=config["salmon"].get("libtype", "A"),
        validate="--validateMappings" if config["salmon"].get("validateMappings", True) else "",
        seqbias="--seqBias" if config["salmon"].get("seqBias", True) else "",
        gcbias="--gcBias" if config["salmon"].get("gcBias", True) else "",
        extra=config["salmon"].get("extra", "")
    threads:
        config["threads"].get("salmon", 8)
    conda:
        "envs/salmon.yaml"
    shell:
        r"""
        mkdir -p results/salmon/{wildcards.sample}
        salmon quant \
          -i {input.index:q} \
          -l {params.libtype:q} \
          {params.reads} \
          -p {threads} \
          {params.validate} \
          {params.seqbias} \
          {params.gcbias} \
          {params.extra} \
          -o results/salmon/{wildcards.sample}
        """

# -------------------------------------------------------------------------
# tximport matrices
# -------------------------------------------------------------------------

rule tximport:
    input:
        quants=expand("results/salmon/{sample}/quant.sf", sample=SAMPLE_NAMES),
        tx2gene=TX2GENE
    output:
        counts="results/matrices/gene_counts.tsv",
        tpm="results/matrices/gene_tpm.tsv",
        length="results/matrices/gene_length.tsv"
    params:
        samples=",".join(SAMPLE_NAMES)
    conda:
        "envs/r_tximport.yaml"
    shell:
        r"""
        mkdir -p results/matrices
        Rscript scripts/tximport.R \
          --samples {params.samples:q} \
          --tx2gene {input.tx2gene:q} \
          --salmon-dir results/salmon \
          --outdir results/matrices
        """

# -------------------------------------------------------------------------
# PREPARE COUNTS FOR DESEQ2
# -------------------------------------------------------------------------

rule prepare_deseq_counts:
    input:
        "results/matrices/gene_counts.tsv"
    output:
        "Gene_Level_Raw_Counts.txt"
    shell:
        r"""
        awk 'BEGIN{{FS=OFS="\t"}} NR==1{{$1="Gene"}} NR>1{{sub(/\.[0-9]+$/, "", $1)}} 1' \
            {input:q} > {output:q}
        """
# -------------------------------------------------------------------------
# DESEQ2
# -------------------------------------------------------------------------

rule deSeq:
    input:
        counts="results/matrices/gene_counts.tsv"
    output:
        done="DESeq2_results/.complete"
    shell:
        r"""
        module load R
        Rscript DGE.r

        touch {output.done}
        """

