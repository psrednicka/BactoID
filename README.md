**BactoID** is a Snakemake workflow for consensus-based identification of bacterial isolates from Oxford Nanopore Technologies (ONT) full-length 16S rRNA amplicon sequencing data.

The workflow processes demultiplexed ONT FASTQ files, generates high-quality consensus sequences, compares them against a local curated 16S rRNA reference database, and produces a structured Excel report containing taxonomic identification, quality metrics and interpretation flags.



## Workflow overview

BactoID performs the following steps:

1. Merge FASTQ files for each barcode.
2. Filter reads by length (default: 1000–2000 bp).
3. Cluster reads and generate consensus sequences using **NGSpeciesID**.
4. Polish consensus sequences using **Racon (3 iterations)**.
5. Search consensus sequences against a local curated **NCBI 16S rRNA database** using **BLAST+**.
6. Retain the best taxonomic matches.
7. Evaluate identification quality and taxonomic ambiguity.
8. Generate a structured Excel report.

## Requirements
BactoID is designed to run on Linux using:
- Snakemake
- Singularity / Apptainer
- Python 3
### The main bioinformatics tools are provided through containers.
Current workflow components include:
- NGSpeciesID
- Racon
- BLAST+
- Python / XlsxWriter

## Input data
BactoID expects demultiplexed ONT FASTQ files arranged in barcode directories:
```text
fastq_pass/
├── barcode01/
│   ├── file1.fastq.gz
│   └── file2.fastq.gz
├── barcode02/
│   ├── file1.fastq.gz
│   └── file2.fastq.gz
├── barcode03/
│   └── ...
└── ...
```

Both .fastq and .fastq.gz files are supported.
Barcode directories are detected automatically.

## Optional sample manifest
An optional tab-separated manifest can be used to replace technical barcode names with user-defined sample names in the final report.

Example manifest.tsv:
```text
Barcode	SampleName
barcode01	KKP_4195
barcode02	KKP_4196
barcode03	Control_1
```

The workflow internally continues to use barcode names, while the final Excel report displays the corresponding sample names.
If no manifest is provided, barcode names are retained.

## Running BactoID
Example:

```bash 
cd /path/to/BactoID

snakemake \
  --cores 8 \
  --use-singularity \
  --config \
  input=/path/to/fastq_pass \
  output=/path/to/BactoID_results

  ```

  With a sample manifest:

 ```bash
 snakemake \
  --cores 8 \
  --use-singularity \
  --config \
  input=/path/to/fastq_pass \
  output=/path/to/BactoID_results \
  manifest=/path/to/manifest.tsv
  ```

   On systems requiring explicit directory binding:

```bash
snakemake \
  --cores 8 \
  --use-singularity \
  --singularity-args "--bind /path/to/data:/path/to/data" \
  --config \
  input=/path/to/fastq_pass \
  output=/path/to/BactoID_results \
  manifest=/path/to/manifest.tsv
  ```

### Read filtering
BactoID is intended for bacterial full-length 16S rRNA amplicons.
Default accepted read length: **1000–2000 bp**. Reads outside this range are removed before consensus generation.
Filtering statistics are retained for each sample.

### Consensus generation
Filtered ONT reads are processed using NGSpeciesID.
```text
Default parameters:
Expected amplicon length: 1500 bp
Allowed deviation:        ±500 bp
Consensus generation:     enabled
Polishing:                Racon
Racon iterations:         3
```
Each supported sequence cluster is retained as an independent consensus and analysed separately. This allows BactoID to identify potential mixed cultures instead of reporting only the largest cluster.

### Taxonomic identification
Consensus sequences are compared against a local curated NCBI 16S rRNA reference database using blastn.
BactoID uses a minimum sequence identity threshold of: **97%**.
Consensus sequences below this threshold are reported as: **NO_IDENTIFICATION**.
The closest database match is nevertheless retained for technical review.
A reference coverage of at least: **95%** is recommended.
Lower values are flagged as: **LOW_REFERENCE_COVERAGE**
```text
Read support
Supporting reads	Status
≥100	adequate support
20–99	LOW_READ_COUNT
1–19	VERY_LOW_READ_COUNT
```

### Multiple clusters
BactoID does not automatically assume that every sample contains a single organism.
Significant clusters are included in the summary when they contain:
- at least 100 supporting reads, and
- at least 10% of the read count of the largest identified cluster.
Several species within the same genus
If multiple supported assignments belong to the same genus, BactoID reports:
SPECIES_UNRESOLVED_BY_16S

Example:
Lactiplantibacillus plantarum / L. pentosus

This indicates insufficient discriminatory power of the 16S rRNA gene and does not by itself prove that the culture is mixed.
Different genera
If substantial independent clusters are assigned to different genera, BactoID reports:
POSSIBLE_MIXED_CULTURE

```text
Output
Main output:
BactoID_results/
├── report/
│   ├── BactoID_results.tsv
│   └── BactoID_results.xlsx
├── intermediate/
│   ├── merged/
│   ├── filtered/
│   ├── ngspeciesid/
│   ├── selected/
│   ├── blast/
│   └── classification/
└── logs/
```
## Important limitations
BactoID performs taxonomic identification based on the bacterial 16S rRNA gene.
Closely related bacterial species may have identical or nearly identical 16S rRNA sequences and therefore cannot always be reliably distinguished at species level.
A species-level result should therefore be interpreted together with:
- sequence identity,
- reference coverage,
- taxonomic ambiguity,
- read support,
- biological context.
For isolates requiring definitive species or strain identification, whole-genome sequencing is recommended.

## Validation
During workflow development, different consensus-polishing strategies were evaluated using experimentally obtained full-length ONT 16S datasets.
For the tested datasets, three iterations of Racon polishing provided more consistent agreement with curated reference sequences than Medaka polishing and was therefore selected as the default BactoID consensus-polishing strategy.

