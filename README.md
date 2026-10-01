# BactoID

<p align="center">
  <img src="resources/logo/BactoID.png" width="600">
</p>

**BactoID** is a Snakemake workflow for consensus-based identification of bacterial isolates from Oxford Nanopore Technologies (ONT) full-length 16S rRNA amplicon sequencing data.

The workflow processes demultiplexed ONT FASTQ files, generates high-quality consensus sequences, compares them against a local curated 16S rRNA reference database, and produces a structured Excel report containing taxonomic identification, quality metrics and interpretation flags. BactoID is distributed as a single Singularity/Apptainer container containing all required software and reference resources.



## Workflow overview

BactoID performs the following steps:

1. automatic detection of barcode directories,
2. merging of FASTQ files for each barcode,
3. read-length filtering,
4. clustering and consensus generation with **NGSpeciesID**,
5. consensus polishing using **Racon ×3**,
6. taxonomic identification using **BLAST+**,
7. identification quality assessment,
8. detection of possible mixed cultures,
9. generation of a structured Excel report,
10. optional replacement of barcode names with user-defined sample names using a manifest file.

## The only external requirement is:

- Linux
- Singularity or Apptainer

All other software and reference resources are included in the BactoID container.

The container includes:

- Snakemake
- NGSpeciesID
- Racon
- BLAST+
- Python
- XlsxWriter
- Pillow
- curated NCBI 16S rRNA BLAST database
- BactoID workflow scripts
- report-generation scripts

## Download

BactoID v0.5.0 is distributed as a self-contained Singularity/Apptainer image.

**Version 0.5.0:**  
https://doi.org/10.5281/zenodo.23076229

**All versions / latest release:**  
https://doi.org/10.5281/zenodo.23076228

Container file:

```text
BactoID_0.5.0.sif
```


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

singularity exec \
  --bind /path/to/data:/data \
  BactoID_0.5.0.sif \
  bactoid \
  --input /data/fastq_pass \
  --output /data/BactoID_results \
  --cores 8
  ```

  With a sample manifest:

 ```bash
singularity exec \
  --bind /path/to/data:/data \
  BactoID_0.5.0.sif \
  bactoid \
  --input /data/fastq_pass \
  --output /data/BactoID_results \
  --manifest /data/fastq_pass/manifest.tsv \
  --cores 8
```

For Apptainer, replace:

```text
singularity exec
```

with:

```text
apptainer exec
```

## Command-line options

```text
--input       Path to the directory containing barcodeXX folders
--output      Path where BactoID results will be created
--manifest    Optional TSV manifest with Barcode and SampleName columns
--cores       Number of CPU cores to use (default: 8)
```

Display help:

```bash
singularity exec BactoID_0.5.0.sif bactoid --help
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
### Reference database

Taxonomic identification is performed against a local database derived from the NCBI RefSeq Targeted Loci 16S rRNA collection for bacterial and archaeal type material.

NCBI RefSeq Targeted Loci:
https://www.ncbi.nlm.nih.gov/refseq/targetedloci/

NCBI 16S ribosomal RNA database (Bacteria and Archaea type strains), snapshot downloaded on 23 July 2026.

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


## Output

```text
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

## Example report
A complete example of the generated BactoID Excel report is available below:

[📊 View example report](examples/BactoID_results.xlsx)


## Citation

If you use BactoID in your research, please cite:

**Średnicka P. BactoID v0.5.0 – containerized ONT full-length 16S bacterial identification workflow. Zenodo.
https://doi.org/10.5281/zenodo.23076229**

BactoID relies on several third-party tools and resources. If you use BactoID in scientific work, please also cite the corresponding software where appropriate.

- **Snakemake**  
  Mölder F, Jablonski KP, Letcher B, et al. Sustainable data analysis with Snakemake. *F1000Research*. 2021;10:33.

- **NGSpeciesID**  
  Sahlin K, Lim MCW, Prost S. NGSpeciesID: DNA barcode and amplicon consensus generation from long-read sequencing data. *Ecology and Evolution*. 2021;11:1392–1398.  
  https://doi.org/10.1002/ece3.7146

- **Racon**  
  Vaser R, Sović I, Nagarajan N, Šikić M. Fast and accurate de novo genome assembly from long uncorrected reads. *Genome Research*. 2017;27(5):737–746.  
  https://doi.org/10.1101/gr.214270.116

- **NCBI BLAST+**  
  Camacho C, Coulouris G, Avagyan V, et al. BLAST+: architecture and applications. *BMC Bioinformatics*. 2009;10:421.  
  https://doi.org/10.1186/1471-2105-10-421

- **Singularity / Apptainer**  
  Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. *PLoS ONE*. 2017;12(5):e0177459.  
  https://doi.org/10.1371/journal.pone.0177459


