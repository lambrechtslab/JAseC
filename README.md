# JASEC
**J**unction-based **A**lternative **S**plicing **E**vent **C**ounter (JASEC) is a bioinformatic program to detect and quantify alternative splicing events (ASEs) from 10x Genome single-cell RNAseq data. JASEC detects reads spanning exon-exon junctions within a gene, while discovering all types of ASEs without relying entirely on transcript annotation, quantifying ASEs on cell level based on barcodes,  adopting the uneven read coverages, and incorporating unique molecular identifier (UMI) information to handle single-cell data (Figure 1d; see Methods). JASEC only considers ASEs supported by at least 10 junction reads, allowing both a reference and alternative AS isoform to be defined.

## System requirements
JASEC was designed to run under the command line of Linux environment compatible with Julia language and with at least 64G memory.

## Installation
Running JASEC relies on Julia language compiler on Linux OS. First of all, the user needs to install Julia (v1.6+). A binary version of Julia program can be downloaded from the [official website](https://julialang.org/downloads/).

JASEC can be installed by running the below code in bash:
```bash
$ julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/lambrechtslab/JASEC.git"))'
```
Or install under julia package environment (press key `]` in Julia REPL):
```julia
(v1.9) pkg> add https://github.com/lambrechtslab/JASEC.git
(v1.9) pkg> build JASEC
```
The dependent packages and the external tool `StringTie` will be installed automatically. A softlink `~/JASEC` will be created at home directory and you can move it to where you want.

## Uninstallation
To uninstall JASEC, you need to remove the JASEC package in julia and clean the softlink.
```
$ julia -e 'using Pkg; Pkg.rm("JASEC")'
$ unlink ~/JASEC
```

## Usage
### Input data preparation
Before analysis by JASEC, fastq files generated from 10x Genome scRNAseq platform should be mapped by `cellranger`, and optionnaly but highly recommended cell types should be identified. Three types of files are needed for JASEC:
1. The bam files output from `cellranger`.
2. A sample sheet tsv file with two columns: sample ID (unique and no blanks), and bam file path. See example in `test/10x/samplesheet.tsv`.
3. A genome annotation file (gtf). This file is suggested to use the same version as in the previous `cellranger`. For example, the Ensembl v93 human annotation gtf file can be downloaded from [here](https://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz). `.gz` compressed file is accepted.
4. (Optional but recommended) A tsv file for the cell group information, which should be composited by three columns: sample ID (identical to file #2), cell barcode (do not include `-1` at the end of barcode), and name of cell group (like cell type as an example). This file can be in compressed in .gz format. If this file is provided, JASEC will calculate PSI and counts for each cell groups. Considering read coverage on single-cell level is limited, appropriately grouping cells could largely facilitate downstream analysis.
5. (Optional) A de-nova assembled GTF file. If missing, JASEC will generate one based on the input bam files and annotation files using `StringTie`.
### Run JASEC program
```bash
$ ~/JASEC -g ANNOTATION.gtf(.gz) -s SAMPLESHEET.tsv [-p THREAD] \
      	     [-c CELL_BARCODE.tsv(.gz)] [-o OUTPUT-DIR] \
             [--assembled-annotation ASSEMBLED_ANNOTATION.gtf(.gz)] \
	     [--stringtie STRINGTIE_PROGRAM] [--temp-dir TEMP_DIR]
```
If `julia` is not in `$PATH` environment, run the below command instead:
```bash
$ /path/to/julia ~/JASEC ...
```
#### Common arguments:
*  `-g`, `--ref-annotation` _ANNOTATION.gtf_ GTF files of public reference annotation.
*  `-s`, `--samplesheet` _SAMPLESHEET.tsv_
                        A tsv file with columns of sample id and bam files. Sample ID should be unique and without internal blanks.
*  `-c`, `--cell-barcode` _CELL_BARCODE.tsv(.gz)_
                        A tsv(.gz) file with columns of sample id, cell barcode and cell group.
*  `-o`, `--output-dir` _OUTPUT_DIR_ Directory of the outputs. JASEC will create one if not exists. (default: ".")
*  `-p`, `--thread` _THREAD_ Number of threads to run JASEC. An efficient thread number is no larger than the number of input bam files. (default: 1)
#### Other arguments:
*  `--assembled-annotation` _ASSEMBLED_ANNOTATION.gtf_  GTF files of de-noval assembled annotation. If missing, JASEC will assemble it from bam files using StringTie.
*  `--stringtie` _STRINGTIE_CMD_ Command of `StringTie` software. When missing, the internal `StringTie` will be used.
*  `--temp-dir` TEMP_DIR Tempary directionary. (default: "_tmp")
#### Tips
* If `_tmp/` already exists, JASEC will consider it as an interrupted task and try to recover from the last calculation. If you want to recalculate from scratch, please remove this folder first.
## Output files explanation
*  `detected_ASE_info.tsv` Features of each detected ASEs, including gene name, splicing type, inclusion/exclusion junction coordinates, etc. The second exon junction `E_junc_B` is only for mutiple-exclusive-exon (MXE).
*  `[SampleID]_count_per_ASE_perCell.tsv.gz`  The junction read numbers per cell. For MXE, the column `E_junc_count` is the sum of two exclusion junction counts.
*  `[SampleID]_count_perASE_perCellGroup.tsv.gz`  The junction read number and PSI per cell type. This file is only available when `--cell-barcode` _CELL_BARCODE.tsv(.gz)_ was given.
#### Below chart indicates in each of splicing types, which reads are counted as exclusion reads (E) and inclusion reads (I), and how the PSI value is calculated.
![](extdata/chart/bitmap.png)
## Demo
This demo shows how to use JASEC to detect and quantify ASEs on single-cell data of 1000 Peripheral Blood Mononuclear Cells (PBMCs) from a Healthy Donor, which is public available on the website of 10x Genomics as an example. For this data, libraries were generated with 10x Genomics 3' Single Cell Gene Expression v3 chemistry and data were pre-analyzed with Cell Ranger version 3.0 for the aligned bam file and the cell clustering result.

First, download the example bam file of 10x single-cell data (4.5G) into working directory:
```
$ wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_possorted_genome_bam.bam
```
And prepare a sample sheet file:
```
$ echo -e "sample1\tpbmc_1k_v3_possorted_genome_bam.bam" > samplesheet.tsv
```
Next, download human genome annotation file into working directory:
```
$ wget https://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz
$ gunzip Homo_sapiens.GRCh38.93.gtf.gz
```
Then we download the cell clustering result and use the major clusters as cell groups:
```
$ wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_analysis.tar.gz
$ tar zxf pbmc_1k_v3_analysis.tar.gz analysis/clustering/graphclust/clusters.csv
$ cat analysis/clustering/graphclust/clusters.csv | sed 1d | awk -F'-1,' '{print "sample1\t"$1"\tcellgroup"$2}' > cell_barcode.tsv
```
At last, run JASEC. It takes approximately one hour.
```
$ ~/JASEC -g Homo_sapiens.GRCh38.93.gtf -s samplesheet.tsv -c cell_barcode.tsv -o JASEC_output
```
The outputs will be in folder `JASEC_output/`.
## Contacts
Jieyi Xiong (jieyi.xiong[at]kuleuven.be); Diether Lambrechts (diether.lambrechts[at]kuleuven.be)
