
![Spectre](./logo.png)
# Spectre - Long read CNV caller

## Installation
Recommended ways of installation are using pip or conda or pip package:
```
pip install ont-spectre
```
(or)
```
conda install nanoporetech::ont-spectre
```

## Run from sources

Setup a conda environment for Spectre (copy and paste the following commands)

```bash
conda create -n spectre python=3.8.5 pip -y
conda activate spectre
pip install -r requirements.txt
```

## How to run
Spectre need as input:
- The result of Mosdepth (directory)
- Reference genome (can be bgzip compressed)
- Window size used in Mosdepth (Make sure the binsize between Mosdepth and Spectre are matching. We suggest a binsize of 1000 base pairs.)
- VCF file containing SNV

>INFO: Make sure to include "/" at the end if you are adding directory paths.


```bash
spectre CNVCaller \
  --bin-size 1000 \
  --coverage mosdepth/sampleid/ \
  --sample-id sampleid \
  --output-dir sampleid_output_directory_path/ \
  --reference reference.fasta.gz \
  --snv sampleid.vcf.gz
```