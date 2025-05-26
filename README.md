# BED file based NGS off-target calculation script
## Bam and BED file based usage


The script use BAM sorted file along with the 4 column format of target BED file as input files.

- Developed in Python3
- 2 input files and 2 other input required
- Calculate target region, extended region and off-target region reads and bases
- Provide top 10 chr based stats as additional data
- The read length data is required to calculate extended panel info and provided in -r option

## Required packages and Installation

Installation of Python packages as provided below :

```sh
pip install pysam pandas numpy intervaltree
```

Before using the provided script - Confirm Python3 and Pip was installed in the machine


## BED file format for input usage

Below is the model of BED file format. The file should be tab sep format with 4 columns.
```sh
chr1    1000    2000    target1
chr1    3000    4000    target2
chr2    5000    6000    target3
```

## Usage of script

### Basic usage
> python off-target-checker_v2.py -b sample.bam -t targets.bed -r 150

### With verbose output to see progress
> python off-target-checker_v2.py -b sample.bam -t targets.bed -r 150 -v

### Save results
> python off-target-checker_v2.py -b sample.bam -t targets.bed -r 150 -o report.txt

### Help

> python off-target-checker_v2.py --help

## Output 

Please check the ouput.txt file for the model ouptut
