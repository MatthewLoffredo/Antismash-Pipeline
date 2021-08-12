# Antismash-Pipeline
Pipeline to find secondary metabolites from antismash in your sequences

Authors: Matt Loffredo

This tool is designed for the user to input a list of accession numbers, which will processed via antismash and output a csv of sencondary metabolites found in the sequences.

# Requirements

* Docker Desktop (https://www.docker.com/products/docker-desktop)
* Antismash standalone docker image (https://hub.docker.com/r/antismash/standalone)
* Python3 and packages used:
  * Biopython (https://biopython.org/wiki/Download)
    * SeqIO
    * Entrez
  * os
  * csv
  * sys
  * urllib
  * gzip
  * shutil
  * argparse

# Instructions

### Docker Desktop

In order to run the tool, Docker Desktop is required, as well as the docker image for antismash. The script executes the image and parses the results. Once you have Docker Desktop installed, run the command
```
docker pull antismash/standalone
```
to pull the required image to your system.

### AntismashScript.py

The input to AntismashScript.py should be a text file with NCBI genome accession numbers (i.e. GCF_001641215.1) listed one per line (see accessions.txt for an example file) using the parameter --input. The script will download the assembly FASTA files for you to a directory called 'sequences' and unzip them. 

**Command for running AntismashScript.py:**
```
python3 AntismashScript.py -e EMAIL -i INPUT_FILE
```

Where EMAIL is your email (used for Entrez), and INPUT is the input file explained above. Your email is a required parameter.

This script may take a long time to run if you use many accession numbers, so be prepared for that.

## Output

A file name AntismashResults.csv will be outputed. The CSV file contains 2 columns, the first being the accession number, the second being the secondary metabolites that were found in that accession. A folder called antismash_output will also be created, which has the antismash output for each accession file.

## Using Test Data:

Test Data (3 accession numbers) is provided in accessions.txt. To run this, run these commands:

```
python3 AntismashScript.py --email EMAIL --input accessions.txt
```