import os
import argparse
from utils import fetch_records
from Bio import SeqIO
from Bio import Entrez
import gzip
import shutil
import csv

# create the command-line parameters for email (for Entrez) and input file
parser = argparse.ArgumentParser()
parser.add_argument(
    '-e', '--email', help="Enter your email for Entrez here", required=True)
parser.add_argument(
    '-i', '--input', help="Input file name (list of NCBI accession numbers, one on each line)")
args = parser.parse_args()

# set Entrez's email to the email provided
Entrez.email = args.email

# read the genome accession #'s from the input file
# (specified in args.input parameter)
ids = []
if args.input != None:
    for line in open(args.input).readlines():
        ids.append(line.strip())

# get current working directory
PATH = os.getcwd()

# fetch all accessions passed in, returns dict of accessions mapped to their fetched files
accessionDict = fetch_records(Entrez, ids)

# make antismash output directory if not already created
if not (os.path.isdir('antismash_output')):
    os.system('mkdir antismash_output')
outpath = f'{PATH}/antismash_output'

# list to store all secondary metabolite results from antismash
allMetabolites = {}

# count the files ran so the user can follow the progress of the script
fileCount = 0

# this for loop goes through each file in the sequences directory, unzips it, and runs antismash on it
for fileName in os.listdir("sequences"):

    # if the file is zipped, unzip it
    if fileName.endswith('.gz'):
        unzippedFileName = fileName[:fileName.index('.gz')]
        # unzip the .gz assembly files
        with gzip.open("sequences/" + fileName, 'rb') as f_in:
            with open("sequences/" + unzippedFileName, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        # delete the original zipped file
        os.system('rm sequences/' + fileName)

    # if the file is already unzipped, set its name to the unzipped variable
    else:
        unzippedFileName = fileName

    fileCount += 1
    print("Running " + unzippedFileName + "..." + "(" + str(fileCount) + ")")

    infile = f'{PATH}/sequences/{unzippedFileName}'

    # run subclusters, knownclusters, active site finder analysis
    options = '--cb-subclusters --cb-knownclusters --asf --genefinding-tool prodigal'

    # run antismash on genome if not run already
    if not os.path.isdir(f'antismash_output/{os.path.splitext(unzippedFileName)[0]}'):
        antismashCmd = f'sh run_antismash.sh {infile} {outpath} {options}'
        print(antismashCmd)
        os.system(antismashCmd)

    # parse secondary metabolites from output and push to master list
    # multiple gbk files are created, skip the first one and parse the rest
    metabolites = []
    output = f'{PATH}/antismash_output/{os.path.splitext(unzippedFileName)[0]}'
    loopNum = 0
    for outFile in os.listdir(output):
        if outFile.endswith('gbk'):
            # skip first gbk file
            if loopNum != 0:
                print('parsing ' + outFile + " for secondary metabolites...")
                record = SeqIO.read(f'{output}/{outFile}', "genbank")
                feats = [
                    feat for feat in record.features if feat.type == "protocluster"]
                for feat in feats:
                    metabolites.append(feat.qualifiers['product'][0])
            else:
                loopNum = loopNum + 1
                continue
        else:
            continue

    allMetabolites[unzippedFileName] = metabolites

# create a dictionary to store each accession and it's metabolites
metaboliteDict = {accessionDict[fileName]: metabolites for (
    fileName, metabolites) in allMetabolites.items()}
print("\nResults:")
print(metaboliteDict)

# #open the output file
csv_output_name = "AntismashResults.csv"
csv_output = open(csv_output_name, 'w')

# write the headers to the output file
csv_output.write('Genome Accession,Secondary Metabolites\n')

# write results to file
for k, v in metaboliteDict.items():
    csv_output.write(k + ",")
    csv_output.write(" ".join(v) + "\n")

exit()