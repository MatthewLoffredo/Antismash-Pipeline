import csv
import os
import urllib

# parse blast input
def parse_blast(filename, headers):
  x=[]
  blast_results=open(filename,'r')
  rows=csv.DictReader(blast_results,headers,delimiter=',')
  for row in rows:
    x.append(row)
  blast_results.close()
  return x

#fetch all accessions passed in and download to ./sequences
def fetch_records(Entrez, ids):

  #make a sequences directory to store all sequence fasta files in (if there isn't already one)
  if not (os.path.isdir('sequences')):
    os.system('mkdir sequences')

  #make dict of each accession and it's retreived record to retrun
  accessionDict = {}

  #this for loop takes each id and downloads it's assembly file from NCBI through FTP
  idsused = 0 #count how many files have been downloaded so the user can follow the progress of the script
  for id in ids:
    
    idsused+=1
    print("Downloading " + id + " from NCBI... (" + str(idsused) + ")")
    
    #search NCBI assembly database by genome accession # (id), and retrieve only 1 record
    handle = Entrez.esearch(db="assembly", term=id+"[id]", retmax='1')
    record = Entrez.read(handle)
    
    #get the assembly search id for the genome accession provided
    searchid=record['IdList'][0]
    
    #get the assembly summary/report using the assembly search id
    esummary_handle=Entrez.esummary(db="assembly",id=searchid, report="full")
    esummary_record = Entrez.read(esummary_handle,validate=False)
    
    #get the FTP link for the assembly file
    try:
      url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
    except IndexError:
      continue
    
    
    #download the assembly file via FTP to the sequences directory
    label = os.path.basename(url)
    if not os.path.isfile('sequences/' + f'{label}.fna'):
      link = os.path.join(url,label+'_genomic.fna.gz')
      link=link.replace('\\','/') #reformat some slashes from the link
      urllib.request.urlretrieve(link, 'sequences/' + f'{label}.fna.gz')
    else:
      print(id + " found in sequences directory.")

    #add id and sequence name to dict
    accessionDict[f'{label}.fna'] = id
  return accessionDict