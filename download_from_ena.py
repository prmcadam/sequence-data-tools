"""
Download paired end Illumina data from ENA (http://www.ebi.ac.uk/ena/) with built in error checking.
Usage python download_from_ena.py -e [ENA ACCESSION] -d [DIRECTORY TO DOWNLOAD TO]

Output directory must already exist.
Script:
1. Connects to ENA
2. Retrieves data for given accession number (run_accession,fastq_md5,fastq_ftp)
3. Checks if .fastq.gz files with same names exist already in output directory
4. Compares md5sums if file already exists
5. If md5sums do not match, file on disk is removed
6. Downloads .fastq.gz if file does not exist, or md5sums do not match
7. Repeats from step 3 to check all are downloaded correctly


Dependencies
wget
gnu parallel
python libraries are part of core python distribution

#To do
#Download progress status
#Possibly works with single end files, need to test properly
#Test with downloading data other than Illumina

"""

import hashlib
import os
import urllib2
from optparse import OptionParser
import sys

def opts():
	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	parser.add_option("-d", "--directory", dest="directory", action="store", help="path to directory of .fastq.gz files", metavar="FILE")
	parser.add_option("-e", "--ENA", action="store", dest="accession", help="ENA project ID.  Can be comma separated list of accessions, or text file with one accession per line")
	return parser.parse_args()

def fetch_ENA_data(accession):
	"""
	Fetch data from ENA for project of interest
	"""
	print "Fetching data from ENA for "+accession
	url='http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession='+str(accession)+'&result=read_run&fields=run_accession,fastq_md5,fastq_ftp&download=text'
#	url='http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession='+str(accession)+'&result=read_run&fields=run_accession,submitted_md5,submitted_ftp&download=text'
	req=urllib2.Request(url)
	ENA_data=urllib2.urlopen(req).read()
	while ENA_data.startswith('Timed out'):
		ENA_data=urllib2.urlopen(req).read()
#	print ENA_data
#	if ENA_data.startswith('Timed out'):
#		sys.exit('!!! Connection to ENA timed out - Try again later !!!')
	return ENA_data.strip().split('\n')
	
def return_md5_dict(ENA_accession):
	"""
	Return dictionary of fastq=[ftp path, expected md5sum]
	"""
	md5_dict={}
	lines=ENA_accession
	if lines[0].startswith('run_accession	fastq_md5	fastq_ftp') != True:
#	if lines[0].startswith('run_accession	submitted_md5	submitted_ftp') != True:
		return False
	else:
		for line in lines[1:]:
			try:
				splitline=line.split()
				run_accession=splitline[0].strip()
				if ';' in line:
					fwd_fastq=splitline[2].split(';')[-2].split('/')[-1].strip()
					rev_fastq=splitline[2].split(';')[-1].split('/')[-1].strip()
					fwd_path=splitline[2].split(';')[-2]
					rev_path=splitline[2].split(';')[-1]
					fwd_md5=splitline[1].split(';')[-2].strip()
					rev_md5=splitline[1].split(';')[-1].strip()
					md5_dict[fwd_fastq]=[fwd_path,fwd_md5]
					md5_dict[rev_fastq]=[rev_path,rev_md5]
				else:
					fwd_fastq=splitline[2].split(';')[0].split('/')[-1].strip()
					fwd_path=splitline[2]
					fwd_md5=splitline[1].strip()
					md5_dict[fwd_fastq]=[fwd_path,fwd_md5]
			except IndexError:
				pass
		return md5_dict

def check_md5(fastq, block_size=2**20):
	"""
	Return md5sum of file on disk
	"""
	f=file(fastq)
	md5 = hashlib.md5()
	while True:
		data = f.read(block_size)
		if not data:
			break
		md5.update(data)
	f.close()
	return md5.hexdigest()

def file_exists(fastq):
	"""
	Check if file exists.  ie. has been downloaded from ENA
	"""
	return os.path.isfile(fastq)
	
def download_missing_files(accession, missing_fastq):
	"""
	Download any missing data from ENA
	"""
	with open('temp_download.txt','w') as out_file:
		out_file.write('run_accession	fastq_md5	fastq_ftp\n')
#		out_file.write('run_accession	submitted_md5	submitted_ftp\n')
		for fastq in missing_fastq:
			out_file.write(accession+'\t'+fastq[1]+'\t'+fastq[0]+'\n')
	sys.stdout.write("\nDownloading "+str(len(missing_fastq))+" fastq files from SRA\n")
	os.system('cat temp_download.txt | cut -f3 | parallel -P 8 wget {} > /dev/null 2>&1')
	return

def fastq_status(fastq,exists,md5_match,correct,missing,md5_dict):
	if exists and md5_match:
		correct.append(fastq)
	elif exists and md5_match != True:
		missing.append(md5_dict[fastq])
		sys.stdout.write("\nRemoving: "+fastq+', incorrect md5sum\n')
		os.remove(fastq)
	elif exists != True:
		missing.append(md5_dict[fastq])
	return correct, missing
		

##########################################################################################

if __name__ == "__main__":
	(options,args) = opts()
	test_accession=options.accession.strip()
	if os.path.isfile(test_accession):
		with open(test_accession,'r') as f:
			accessions=[x.strip() for x in f]
	else:
		accessions=test_accession.split(',')
		
	directory=options.directory
	os.chdir(directory)

	for accession in accessions:
		ENA_data=fetch_ENA_data(accession)
		md5_dict=return_md5_dict(ENA_data)
		if md5_dict:
			count=0
			correct,absent=[],[]
			print 'Fetched '+str(len(md5_dict))+' fastq records'
			for fastq in md5_dict:
				count+=1
				sys.stdout.write("\rProcessing: "+str(count)+'/'+str(len(md5_dict)))
				sys.stdout.flush()
				exists=file_exists(fastq)
				if exists:
					file_md5=check_md5(fastq)
					md5_match = file_md5==md5_dict[fastq][1]
					correct,absent=fastq_status(fastq, exists,md5_match,correct,absent,md5_dict)
				else:
					absent.append(md5_dict[fastq])
			download_missing_files(accession, absent)
		else:
			print 'Incorrect input ENA format\n'
			sys.exit()
	
				
		with open('temp_download.txt','r') as new_downloads:
			lines=new_downloads.readlines()
			md5_dict=return_md5_dict(lines)
			count=0
			correct,absent=[],[]
			if len(md5_dict)==0:
				print 'All files successfully downloaded\n'
			else:
				print '\nSecond attempt\nFetched '+str(len(md5_dict))+' fastq records'
				for fastq in md5_dict:
					count+=1
					sys.stdout.write("\rProcessing: "+str(count)+'/'+str(len(md5_dict)))
					sys.stdout.flush()
					exists=file_exists(fastq)
					if exists:
						file_md5=check_md5(fastq)
						md5_match = file_md5==md5_dict[fastq][1]
						correct,absent=fastq_status(fastq,exists,md5_match,correct,absent,md5_dict)
					else:
						absent.append(md5_dict[fastq])
	
				download_missing_files(accession, absent)
