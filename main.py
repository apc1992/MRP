import os
import argparse
from Bio import SeqIO
import time
import shutil
import csv
import pandas as pd

parser = argparse.ArgumentParser(description='The tool will generate n-mers library and align into human proteome', add_help=True)
parser.add_argument('-i', '--input', help="input fasta file", dest="INPUT")
parser.add_argument('-min_length', help="minimum peptide length", dest="MIN", default=8)
parser.add_argument('-max_length', help="maximum peptide length", dest="MAX", default=11)
parser.add_argument('-n_alignments', help="number of alignments retrieved for each peptide", dest="NALIGN", default=5)
parser.add_argument('-r', '--reference', help="human proteome", dest="REF")
parser.add_argument('--full', help='Display full results', action='store_true')
args = parser.parse_args()


def check_fasta(file):
	"""
	Check if a file is in FASTA format
	:param file: input file to check
	:return: boolean
	"""
	print("Checking if input file is in FASTA format...")
	with open(file) as check:
		entries = SeqIO.parse(check, "fasta")
		return any(entries)


def check_peps(file):
	"""
	Check if a file is in csv/tsv format and contains pep column
	:param file: input file to check
	:return: boolean
	"""
	print("Checking format for peptides file...")
	# Check header of the file #
	with open(file, newline='') as f:
		header = f.readline()
		if "pep" not in header.lower():
			print("Error: Peptide column name should contain the string 'pep'")
			return False
	with open(file, newline='') as f:
		# Try to read the first few lines of the file using the csv module
		try:
			dialect = csv.Sniffer().sniff(f.read(1024))
		except csv.Error:
			print("Only one column detected")
			return 'single'
		# Return the format based on the delimiter of the dialect
		if dialect.delimiter == ',':
			print("Valid csv format for peptide file")
			return 'csv'
		elif dialect.delimiter == '\t':
			print('Valid tsv format for peptide file')
			return 'tsv'
		else:
			print("Peptide file not in csv or tsv format")
			return None


def prepare_database(fasta):
	print("Generating database for input reference...")
	cwd = os.getcwd()
	filename = os.path.basename(fasta)
	db_path = os.path.join(cwd, "db")
	if not os.path.exists(db_path):
		os.mkdir(db_path)
	fastafile = os.path.join(db_path, filename)
	shutil.copy(fasta, fastafile)
	instr = ["makeblastdb", "-in", fastafile, "-dbtype", "prot"]
	instr = " ".join(instr)
	os.system(instr)
	os.remove(fastafile)
	print("Database generated in " + db_path)
	index = os.path.join(db_path, filename)
	return index


def generate_nmers(fasta):
	print("Generating peptides from length " + str(args.MIN) + " to " + str(args.MAX))
	filename = os.path.basename(fasta)
	filename = filename.split('.fa')[0]
	print("Writing peptides as {}_nmers.txt".format(filename))
	with open('{}_nmers.txt'.format(filename), 'w') as out:
		out.write("peptide,id,length\n")
		for l in range(args.MIN, args.MAX + 1):
			with open(fasta) as inp:
				for line in inp:
					line = line.rstrip()
					if line.startswith('>'):
						id = line.replace('>', '')
					if not line.startswith('>'):
						for i in range(len(line)):
							seq = line[i:i + l]
							if len(seq) == l:
								out.write(seq + "," + id + "," + str(l) + "\n")
	return os.path.abspath("{}_nmers.txt".format(filename))


def align_nmers(row, db):

	## WHAT ABOUT USING BLASTP-short instead of blastp
	cwd = os.getcwd()
	align_path = os.path.join(cwd, "alignments")
	n_align = args.NALIGN
	if not os.path.exists(align_path):
		os.mkdir(align_path)
	outfile = os.path.join(cwd, align_path, row.peptide + "_alignments.txt")
	pep = '"' + row.peptide + '"'
	blast_args = ["blastp", "-db", db, "-query", "-",
				"-outfmt", "'6 sseqid sstart send pident score gaps'",
				"-gapopen", "32767", "-gapextend", "32767",
				"-evalue", "1e6", "-max_hsps", "1", "-matrix", "BLOSUM62",
				"-num_alignments", str(n_align), "-out", outfile]
	blast_args = " ".join(blast_args)
	cmd = "echo " + pep + "|" + blast_args
	os.system(cmd)
	try:
		alignments = pd.read_csv(outfile, sep="\t", header=None)
		alignments = alignments.rename(columns={0: 'ref_id', 1: 'start_pos', 2: 'end_pos', 3: 'identity', 4: 'score', 5: 'n_gaps'})
		best_match = alignments.loc[alignments['identity'].idxmax()]
		out = best_match['ref_id'] + ";" + str(best_match['identity']) + ";" + str(best_match['start_pos']) + ";" + str(best_match['end_pos']) + ";" + str(best_match['score']) + ";" + str(best_match['n_gaps'])
	except pd.errors.EmptyDataError:
		out = "none;0;-;-;none;-"
	return out


def main():
	start_time = time.time()
	cwd = os.getcwd()
	filename = os.path.basename(args.INPUT)
	filename = filename.split('.')[0]
	# Creating database from the reference provided
	if args.REF:
		db = prepare_database(args.REF)
	# check input file format
	if check_fasta(args.INPUT):
		print("Valid FASTA format!")
		peps = generate_nmers(args.INPUT)
		df = pd.read_csv(peps)
	else:
		print("Not in FASTA format!")
		fmt = check_peps(args.INPUT)
		if fmt == 'single' or fmt == 'csv':
			df = pd.read_csv(args.INPUT)
		elif fmt == 'tsv':
			df = pd.read_csv(args.INPUT, sep='\t')
	# Process the peptides
	n_peps = len(df)
	print(str(n_peps) + " peptides before filtering!")
	print("Aligning peptides into the reference...")
	df["Alignments"] = df.apply(lambda row: align_nmers(row, db), axis=1)
	df[["ref_id", "identity", "start_pos", "end_pos", "score", "n_gaps"]] = df["Alignments"].str.split(';', expand=True)
	del df["Alignments"]
	if not args.full:
		shutil.rmtree(os.path.join(cwd, "alignments"))
		shutil.rmtree(os.path.join(cwd, "db"))
	if args.full:
		df.to_csv('{}_full_alignments.csv'.format(filename))
	df_filtered = df[df['identity'] != "100.0"]
	df_filtered.to_csv('{}_filtered_alignments.csv'.format(filename))
	n_peps = len(df_filtered)
	print(str(n_peps) + " peptides after filtering!")
	end_time = time.time()
	total_time = end_time - start_time
	print(f"Total running time: {total_time:.2f} seconds")


if __name__ == "__main__":
	main()
