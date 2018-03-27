'''
Written by Kristen Wells

Takes a fastq file and aligns to the mouse genome (including ERCCs) using STAR

The STAR output is pipped into FeatureCounts to count the number of times
each gene appears.

Depends on subprocess, sys, os, collections, sys, getopt
'''

import subprocess
import os
import copy
from collections import defaultdict
import sys
import getopt

def main():
	'''
	Uses getopt to read in options from the command line
	'''
	try:
		opts, args = getopt.getopt(sys.argv[1:], "ha:c:f:p:g:d:i:", ["align=", "count=",
			"file=", "output=", "GTF=", "index_dir=", "input=", "help",
			"output_dir=", "paired_end="])
	except getopt.GetoptError as err:
		print(err)
		usage()
		sys.exit(2)

	# set default values
	count_files = True
	align_files = True
	make_one_file = True
	paired_end = False
	output_dir = os.getcwd() + "/"
	output_file = "counts.txt"

	opt_dict = dict(opts)

	# Check to ensure all necessary arguments are present
	if not ('-i' in opt_dict or '--input' in opt_dict):
		usage()
		sys.exit(2)
	elif not ('-g' in opt_dict or '--GTF' in opt_dict):
		usage()
		sys.exit(2)

	elif not ('-d' in opt_dict or '--index_dir' in opt_dict):
		usage()
		sys.exit(2)

	# Go through each argument and set appropriate values
	for o, a in opts:
		if o in ("-h", "--help"):
			usage()
			sys.exit(2)
		elif o in ("-a", "--align"):
			align_files = (a=='True')
		elif o in ("-c" , "--count"):
			count_files = (a=='True')
		elif o in ("-f", "--file"):
			make_one_file = (a=="True")
		elif o == "--output":
			output_file = a
		elif o in ("-g", "--GTF"):
			gtf_path = a
		elif o in ("-d", "--index_dir"):
			index_dir = a
		elif o in ("-i", "--input"):
			input_files = a
		elif o == "--output_dir":
			output_dir = a
		elif o in ("-p" or "--paired_end"):
			paired_end = (a=="True")
		else:
			usage()
			sys.exit(2)

	# Return all values
	return([align_files, count_files, make_one_file, paired_end, output_file,
		gtf_path, index_dir, input_files, output_dir])

def usage():
	'''
	Prints out the usage statement for the script
	'''

	usage = """
	-h     --help        Get usage informtion

=============================================================================================================
	Required arguments:
	-g     --GTF         GTF file path
	-d     --index_dir   Index directory path
	-i     --input       Name of input file with sample names and file locations (tab delimited, 2 columns)


=============================================================================================================
    Optional arguments
	-a    --align        Align reads (True/False)                default = True
	-c    --count        Count reads (True/False)                default = True
	-f    --file         Make alignment file (True/False)        default = True
	-p    --paired_end   Input fastq is paired end (True/False)  default = False
	--output             Name of output file                     default = counts.txt
	--output_dir         Output file directory                   default = working directory
	"""
	print(usage)

def align_reads(fastq_input, outfile_prefix):
	'''
	Calls STAR aligner
	'''

	star_list = ["STAR", "--runThreadN", "10", "--genomeDir",
		index_dir, " --sjdbGTFfile", gtf_path,
		"--readFilesIn"] + fastq_input + ["--readFilesCommand",
		"zcat", "--outSAMtype", "BAM", "SortedByCoordinate",
		"--outFileNamePrefix", outfile_prefix]

	star_call = (star_list)

	subprocess.call(star_call)

def featCount(input_list):
	'''
	Calls featureCounts
	'''

	featureCounts_list = ["featureCounts","-a", gtf_path, "-o", "countsOutput",
		"-s", "0", "-R", "CORE", "--primary"]

	featureCounts_list.extend(input_list)
	
	featcounts_call = (featureCounts_list)

	subprocess.call(featcounts_call)


# As long as this is the main script, read from the command line
if __name__ == "__main__":
	(align_files, count_files, make_one_file, paired_end, output_file,
		gtf_path, index_dir, input_file_sheet, output_dir) = main()

# initialize lists and dictionaries
bam_list = []
featCount_dict = {}
gene_dict = {}
sample_dict = {}
gene_list = []
sample_list = []
header = ["chr", "start", "end", "name", "strand"]
aligned_files = output_dir + "aligned_files.txt"

# Only align reads if align_file = True
if align_files:

	print("############")
	print("# Aligning #")
	print("############")
	print()

	# Open an output file to hold the paths to each aligned file
	with open(aligned_files, "w") as file_output:

		# Open the file containing paths to each fastq file
		with open(input_file_sheet, "r") as sample_sheet:

			# Read each line of input sheet (read each file path and sample name)
			for line in sample_sheet:
				line = line.strip().split("\t")

				# Name the fastq file and sample name based on the input sheet
				if paired_end:
					fastq_file = line[0:2]
					sample_name = line[2]

				else:
					fastq_file = [line[0]]
					sample_name = line[1]

				#name star output file
				bam_file = (output_dir + sample_name + "Aligned.sortedByCoord.out.bam")
			
				#call the star aligner
				align_reads(fastq_file, sample_name)

				# write the full file path to the output file
				file_output.write(bam_file + "\t" + sample_name + "\n")

# Open the file containing paths to aligned files
with open(aligned_files, "r") as file_names:
	for line in file_names:
		bam_file, sample_name = line.strip().split("\t")

		# Add each file path to the bam list
		bam_list.append(bam_file)

		# Add bam file and sample name to dictionary (this is to keep the two easily linked)
		featCount_dict[bam_file] = sample_name

		# Add the sample name to the list of samples. This is for the final output file
		sample_list.append(sample_name)

# Only call feature counts if count_gene = True
if count_files:

	print("##################")
	print("# Counting reads #")
	print("##################")
	print()

	# call featureCounts
	featCount(bam_list)

# Only make the final file if make_one_file = True
if make_one_file:

	print("#####################")
	print("# Making count file #")
	print("#####################")
	print()

	# Open the GTF file and make a dictionary of all the genes in the file
	with open(gtf_path, "r") as gene_reads:
		for line in gene_reads:
			line = line.strip().split("\t")

			# Lazy way of only including lines that are long enough
			if len(line) > 2 and (line[2] == "gene" or line[1] == "ERCC"):

				# Name each part of the line of the GTF file, you have to trick it a bit
				chromosome = line[0]
				start = line[3]
				end = line[4]
				strand = line[6]
				gene_info = line[8].split(" ")
				gene_name = gene_info[1]
				gene_name = gene_name[1:len(gene_name) - 2]

				# Add all information to a dictionary. This makes it so you can easily
				# Find any information about a gene given its name
				gene_dict[gene_name] = [chromosome, start, end, gene_name, strand]

				# Also make a list of genes. This is necessary to keep ordering consistent
				gene_list.append(gene_name)

	# Go though each file in your bam list. This is every STAR output file.
	for bam_file in bam_list:

		# Name the file to match the featureCounts naming scheeme
		file_name = bam_file + ".featureCounts"

		# Now that the file is properly named, open the file
		with open(file_name, "r") as count_file:

			# Print so you can keep track of where you are
			print(file_name)

			# Make a dictionary of counts. This is an empty dictionary
			count_dict = defaultdict(int)

			# read each line of the featureCount file
			for line in count_file:
				line = line.strip().split("\t")

				# count each gene in featureCount file
				if line[2] == "1":

					# The gene name is the fourth column
					gene = line[3]
					count_dict[gene] += 1

		# Name the sample
		sample_id = featCount_dict[bam_file]

		# Add the dictionary for the sample into a new dictionary (this is a dictionary
			# of dictionaries)
		sample_dict[sample_id] = count_dict


	# Open the output file (this is currently an empty file)
	with open(output_file, "w") as final_output:

		# Make the header. This includes naming each column by the sample ids
		header.extend(sample_list)
		final_output.write("\t".join(header) + "\n")

		# Go through each gene in the order of the GTF file
		for gene in gene_list:

			# Write all of the information for that gene for each sample in the output file
			# This includes name, start, stop, strand, ect.)
			final_output.write("\t".join(gene_dict[gene]))

			# Now write the count for each sample. Use the sample_list to ensure the order
			# is the same as the header
			for sample in sample_list:
				final_output.write("\t" + str(sample_dict[sample][gene]))

			# Write a newline to the output file.
			final_output.write("\n")
