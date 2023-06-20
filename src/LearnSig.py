#!/usr/bin/env python3

# ------------------------------
#  python modules
# ------------------------------

import os,sys
import pandas as pd
import numpy as np
import subprocess
from scipy import stats

# ------------------------------
#  own python modules
# ------------------------------

import CondSig
from CondSig.BasicSetting import *
from CondSig import Utility

# ------------------------------
#  Misc functions
# ------------------------------

class SignatureLearning():
	"""Learn signature using biterm topic model"""

	def __init__(self, df_peakov_raw, df_peakov_filtered, args, output_prefix):
		self.df_peakov_raw = df_peakov_raw
		self.df_peakov_filtered = df_peakov_filtered
		self.args = args

		self.output_prefix = output_prefix
		

	def learn_focus(self, focus_CAP, count, total_jobs):
		"""Learn signature on each focus context"""
		
		info("2-{0}. Learn co-occupancy signatures on focus context of {1} ({0}/{2})... " . format(count, focus_CAP, total_jobs))

		
		info("2-{0}-1. Select focus context of {1} ..." . format(count, focus_CAP))
		
		# calculate co-occupancy score to evaluate co-occupancy level between focus CAP and all CAPs
		df_gF1s = Utility.GF1Calculation(self.df_peakov_filtered, self.args).cal_focus(focus_CAP)

		# determine potential combinatorial CAPs of focus CAP (including focus CAP)
		CAP_remained_number = int(self.args.focus_number)
		potentail_combinatorial_CAPs = df_gF1s.loc[focus_CAP,:].sort_values(ascending = False).index.values[:CAP_remained_number]
		df_peakov_filtered_pcCAPs = self.df_peakov_filtered.loc[:,potentail_combinatorial_CAPs]

		# determine focus context and control context
		df_peakov_filtered_transformed = df_peakov_filtered_pcCAPs * stats.zscore(df_gF1s.loc[focus_CAP,potentail_combinatorial_CAPs])
		focus_bins = df_peakov_filtered_transformed.loc[df_peakov_filtered_transformed.sum(axis = 1) > 0,:].index.values
		# control_bins = df_peakov_filtered_transformed.loc[df_peakov_filtered_transformed.sum(axis = 1) <= 0,:].index.values

		df_peakov_filtered_focus = self.df_peakov_filtered.loc[focus_bins, potentail_combinatorial_CAPs]
		df_peakov_filtered_focus = df_peakov_filtered_focus.loc[df_peakov_filtered_focus.sum(axis = 1) >= 2,:]

		# build biterm topic model on focus context
		df_valid_promoter_topics, df_valid_nonpromoter_topics = Utility.BuildBTM(focus_CAP, count, self.df_peakov_raw, self.df_peakov_filtered, df_peakov_filtered_focus, self.args, self.output_prefix).run()

		return(df_valid_promoter_topics, df_valid_nonpromoter_topics)


class DataProcessing():
	"""Process data provided by annotation file to generate genome-wide occupancy
	and select focus context for each CAP iteratively
	"""

	def __init__(self, 
				args
				):

		self.args = args
		self.df_dataset = pd.read_csv(self.args.data_annotation, header = None, sep = "\t").iloc[:, 0:4]
		self.df_dataset.columns = ["factor", "label", "file", "uniprot_id"]
		generate_abspath = lambda x:os.path.join(self.args.data_dir, x) 
		self.df_dataset.loc[:, "file"] = list(map(generate_abspath, self.df_dataset.loc[:, "file"].values))
		uni_uniprot_id, count_uniprot_id = np.unique(self.df_dataset.loc[:, "uniprot_id"].values, return_counts = True)
		if sum(count_uniprot_id > 1) > 0:
			error("More than 1 CAPs are assigned to the same uniprot id ({0}), please check it and don't use redudant data for each CAP." . format(uni_uniprot_id[count_uniprot_id > 1]))
			sys.exit(1)

		self.chromsize_file = self.args.genome_annotation.chromsize_file


	def run(self):

		os.chdir(self.args.out_dir) # change work directory to the output directory

		info("1. Generate genome-wide occupancy matrix for all CAPs ...")
		
		if self.df_dataset.shape[0] < 50:
			error("Too few CAPs input, the number of CAPs should be >= 50.")
			sys.exit(1)
		
		genome_bins_file = self.split_genome_bins() # genome was divided into consecutive bins
		output_prefix = "{0}_1kb_bins" . format(self.args.name)
		genome_occupancy_file = Utility.OccupancyMatrix(genome_bins_file, self.df_dataset, output_prefix, self.args).intervals_intersect_genome_bins()
		df_peakov_raw, df_peakov_filtered = self.filter_occupancy_matrix(genome_occupancy_file)
		return(df_peakov_raw, df_peakov_filtered, output_prefix)

	def split_genome_bins(self, res = 1000):
		"""split genome to consecutive bins with the given resolution"""
		
		df_gsize = pd.read_csv(self.chromsize_file, header = None, sep = "\t")
		genome_bins_file = "{0}_1kb_bins.bed" . format(self.args.genome_version)
		
		with open("{0}_1kb_bins_all.bed" . format(self.args.genome_version), "w") as outf:
			count = 1
			for index, row in df_gsize.iterrows():
				chrom = row[0]
				size = row[1]
				left = 0
				right = int(res)
				while right < size:
					outf.write("{0}\t{1}\t{2}\tbin{3}\n".format(chrom, left, right, count))
					left += int(res)
					right += int(res)
					count += 1
				outf.write("{0}\t{1}\t{2}\tbin{3}\n".format(chrom, left, size, count))
				count += 1

		cmd_remove_blacklist = "bedtools intersect -wa -v -e -f 0.5 -F 0.5 -a {0} -b {1} > {2}" . format(
			"{0}_1kb_bins_all.bed" . format(self.args.genome_version),
			self.args.genome_annotation.blacklist_file,
			genome_bins_file
			)
		subprocess.run(cmd_remove_blacklist, shell = True, check = True)

		return(genome_bins_file)


	def filter_occupancy_matrix(self, genome_occupancy_file):
		"""filter occupancy matrix to remove TR with too few occupancy events and bins with too few or too many occupancy events"""
		
		df_peakov_raw = pd.read_csv(genome_occupancy_file, header = 0, sep = "\t", index_col = "name")
		df_peakov = df_peakov_raw.iloc[:, 3:] # remove the columns contains chrom, start and end

		# remove features with fewer than 500 occupancy events
		feature_sum = df_peakov.sum(axis = 0)
		features_removed = feature_sum[feature_sum <= 500].index
		df_peakov_filtered = df_peakov.loc[:, np.setdiff1d(df_peakov.columns, features_removed)]

		# remove bins with too few or too many co-occupancy events
		df_peakov_filtered = df_peakov_filtered.loc[((df_peakov_filtered.sum(axis = 1) >= 2) & (df_peakov_filtered.sum(axis = 1) <= (0.9 * df_peakov_filtered.shape[1]))), :]
		df_peakov_filtered.to_csv("{0}_1kb_bins_peakov_filtered.bed" . format(self.args.name), header = True, sep = "\t", index = True)

		return(df_peakov_raw, df_peakov_filtered)

