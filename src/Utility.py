#!/usr/bin/env python3

# ------------------------------
#  python modules
# ------------------------------

import pandas as pd
import numpy as np
import subprocess
from scipy import stats
from itertools import combinations
import multiprocessing
from multiprocessing import Pool
import pyBigWig

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.cm import ScalarMappable
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

matplotlib.rcParams['font.sans-serif'] = "Helvetica"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


# ------------------------------
#  own python modules
# ------------------------------

import CondSig
from CondSig.BasicSetting import *

# ------------------------------
#  Misc functions
# ------------------------------

class BuildBTM():
	"""Build biterm topic model on focus context"""

	def __init__(self, focus_TR, count, df_peakov_raw, df_peakov_filtered, df_peakov_filtered_focus, args, output_prefix):
		self.focus_TR = focus_TR
		self.count = count
		self.df_peakov_raw = df_peakov_raw
		self.df_peakov_filtered = df_peakov_filtered # whole context
		self.df_peakov_filtered_focus = df_peakov_filtered_focus # focus context
		self.args = args
		self.promoter_bins = pd.read_csv("../{0}_peakov_promoter.bed" . format(output_prefix), header = None, sep = "\t").iloc[:,3].values
		self.outpath = self.focus_TR

		self.topNword_show = 20

	def run(self):

		self.run_biterm()
		df_filtered_promoter_topics, df_filtered_nonpromoter_topics = self.interpret_biterm()

		return(df_filtered_promoter_topics, df_filtered_nonpromoter_topics)

	def run_biterm(self):
		"""learn topic number on focus context of focus TR"""
		
		if not os.path.isdir(self.outpath):
			os.mkdir(self.outpath)

		# split promoter and non-promoter whole context
		df_peakov_filtered_promoter = self.df_peakov_filtered.loc[self.df_peakov_filtered.index.isin(self.promoter_bins),:]
		df_peakov_filtered_nonpromoter = self.df_peakov_filtered.loc[~self.df_peakov_filtered.index.isin(self.promoter_bins),:]		

		# split promoter and non-promote focus context
		df_peakov_filtered_focus_promoter = self.df_peakov_filtered_focus.loc[self.df_peakov_filtered_focus.index.isin(self.promoter_bins),:]
		df_peakov_filtered_focus_nonpromoter = self.df_peakov_filtered_focus.loc[~self.df_peakov_filtered_focus.index.isin(self.promoter_bins),:]

		# determine the vocab size (CAP count) of promoter and non-promoter context
		vocab_size_promoter = len(df_peakov_filtered_focus_promoter.loc[:, df_peakov_filtered_focus_promoter.sum(axis = 0) >= 1].columns)
		vocab_size_nonpromoter = len(df_peakov_filtered_focus_nonpromoter.loc[:, df_peakov_filtered_focus_nonpromoter.sum(axis = 0) >= 1].columns)

		# build the corpus
		with open("{0}/{1}_promoter_doc.txt" . format(self.outpath, self.focus_TR), "w") as doc_promoter:
			for index,row in df_peakov_filtered_focus_promoter.iterrows():
				doc_promoter.write("\t".join(row[row>0].index) + "\n")
		with open("{0}/{1}_nonpromoter_doc.txt" . format(self.outpath, self.focus_TR), "w") as doc_nonpromoter:
			for index,row in df_peakov_filtered_focus_nonpromoter.iterrows():
				doc_nonpromoter.write("\t".join(row[row>0].index) + "\n")

		# run BTM
		info("2-{0}-2. Run biterm topic model on focus context (promoter) ... " . format(self.count))
		btm_promoter = "biterm_run {0} {1} {2} {3} {4} {5} > {4}/{0}_promoter.log" . format(
			self.focus_TR, self.args.min_signatures, self.args.max_signatures, vocab_size_promoter, self.outpath, "promoter")
		subprocess.run(btm_promoter, shell = True, check = True)

		info("2-{0}-3. Run biterm topic model on focus context (non-promoter) ... " . format(self.count))
		btm_nonpromoter = "biterm_run {0} {1} {2} {3} {4} {5} > {4}/{0}_nonpromoter.log" . format(
			self.focus_TR, self.args.min_signatures, self.args.max_signatures, vocab_size_nonpromoter, self.outpath, "nonpromoter")
		subprocess.run(btm_nonpromoter, shell = True, check = True)

		self.df_peakov_filtered_focus_promoter = df_peakov_filtered_focus_promoter
		self.df_peakov_filtered_focus_nonpromoter = df_peakov_filtered_focus_nonpromoter

		self.df_peakov_filtered_promoter = df_peakov_filtered_promoter
		self.df_peakov_filtered_nonpromoter = df_peakov_filtered_nonpromoter

	
	def interpret_biterm(self):
		"""interpret learned biterm topic model to co-occupancy signatures"""

		info("2-{0}-4. Interpret biterm topic model ... " . format(self.count))
		promoter_topic_number = self.determine_optimal_topic_number("promoter")
		nonpromoter_topic_number = self.determine_optimal_topic_number("nonpromoter")

		# raw and filtered topic contents
		df_promoter_topics, df_filtered_promoter_topics = self.interpret_biterm_topic_contents("promoter", promoter_topic_number)
		df_nonpromoter_topics, df_filtered_nonpromoter_topics = self.interpret_biterm_topic_contents("nonpromoter", nonpromoter_topic_number)

		# plot raw topic contents
		self.plot_biterm_topic_contents(df_promoter_topics, "promoter")
		self.plot_biterm_topic_contents(df_nonpromoter_topics, "nonpromoter")

		return(df_filtered_promoter_topics, df_filtered_nonpromoter_topics)


	def interpret_biterm_topic_contents(self, region, topic_number):
		"""interpret biterm topic model to topics containing topN words"""

		# all vocab in corpus
		df_vocab = pd.read_csv("{0}/{1}_{2}_voca.txt".format(self.outpath, self.focus_TR, region), header = None, sep = "\t")

		# word distribution of each topic
		df_pw_z = pd.read_csv("{0}/{1}/k{2}.pw_z".format(self.outpath, region, topic_number), header = None, sep = " ")
		df_pw_z = df_pw_z.iloc[:,0:df_pw_z.shape[1]-1] # delete the NaN column created by reading csv
		df_pw_z.columns = df_vocab.iloc[:,1].values

		# focus context and whole context
		if region == "promoter":
			df_peakov_filtered_focus = self.df_peakov_filtered_focus_promoter
			df_peakov_filtered_whole = self.df_peakov_filtered_promoter
		else:
			df_peakov_filtered_focus = self.df_peakov_filtered_focus_nonpromoter
			df_peakov_filtered_whole = self.df_peakov_filtered_nonpromoter

		topics = []
		filtered_topics = []
		for index,row in df_pw_z.iterrows():
			
			# get the raw content of each topic
			topic_name = "topic{0}".format(index + 1)
			words_topN_index = np.argsort(row)[::-1][:self.topNword_show] # descending order of word probability and select top N word to show in the plot
			words_topN = row.index[words_topN_index] # word
			probs_topN = np.round(row[words_topN_index],3) # word probability
			zscores = stats.zscore(row) # word probability (z-score)
			zsocres_indexed = pd.Series(zscores, index = row.index)
			zscores_topN = np.round(zscores[words_topN_index],3)

			array = pd.Series([
				topic_name,
				",".join(words_topN), 
				",".join(map(str, probs_topN)),
				",".join(map(str, zscores_topN))
				])
			array.index = [
							"topic_name", 
							"word_top{0}" . format(self.topNword_show), 
							"word_prob_top{0}" . format(self.topNword_show), 
							"word_zscore_top{0}" . format(self.topNword_show)
							]

			topics.append(array)

			
			# get component words of topic (z-score >= z-score threshold)
			topic_pos_word = words_topN[zscores_topN >= self.args.zscore].values[:10] # select top-ranked components with z-score higher than z-score threshold but no more than 10 components were selected
			
			## topic positive sites
			topic_pos_peaks = df_peakov_filtered_whole.loc[df_peakov_filtered_whole[topic_pos_word].sum(axis = 1) >= (0.8 * len(topic_pos_word)), :].index.values

			## write topic positive sites to bed file
			if len(topic_pos_peaks) >= 500 and len(topic_pos_word) >= 3:

				df_peakov_raw_topic_pos = self.df_peakov_raw.loc[topic_pos_peaks, ["chrom", "start", "end"]]
				df_peakov_raw_topic_pos.loc[:, "name"] = topic_pos_peaks
				df_peakov_raw_topic_pos_sorted = df_peakov_raw_topic_pos.sort_values(by = ["chrom", "start"])
				df_peakov_raw_topic_pos_sorted.to_csv("{0}/{1}/{2}_{1}_{3}_pos_sites.bed".format(self.outpath, region, self.focus_TR, topic_name), header = False, sep = "\t", index = False)

				filtered_topics.append(["{0}_{1}" . format(self.focus_TR, topic_name), ",".join(topic_pos_word)])


		df_topics = pd.DataFrame(topics)
		df_topics.to_csv("{0}/{1}_{2}_raw_topics_contents.txt".format(self.outpath, self.focus_TR, region), header = True, index = False, sep = "\t")

		df_filtered_topics = pd.DataFrame(filtered_topics)
		if not df_filtered_topics.empty:
			df_filtered_topics.columns = ["topic_name", "component_word"]
			df_filtered_topics.to_csv("{0}/{1}_{2}_valid_topics_contents.txt".format(self.outpath, self.focus_TR, region), header = True, index = False, sep = "\t")

		return(df_topics, df_filtered_topics)

	def topic_coherence_cal(self, topic_pos_word):
		"""calculate the topic coherence using given topic component words on whole context"""
		
		PMI_topic = 0
		N = len(topic_pos_word)
		for word_j, word_k in combinations(topic_pos_word, 2):
			PMI = 0
			p_jk = float(sum((self.df_peakov_filtered[word_j]>0) & (self.df_peakov_filtered[word_k]>0))) / float(self.df_peakov_filtered.shape[0])
			p_j = float(sum((self.df_peakov_filtered[word_j]>0))) / float(self.df_peakov_filtered.shape[0])
			p_k = float(sum((self.df_peakov_filtered[word_k]>0))) / float(self.df_peakov_filtered.shape[0])
			PMI = np.log(p_jk / (p_j * p_k))
			PMI_topic += PMI
		coherence = 2 * PMI_topic / (N * (N - 1))
		return(coherence)


	def plot_biterm_topic_contents(self, df_topics, region):

		logging.getLogger().setLevel(logging.ERROR)
		fig, axs = plt.subplots(figsize = (4 * df_topics.shape[0], 4 ), nrows = 1, ncols = df_topics.shape[0])
		mycmap = cm.get_cmap("viridis")
		
		for i in range(df_topics.shape[0]):
			topic_name = df_topics.loc[i, "topic_name"]
			TRs = list(map(lambda x:x.split("_")[-1], df_topics.loc[i, "word_top{0}".format(self.topNword_show)].split(",")))
			probs = list(map(float, df_topics.loc[i, "word_zscore_top{0}".format(self.topNword_show)].split(",")))
			probs_normalize = (probs - np.min(probs)) / (np.max(probs) - np.min(probs))
			my_color = mycmap(probs_normalize)
			ax = axs[i]
			_bar = ax.barh(TRs[::-1], probs[::-1], color = my_color[::-1])
			ax.axvline(x = self.args.zscore, linestyle = "--", color = "grey", linewidth = 1.5)
			ax.set_xlabel("Probability of CAPs (z-score)")
			ax.set_title("{0}-{1}" . format(region, topic_name))

		plt.subplots_adjust(wspace = 0.4, left = 0.15, right = 0.85, bottom = 0.15)
		plt.savefig("{0}/{1}_{2}_raw_topics_contents.pdf".format(self.outpath, self.focus_TR, region))
		plt.close()
		logging.getLogger().setLevel(logging.INFO)


	def determine_optimal_topic_number(self, region):
		"""determine the optimal topics number for topics learned from focus context"""
		list_topic_cs = []

		# promoter
		for topic_number in range(int(self.args.min_signatures), int(self.args.max_signatures), 1):
			pz_d_file = "{0}/{1}/k{2}.pz_d".format(self.outpath, region, topic_number) # topic proportions for documents
			df_pz_d = pd.read_csv(pz_d_file, header = None, sep = " ")
			df_pz_d = df_pz_d.iloc[:,0:df_pz_d.shape[1]-1] # delete the NaN column created by reading csv
			variances = df_pz_d.var(axis = 0).values
			means = df_pz_d.mean(axis = 0).values
			SSN = np.log( 
						sum([variances[i] / np.square(means[i]) for i in range(len(variances))]) / len(variances)
							) # specificity score
			PSN = np.log(df_pz_d.var(axis = 1).sum() / df_pz_d.shape[0]) # purity score
			alpha = PSN / (SSN + PSN) # a parameter to balance specificity score and purity score
			CS = alpha * SSN + (1 - alpha) * PSN # combination score
			list_topic_cs.append([topic_number, SSN, PSN, CS])

		# select topic number with the maximum combination score
		df_list_topic_cs = pd.DataFrame(list_topic_cs)
		optimal_topic_numer = int(df_list_topic_cs.iloc[df_list_topic_cs[3].argmax(), 0])

		return(optimal_topic_numer)


class GF1Calculation():
	"""Calculate generalized F1 score for CAP pairs"""
	
	def __init__(self, df_peakov_filtered, args):
		self.df_peakov_filtered = df_peakov_filtered
		self.n_cpus = args.threads

	def gF1_cal(self, word_j, word_k):
		"""calculate generalized F1 score for single CAP pairs"""
		TP = sum((self.df_peakov_filtered[word_j]>0) & (self.df_peakov_filtered[word_k]>0))
		FP = sum((self.df_peakov_filtered[word_j]==0) & (self.df_peakov_filtered[word_k]>0))
		FN = sum((self.df_peakov_filtered[word_j]>0) & (self.df_peakov_filtered[word_k]==0))
		precision = TP / (TP + FP)
		recall = TP / (TP + FN)
		if TP > 0:
			gF1 = 2 * precision * recall / (precision + recall)
		else:
			gF1 = 0
		return(gF1)

	def cal_focus(self, focus_CAP):
		"""calculate generalized F1 score for pairs between focus CAP and all CAPs. 
		TR name should be 'label' but not 'factor'."""
		
		gF1_list = []
		gF1s = [self.gF1_cal(focus_CAP, CAP) for CAP in self.df_peakov_filtered.columns.values]
		gF1_list.append(gF1s)
		df_gF1s = pd.DataFrame(gF1_list)
		df_gF1s.index = [focus_CAP]
		df_gF1s.columns = self.df_peakov_filtered.columns.values
		return(df_gF1s)

class OccupancyMatrix():
	"""Generate occupancy matrix for all CAPs at genome-wide bins"""
	
	def __init__(self, genome_bins_file, df_dataset, output_prefix, args):
		
		self.genome_bins_file = genome_bins_file
		self.df_dataset = df_dataset
		self.all_factor_label = df_dataset.loc[:, "label"].values
		self.output_prefix = output_prefix
		self.args = args
		self.promoter_annotation_file = self.args.genome_annotation.promoter_file

	def intervals_intersect_split_bins(self, split_bins_file_index, split_bins_file):
		"""generate occupancy matrix (0/1, 0 denotes non-occupied and 1 denotes occupied) for split bins"""
		
		cmd = ""
		cmd += "bedtools intersect -c -e -f 0.5 -F 0.5 -a {0} -b {1} | " . format(split_bins_file, self.df_dataset.loc[self.df_dataset["label"] == self.all_factor_label[0], "file"].values[0])
		for i in range(1, len(self.all_factor_label) - 1):
			cmd += "bedtools intersect -c -e -f 0.5 -F 0.5 -a - -b {0} | " . format(self.df_dataset.loc[self.df_dataset["label"] == self.all_factor_label[i], "file"].values[0])
		cmd += "bedtools intersect -c -e -f 0.5 -F 0.5 -a - -b {0} | " . format(self.df_dataset.loc[self.df_dataset["label"] == self.all_factor_label[-1], "file"].values[0])
		cmd += """awk '{OFS=FS="\\t"}{for(i=5;i<=NF;i++){$i=($i>1?1:$i)};print $0}' - > %s_chunk%s_raw.bed && """ %(self.output_prefix, split_bins_file_index) # convert overlap count > 1 to 1
		cmd += """awk '{OFS=FS="\\t"}{sum=0;for(i=5; i<=NF; i++){sum+=$i};if(sum>=1)print $0}' %s_chunk%s_raw.bed > %s_chunk%s.bed && rm %s_chunk%s_raw.bed & \n""" %(self.output_prefix, split_bins_file_index, self.output_prefix, split_bins_file_index, self.output_prefix, split_bins_file_index) # remove bins without any CAP occupancy
		# info(cmd)
		
		return(cmd)

	def intervals_intersect_genome_bins(self):
		"""generate occupancy matrix for all genome-wide bins (in parallel)"""
		
		header = "chrom\tstart\tend\tname\t{0}\n" . format("\t".join(self.all_factor_label))
		outf_peakov = open("{0}_peakov.bed".format(self.output_prefix), "w")
		outf_peakov.write(header)
		outf_peakov.close()

		# split genome-wide bins to sub files for calculating occupancy matrix in parallel (not work in macOS)
		cmd_split_bins = "split --additional-suffix _{0} -n l/{1} -d {2}" . format(self.output_prefix, self.args.threads, self.genome_bins_file, self.args.out_dir) # split genome-wide bins to multiple files
		subprocess.run(cmd_split_bins, shell = True, check = True)
		split_bins_files = [i for i in os.listdir(os.getcwd()) if i.endswith(self.output_prefix)]
		outf = open("{0}_peakov.sh".format(self.output_prefix), "w")
		for i in range(len(split_bins_files)):
			cmd_i = self.intervals_intersect_split_bins(i, split_bins_files[i])
			outf.write(cmd_i)
		outf.write("wait\n")
		outf.write("cat {0}_chunk*.bed >> {0}_peakov.bed\nwait".format(self.output_prefix))
		outf.close()

		# execute shell script
		subprocess.run("bash {0}_peakov.sh && rm {0}_chunk*.bed && rm x*_{0}".format(self.output_prefix), shell=True, check=True)

		# split promoter and non-promoter regions
		cmd_split_promoter = "bedtools intersect -wa -u -e -f 0.5 -F 0.5 -a {0}_peakov.bed -b {1} > {0}_peakov_promoter.bed" . format(self.output_prefix, self.promoter_annotation_file)
		cmd_split_nonpromoter = "bedtools intersect -wa -v -e -f 0.5 -F 0.5 -a {0}_peakov.bed -b {1} > {0}_peakov_nonpromoter.bed" . format(self.output_prefix, self.promoter_annotation_file)
		subprocess.run(cmd_split_promoter, shell = True, check = True)
		subprocess.run(cmd_split_nonpromoter, shell = True, check = True)

		return("{0}_peakov.bed".format(self.output_prefix)) # file of genome-wide occupancy matrix

class BWMatrix():
	"""scan signals from big wiggle track file"""
	
	def __init__(self, bw_file, point_number = 1, process_number = 10, anno = True):
		self.bw_file = bw_file
		self.point_number = point_number
		self.process_number = min(process_number, multiprocessing.cpu_count() - 1)
		self.anno = anno

	def stat_unit(self, array, bw):
		try:
			score = pd.Series(bw.stats(array[0], array[1], array[2], nBins = self.point_number)).astype("float")
		except RuntimeError:
			# error("Invalid interval bounds, the interval is {0} for {1}" . format(list(array), self.bw_file))
			score = pd.Series([np.nan] * self.point_number)
		return(score)

	def stat_matrix(self, df_split):
		bw = pyBigWig.open(self.bw_file, "r")
		if not bw.isBigWig():
			error("Input {0} is not a big wiggle file, exit!" . format(self.bw_file))
			sys.exit(1)
		scores = df_split.apply(self.stat_unit, args = (bw, ), axis = 1, raw = False)
		bw.close()
		if self.anno:
			return(pd.concat([df_split, scores], axis = 1))
		else:
			return(scores)

	def stat_matrix_parallel(self, df):
		"""scan signals on the track provided by pandas dataframe format"""
		df_splits = np.array_split(df, self.process_number)
		pool = Pool(self.process_number)
		df_out = pd.concat(pool.map(self.stat_matrix, df_splits))
		pool.close()
		pool.join()
		return(df_out)

class Motif_OverlapMatrix():
	"""Motif processing to generate overlap matrix at genome-wide bins"""
	
	def __init__(self, df_motif_factor_uid_dataset, args):
		
		self.args = args
		self.output_prefix = "{0}_1kb_bins" . format(self.args.name)
		self.genomewide_bins_file = "{0}/{1}_1kb_bins_peakov.bed" . format(self.args.Sig_InputPath, self.args.Sig_InputName)
		self.genomewide_bins_loci_file = "{0}_1kb_bins_peakov_noheader.bed" . format(self.args.Sig_InputName)
		
		self.df_motif_file = df_motif_factor_uid_dataset
		self.all_motif_label = self.df_motif_file.loc[:, "label"].values
		
		self.n_cpus = self.args.threads

	def intervals_intersect_split_bins(self, split_bins_file_index, split_bins_file):
		"""generate motif overlap matrix for split bins"""
		cmd = ""
		cmd += "bedtools intersect -c -e -f 0.5 -F 0.5 -a {0} -b {1} | " . format(split_bins_file, self.df_motif_file.loc[self.df_motif_file["label"] == self.all_motif_label[0], "file"].values[0])
		for i in range(1, len(self.all_motif_label) - 1):
			cmd += "bedtools intersect -c -e -f 0.5 -F 0.5 -a - -b {0} | " . format(self.df_motif_file.loc[self.df_motif_file["label"] == self.all_motif_label[i], "file"].values[0])
		cmd += "bedtools intersect -c -e -f 0.5 -F 0.5 -a - -b {0} > {1}_chunk{2}.bed & \n" . format(self.df_motif_file.loc[self.df_motif_file["label"] == self.all_motif_label[len(self.all_motif_label) - 1], "file"].values[0], self.output_prefix, split_bins_file_index)
		return(cmd)

	def intervals_intersect_genomewide_bins(self):
		"""generate overlap matrix for all genome-wide bins (in parallel)"""
		header = "chrom\tstart\tend\tname\t{0}\n" . format("\t".join(self.all_motif_label))
		outf_peakov = open("{0}_motifov.bed".format(self.output_prefix), "w")
		outf_peakov.write(header)
		outf_peakov.close()

		cmd_noheader = "cut -f 1-4 {0} > {1}" . format(self.genomewide_bins_file, self.genomewide_bins_loci_file)
		subprocess.run(cmd_noheader, shell = True, check = True)
		
		cmd_split_bins = "split --additional-suffix _{0} -n l/{1} -d {2}" . format(self.output_prefix, self.n_cpus, self.genomewide_bins_loci_file) # split genome-wide bins to multiple files
		subprocess.run(cmd_split_bins, shell = True, check = True)
		split_bins_files = [i for i in os.listdir(os.getcwd()) if i.endswith(self.output_prefix)]
		outf = open("{0}_motifov.sh".format(self.output_prefix), "w")
		for i in range(len(split_bins_files)):
			cmd_i = self.intervals_intersect_split_bins(i, split_bins_files[i])
			outf.write(cmd_i)
		outf.write("wait\n")
		outf.write("cat {0}_chunk*.bed >> {0}_motifov.bed\nwait".format(self.output_prefix))
		outf.close()

		subprocess.run("bash {0}_motifov.sh && rm *_{0} && rm {0}_chunk*.bed".format(self.output_prefix), shell = True, check = True)
		# subprocess.run("bash {0}_motifov.sh".format(self.output_prefix), shell = True, check = True)
		return("{0}_motifov.bed" . format(self.output_prefix))



