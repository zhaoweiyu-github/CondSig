#!/usr/bin/env python3


# ------------------------------
#  python modules
# ------------------------------

import os,sys
import pandas as pd
import numpy as np
import subprocess
from scipy import stats
from functools import partial
import multiprocessing
from itertools import combinations

import shap
from sklearn import metrics
from xgboost.sklearn import XGBClassifier

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
from CondSig import Utility

# ------------------------------
#  Misc functions
# ------------------------------

class SigFilter():
	"""Filter condensate-like signatures and remove redudant signatures"""

	def __init__(self, args):
		self.args = args
		self.individual_auc_threshold = 0.6
		self.mean_auc_threshold = 0.65
		self.mean_shap_threshold = 1
		self.JI_TR_top5_threshold = 1.0 / 4

	def filter(self, df_filtered_topics, region):

		CLsig_features_list= []
		# filter condensate-like features
		for index, row in df_filtered_topics.iterrows():
			SigName = row["topic_name"]
			focus_CAP = "_".join(SigName.split("_")[:-1])
			sig_idx = SigName.split("_")[-1]
			
			ROC_file = "ROC/{0}_{1}_{2}_AUROC.txt" . format(
				focus_CAP,
				region,
				sig_idx
			)

			df_ROC = pd.read_csv(ROC_file, header = 0, sep = "\t", index_col = "feature")
			# df_XGB = pd.read_csv(XGB_file, header = 0, sep = "\t", index_col = "feature")

			condensate_features = ["PPI", "RBS", "LLPS", "MLO", "IDR", "RBP"]
			condensate_features = np.intersect1d(df_ROC.index.values, condensate_features) # to avoid unavailable RNA binding strength
			
			roc_status = df_ROC.loc[condensate_features, "AUROC"].astype("float64")
			# xgb_status = df_XGB.loc[condensate_features, "mean_shap"].astype("float64")
			
			# condensate-like feature positive corresponds to signatures
			condensate_features_pos = roc_status[roc_status >= self.individual_auc_threshold].index.values
			if len(condensate_features_pos) >= 3:
				# top 3 condensate-like features positively respond to signature 
				condensate_features_pos_top3 = roc_status[condensate_features_pos].nlargest(3).index.values
				if roc_status[condensate_features_pos_top3].mean() >= self.mean_auc_threshold:
					CLsig_features_list.append([SigName, ",".join(condensate_features_pos), len(condensate_features_pos), roc_status[condensate_features_pos].mean()])

		df_CLsig_features = pd.DataFrame(CLsig_features_list)
		# df_CLsig_features.columns = ["signature_name", "qualified_CL_features", "qualified_CL_features_count", "mean_AUROC", "mean_feature_importance"]
		df_CLsig_features.columns = ["signature_name", "qualified_CL_features", "qualified_CL_features_count", "mean_AUROC"]
		df_CLsig_features.index = df_CLsig_features.loc[:, "signature_name"].values
		df_CLsig_features.to_csv("CondSig/{0}_{1}_nonunique_condensate_like_signatures.txt" . format(self.args.name, region), header = True, index = False, sep = "\t")
		CLsigs = df_CLsig_features.loc[:, "signature_name"].values

		# remove redudant condensate-like signatures
		df_filtered_topics.index = df_filtered_topics.loc[:, "topic_name"].values
		df_filtered_topics_CLsigs = df_filtered_topics.loc[CLsigs, :]
		# sort signatures based on count and level of enriched condensation features
		df_CLsig_features_sorted = df_CLsig_features.sort_values(by = ["qualified_CL_features_count", "mean_AUROC"], ascending = False)

		list_JI_sigs = []
		for sig_i, sig_j in combinations(df_CLsig_features_sorted.index.values, 2):
			TR_sig_i = df_filtered_topics_CLsigs.loc[sig_i, "component_word"].split(",")
			TR_sig_j = df_filtered_topics_CLsigs.loc[sig_j, "component_word"].split(",")

			# calculate JI between top 5 sig-i-pos TRs and sig-j-pos TRs
			JI_TR_top5 = len(np.intersect1d(TR_sig_i[:5], TR_sig_j[:5])) / len(np.union1d(TR_sig_i[:5], TR_sig_j[:5])) # v.1.2.2
			
			# calculate JI between sig-i-pos and sig-j-pos sites
			focus_TR_i = "_".join(sig_i.split("_")[:-1])
			sig_idx_i = sig_i.split("_")[-1]
			focus_TR_j = "_".join(sig_j.split("_")[:-1])
			sig_idx_j = sig_j.split("_")[-1]
			
			list_JI_sigs.append([sig_i, sig_j, JI_TR_top5])

		df_JI_sigs = pd.DataFrame(list_JI_sigs)
		df_JI_sigs.columns = ["primary", "secondary", "JI_TR_top5"]

		df_JI_sigs_grouped = df_JI_sigs.groupby("primary")
		redundant_sigs = []
		for group, df_grpup in df_JI_sigs_grouped:
			if group not in redundant_sigs:
				redundant_sigs += list(df_grpup.loc[df_grpup["JI_TR_top5"] > self.JI_TR_top5_threshold, "secondary"].values)
		df_unique_CLsig = df_filtered_topics_CLsigs.loc[np.setdiff1d(df_filtered_topics_CLsigs.index.values, redundant_sigs), :]
		df_unique_CLsig_featues = df_CLsig_features.loc[df_unique_CLsig.index.values, ["qualified_CL_features", "qualified_CL_features_count", "mean_AUROC"]]

		df_unique_CLsig_anno = pd.concat([df_unique_CLsig, df_unique_CLsig_featues], axis = 1)
		df_unique_CLsig_anno = df_unique_CLsig_anno.sort_values(by = ["qualified_CL_features_count", "mean_AUROC"], ascending = False)
		df_unique_CLsig_anno.columns = ['signature_name', 'component_CAP', 'qualified_CL_features', 'qualified_CL_features_count', 'mean_AUROC']
		
		df_unique_CLsig_anno.to_csv("CondSig/{0}_{1}_unique_condensate_like_signatures.txt" . format(self.args.name, region), header = True, index = False, sep = "\t")

		# summary fraction
		df_CondSig = df_unique_CLsig_anno.copy()
		# move signature name to CondSig name
		df_CondSig.loc[:, 'signature_name'] = ["{0}_{1}_CondSig_{2}" . format(self.args.name, region, i+1) for i in range(df_CondSig.shape[0])]
		df_CondSig.columns = ['CondSig', 'component_CAP', 'qualified_CL_features', 'qualified_CL_features_count', 'mean_AUROC']
		df_CondSig.to_csv("Summary/{0}_{1}_CondSigs.txt" . format(self.args.name, region), header = True, index = False, sep = "\t")

		os.makedirs('Summary/CondSig_pos_sites', exist_ok=True)
		for i in range(df_CondSig.shape[0]):
			SigName = df_unique_CLsig_anno.loc[:,"signature_name"].values[i]
			focus_CAP = "_".join(SigName.split("_")[:-1])
			sig_idx = SigName.split("_")[-1]
			peaks_sig_pos = "Sites/{0}_{1}_{2}_pos_sites.bed" . format(
				focus_CAP,
				region,
				sig_idx)
			CondSigName=df_CondSig.iloc[i,:]['CondSig']
			peaks_Condsig_pos = "Summary/CondSig_pos_sites/{0}_pos_sites.bed" . format(CondSigName)
			cmd_copy = 'cp {0} {1}' . format(peaks_sig_pos, peaks_Condsig_pos)
			subprocess.run(cmd_copy, shell = True, check = True)
			
class SigEval():
	"""Evaluation of condensate-like capacity for each signature was performed by checking 
	whether signature-positive sites showed higher condensate-like features than signature-negative sites"""

	def __init__(self, args):

		self.args = args
		self.sampling_times = 100

	def run(self, SigName, region):

		focus_CAP = "_".join(SigName.split("_")[:-1])
		sig_idx = SigName.split("_")[-1]

		df_sig_pos_FE = pd.read_csv("Features/{0}_{1}_{2}_pos_FE.txt" . format(
				focus_CAP,
				region,
				sig_idx), header = 0, sep = "\t")

		df_sig_neg_FE = pd.read_csv("Features/{0}_{1}_{2}_neg_FE.txt" . format(
				focus_CAP,
				region,
				sig_idx), header = 0, sep = "\t")

		df_sig_pos_FE = df_sig_pos_FE.loc[:, np.setdiff1d(df_sig_pos_FE.columns, ["chrom", "start", "end"])]
		df_sig_neg_FE = df_sig_neg_FE.loc[:, np.setdiff1d(df_sig_neg_FE.columns, ["chrom", "start", "end"])]

		df_sig_pos_FE.loc[:, "label"] = True
		df_sig_neg_FE.loc[:, "label"] = False

		func_roc_analysis_balanced = partial(roc_analysis_balanced, df_sig_pos_FE, df_sig_neg_FE)
		pool = multiprocessing.Pool(self.args.threads)
		df_auc_list = pool.map(func_roc_analysis_balanced, np.arange(self.sampling_times))
		pool.close()
		pool.join()

		df_auc_merged = pd.concat(df_auc_list, axis = 1)
		df_auc_mean = pd.DataFrame([df_auc_merged.loc[:, "feature"].iloc[:, 0].values,
											df_auc_merged.loc[:, "AUROC"].mean(axis = 1).values
											]).T
		df_auc_mean.columns = ["feature", "AUROC"]
		df_auc_mean.loc[:, "AUROC"] = df_auc_mean.loc[:, "AUROC"].astype("float")
		df_auc_mean.loc[:, "AUROC"] = df_auc_mean.loc[:, "AUROC"].round(decimals = 3)

		self.roc_output(df_auc_mean, focus_CAP, region, sig_idx)


	def roc_output(self, df_auc, focus_CAP, region, sig_idx):
		"""outplot ROC analysis result"""

		df_auc = df_auc.sort_values(by = "AUROC", ascending = True)
		df_auc.to_csv("ROC/{0}_{1}_{2}_AUROC.txt" . format(
				focus_CAP,
				region,
				sig_idx
			), header = True, sep = "\t", index = False)

		# plot bar plot for auc
		logging.getLogger().setLevel(logging.ERROR)
		fig, ax1 = plt.subplots(figsize = (4.5,4), nrows = 1, ncols = 1)
		mycmap = cm.get_cmap("viridis")
		_bar = ax1.barh(df_auc.loc[:, "feature"], df_auc.loc[:, "AUROC"], color = mycmap(df_auc.loc[:, "AUROC"]), height = 0.7)
		ax1.axvline(x = 0.6, linestyle = "--", color = "grey", linewidth = 1.5)
		ax1.set_xlabel("AUROC")
		ax1.set_title("AUROC")
		ax1.set_xlim(0, 1)

		plt.subplots_adjust(left = 0.45, right = 0.95, bottom = 0.15)
		plt.savefig("ROC/{0}_{1}_{2}_AUROC.pdf" . format(
				focus_CAP,
				region,
				sig_idx))
		plt.close()
		logging.getLogger().setLevel(logging.INFO)


class SigFE():
	"""Feature enginerring for learned co-occupancy signatures"""

	def __init__(self, args):

		self.args = args
		
		self.df_peakov = pd.read_csv("{0}/{1}_1kb_bins_peakov.bed" . format(self.args.Sig_InputPath, self.args.Sig_InputName), header = 0, sep = "\t", index_col = "name")
		self.df_LLPS = pd.read_csv(self.args.LLPS, header = 0, sep = "\t").dropna()
		self.df_MLO = pd.read_csv(self.args.MLO, header = 0, sep = "\t").dropna()
		self.df_PPI_merged = pd.read_csv(self.args.PPI, header = 0, sep = "\t")

		# self.tss_file = self.args.genome_annotation.TSS_file
		
		self.df_dataset = pd.read_csv(self.args.data_annotation, header = None, sep = "\t").iloc[:,0:4]
		self.df_dataset.columns = ["factor", "label", "file", "uniprot_id"]
		self.df_dataset.index = self.df_dataset.loc[:, "label"].values
		uni_uniprot_id, count_uniprot_id = np.unique(self.df_dataset.loc[:, "uniprot_id"].values, return_counts = True)
		if sum(count_uniprot_id > 1) > 0:
			error("More than 1 CAPs are assigned to the same uniprot id ({0}), please check it and don't use redudant data for each CAP." . format(uni_uniprot_id[count_uniprot_id > 1]))
			sys.exit(1)

		df_peakov_filtered = pd.read_csv("{0}/{1}_1kb_bins_peakov_filtered.bed" . format(self.args.Sig_InputPath, self.args.Sig_InputName), header = 0, sep = "\t", index_col = "name")
		promoter_bins = pd.read_csv("{0}/{1}_1kb_bins_peakov_promoter.bed" . format(self.args.Sig_InputPath, self.args.Sig_InputName), header = None, sep = "\t").iloc[:,3].values
		self.df_peakov_filtered_promoter = df_peakov_filtered.loc[df_peakov_filtered.index.isin(promoter_bins),:]
		self.df_peakov_filtered_nonpromoter = df_peakov_filtered.loc[~df_peakov_filtered.index.isin(promoter_bins),:]		
		

	def FE_preprocess(self):

		"""preprocessing input feature file"""

		# protein-protein interaction between proteins of dataset (to reduce the execution time of FE run)
		df_PPI_dataset = self.df_PPI_merged.loc[
						((
							self.df_PPI_merged["SWISS-PROT Accessions Interactor A"].isin(self.df_dataset.loc[:, "uniprot_id"].values)
							) & (
								self.df_PPI_merged["SWISS-PROT Accessions Interactor B"].isin(self.df_dataset.loc[:, "uniprot_id"].values)
								)), :]
		self.df_PPI_dataset = df_PPI_dataset

		# LLPS capacity
		df_LLPS_status = self.df_dataset.loc[:, ["label", "uniprot_id"]]
		df_LLPS_status.loc[:, "status"] = 0
		df_LLPS_status.loc[df_LLPS_status["uniprot_id"].isin(self.df_LLPS.loc[:, "UniProt ID"].dropna().values), "status"] = 1
		df_LLPS_status.index = df_LLPS_status.loc[:, "uniprot_id"].values
		self.df_LLPS_status = df_LLPS_status

		# disordered domain content
		disorder_content_cutoff = 0.153 # determined by IDR percentage of LLPS proteins in human
		df_mobidb_lite = pd.read_csv(self.args.IDR, header = None, sep = "\t")
		df_mobidb_lite.columns = ["uniprot_id", "disorder_content"]
		df_mobidb_lite.loc[:, "uniprot_id"] = df_mobidb_lite.loc[:, "uniprot_id"].str.upper()
		df_mobidb_lite.loc[:, "status"] = 0
		df_mobidb_lite.loc[df_mobidb_lite["disorder_content"] >= disorder_content_cutoff, "status"] = 1
		df_mobidb_lite.index = df_mobidb_lite.loc[:, "uniprot_id"].values
		self.df_IDR_status = df_mobidb_lite

		# RNA binding domain content
		df_RBD = pd.read_csv(self.args.RBD, header = 0, sep = "\t")
		RBD_uids = df_RBD.loc[:, "uid"].dropna().values
		# self.RBD_uids = RBD_uids
		df_RBD_status = self.df_dataset.loc[:, ["label", "uniprot_id"]]
		df_RBD_status.loc[:, "status"] = 0
		df_RBD_status.loc[df_RBD_status["uniprot_id"].isin(RBD_uids), "status"] = 1
		df_RBD_status.index = df_RBD_status.loc[:, "uniprot_id"].values
		self.df_RBD_status = df_RBD_status

	def FE_run(self, SigName, component_CAPs, region):

		focus_CAP = "_".join(SigName.split("_")[:-1])
		sig_idx = SigName.split("_")[-1]

		potential_combinatorial_TRs = pd.read_csv(
			"{0}/LearnSig/{1}/{1}_{2}_voca.txt" . format(
				self.args.Sig_InputPath,
				focus_CAP,
				region
				),
			header = None, sep = "\t"
			).iloc[:, 1].values

		if region == "promoter":
			df_peakov_sites_region = self.df_peakov_filtered_promoter.loc[:, potential_combinatorial_TRs]
		else:
			df_peakov_sites_region = self.df_peakov_filtered_nonpromoter.loc[:, potential_combinatorial_TRs]

		pos_sites = df_peakov_sites_region.loc[
			df_peakov_sites_region.loc[:, component_CAPs].sum(axis = 1) >= (0.8 * len(component_CAPs)), :].index.values
		
		neg_sites = df_peakov_sites_region.loc[
		(df_peakov_sites_region.sum(axis = 1) >= (0.8 * len(component_CAPs))) & (df_peakov_sites_region.loc[:, component_CAPs].sum(axis = 1) < 2), :].index.values

		df_peaks_sig_pos = self.df_peakov.loc[pos_sites, ["chrom", "start",	"end"]].copy()
		df_peaks_sig_neg = self.df_peakov.loc[neg_sites, ["chrom", "start",	"end"]].copy()

		df_peaks_sig_pos.to_csv("Sites/{0}_{1}_{2}_pos_sites.bed" . format(
				focus_CAP,
				region,
				sig_idx), header = False, sep = "\t", index = False)

		df_peakov_sig_pos = self.df_peakov.loc[pos_sites, potential_combinatorial_TRs].copy()
		df_peakov_sig_neg = self.df_peakov.loc[neg_sites, potential_combinatorial_TRs].copy()
		# info(df_peakov_sig_pos.shape)
		# info(df_peakov_sig_neg.shape)
		# info(np.mean(df_peakov_sig_pos.sum(axis = 1)))
		# info(np.mean(df_peakov_sig_neg.sum(axis = 1)))

		# return(None)

		# convert factor label to uniprot ids
		df_peakov_sig_pos_uids = df_peakov_sig_pos.copy()
		df_peakov_sig_neg_uids = df_peakov_sig_neg.copy()
		df_peakov_sig_pos_uids.columns = self.df_dataset.loc[potential_combinatorial_TRs, "uniprot_id"].values
		df_peakov_sig_neg_uids.columns = self.df_dataset.loc[potential_combinatorial_TRs, "uniprot_id"].values

		df_peakov_sig_pos_uids_splits = np.array_split(df_peakov_sig_pos_uids, self.args.threads)
		df_peakov_sig_neg_uids_splits = np.array_split(df_peakov_sig_neg_uids, self.args.threads)

		# Protein-protein interaction frequency
		info("protein-protein interaction")
		func_PPI = partial(PPI_frequency_cal, self.df_PPI_dataset)
		pool = multiprocessing.Pool(self.args.threads)
		PPI_frequency_pos = np.concatenate(pool.map(func_PPI, df_peakov_sig_pos_uids_splits))
		pool.close()
		pool.join()

		pool = multiprocessing.Pool(self.args.threads)
		PPI_frequency_neg = np.concatenate(pool.map(func_PPI, df_peakov_sig_neg_uids_splits))
		pool.close()		
		pool.join()

		# LLPS frequency
		info("LLPS")
		func_LLPS = partial(LLPS_frequency_cal, self.df_LLPS_status)
		pool = multiprocessing.Pool(self.args.threads)
		LLPS_frequency_pos = np.concatenate(pool.map(func_LLPS, df_peakov_sig_pos_uids_splits))
		pool.close()
		pool.join()

		pool = multiprocessing.Pool(self.args.threads)
		LLPS_frequency_neg = np.concatenate(pool.map(func_LLPS, df_peakov_sig_neg_uids_splits))
		pool.close()
		pool.join()

		# MLO frequency
		info("MLO")
		func_MLO = partial(MLO_frequency_cal, self.df_MLO)
		pool = multiprocessing.Pool(self.args.threads)
		MLO_frequency_pos = np.concatenate(pool.map(func_MLO, df_peakov_sig_pos_uids_splits))
		pool.close()
		pool.join()

		pool = multiprocessing.Pool(self.args.threads)
		MLO_frequency_neg = np.concatenate(pool.map(func_MLO, df_peakov_sig_neg_uids_splits))
		pool.close()
		pool.join()

		# disorder frequency
		info("disorder frequency")
		func_disorder = partial(disorder_frequency_cal, self.df_IDR_status)
		pool = multiprocessing.Pool(self.args.threads)
		disorder_pos = np.concatenate(pool.map(func_disorder, df_peakov_sig_pos_uids_splits))
		pool.close()
		pool.join()
		
		pool = multiprocessing.Pool(self.args.threads)
		disorder_neg = np.concatenate(pool.map(func_disorder, df_peakov_sig_neg_uids_splits))
		pool.close()
		pool.join()

		# RNA-binding domain frequency
		info("RNA-binding domain")
		func_RBD = partial(RBD_frequency_cal, self.df_RBD_status)
		pool = multiprocessing.Pool(self.args.threads)
		RBD_pos = np.concatenate(pool.map(func_RBD, df_peakov_sig_pos_uids_splits))
		pool.close()
		pool.join()

		pool = multiprocessing.Pool(self.args.threads)
		RBD_neg = np.concatenate(pool.map(func_RBD, df_peakov_sig_neg_uids_splits))
		pool.close()
		pool.join()		

		# with RNA binding strength
		if self.args.RBS != None:
			# RNA binding strength
			RNA_pos = self.RNA_FE(df_peaks_sig_pos, self.args.RBS, self.args.threads)
			RNA_neg = self.RNA_FE(df_peaks_sig_neg, self.args.RBS, self.args.threads)

			df_sig_pos_FE = pd.concat([df_peaks_sig_pos, 
				pd.DataFrame(PPI_frequency_pos, index = df_peakov_sig_pos_uids.index.values),
				pd.DataFrame(RNA_pos),
				pd.DataFrame(LLPS_frequency_pos, index = df_peakov_sig_pos_uids.index.values),
				pd.DataFrame(MLO_frequency_pos, index = df_peakov_sig_pos_uids.index.values),
				pd.DataFrame(disorder_pos, index = df_peakov_sig_pos_uids.index.values),
				pd.DataFrame(RBD_pos, index = df_peakov_sig_pos_uids.index.values)
				], axis = 1)

			df_sig_pos_FE.index = df_peaks_sig_pos.index.values			
			# df_sig_pos_FE.columns = ["chrom", "start", "end", "PPI", "RNA binding strength", "LLPS capacity", "MLO", "Disordered domain", "RNA binding domain"]
			df_sig_pos_FE.columns = ["chrom", "start", "end", "PPI", "RBS", "LLPS", "MLO", "IDR", "RBP"]

			df_sig_neg_FE = pd.concat([df_peaks_sig_neg, 
				pd.DataFrame(PPI_frequency_neg, index = df_peakov_sig_neg_uids.index.values),
				pd.DataFrame(RNA_neg),
				pd.DataFrame(LLPS_frequency_neg, index = df_peakov_sig_neg_uids.index.values),
				pd.DataFrame(MLO_frequency_neg, index = df_peakov_sig_neg_uids.index.values),
				pd.DataFrame(disorder_neg, index = df_peakov_sig_neg_uids.index.values),
				pd.DataFrame(RBD_neg, index = df_peakov_sig_neg_uids.index.values)
				], axis = 1)

			df_sig_neg_FE.index = df_peaks_sig_neg.index.values
			df_sig_neg_FE.columns = ["chrom", "start", "end", "PPI", "RBS", "LLPS", "MLO", "IDR", "RBP"]

			df_sig_pos_FE.to_csv("Features/{0}_{1}_{2}_pos_FE.txt" . format(
				focus_CAP,
				region,
				sig_idx), header = True, sep = "\t", index = False)

			df_sig_neg_FE.to_csv("Features/{0}_{1}_{2}_neg_FE.txt" . format(
				focus_CAP,
				region,
				sig_idx), header = True, sep = "\t", index = False)

		# without RNA binding strength
		else:
			df_sig_pos_FE = pd.concat([df_peaks_sig_pos, 
				# pd.DataFrame(CA_pos), 
				# pd.DataFrame(motif_frequency_pos, index = df_peakov_sig_pos_uids.index.values), 
				# pd.DataFrame(distance_pos_log10, index = df_peakov_sig_pos_uids.index.values), 
				pd.DataFrame(PPI_frequency_pos, index = df_peakov_sig_pos_uids.index.values),
				pd.DataFrame(LLPS_frequency_pos, index = df_peakov_sig_pos_uids.index.values),
				pd.DataFrame(MLO_frequency_pos, index = df_peakov_sig_pos_uids.index.values),
				pd.DataFrame(disorder_pos, index = df_peakov_sig_pos_uids.index.values),
				pd.DataFrame(RBD_pos, index = df_peakov_sig_pos_uids.index.values)
				], axis = 1)

			df_sig_pos_FE.index = df_peaks_sig_pos.index.values
			df_sig_pos_FE.columns = ["chrom", "start", "end", "PPI", "LLPS", "MLO", "IDR", "RBP"]

			df_sig_neg_FE = pd.concat([df_peaks_sig_neg, 
				# pd.DataFrame(CA_neg), 
				# pd.DataFrame(motif_frequency_neg, index = df_peakov_sig_neg_uids.index.values), 
				# pd.DataFrame(distance_neg_log10, index = df_peakov_sig_neg_uids.index.values), 
				pd.DataFrame(PPI_frequency_neg, index = df_peakov_sig_neg_uids.index.values),
				pd.DataFrame(LLPS_frequency_neg, index = df_peakov_sig_neg_uids.index.values),
				pd.DataFrame(MLO_frequency_neg, index = df_peakov_sig_neg_uids.index.values),
				pd.DataFrame(disorder_neg, index = df_peakov_sig_neg_uids.index.values),
				pd.DataFrame(RBD_neg, index = df_peakov_sig_neg_uids.index.values)
				], axis = 1)

			df_sig_neg_FE.index = df_peaks_sig_neg.index.values
			df_sig_neg_FE.columns = ["chrom", "start", "end", "PPI", "LLPS", "MLO", "IDR", "RBP"]

			df_sig_pos_FE.to_csv("Features/{0}_{1}_{2}_pos_FE.txt" . format(
				focus_CAP,
				region,
				sig_idx), header = True, sep = "\t", index = False)

			df_sig_neg_FE.to_csv("Features/{0}_{1}_{2}_neg_FE.txt" . format(
				focus_CAP,
				region,
				sig_idx), header = True, sep = "\t", index = False)


	def CA_FE(self, df_sites, ca_bw, process_number):
		"""calculate chromatin accessibility at the given site"""
		
		df_sites_window = df_sites.copy()
		bwMatrix = Utility.BWMatrix(ca_bw, point_number = 1, process_number = process_number, anno = False)
		df_sites_ca = bwMatrix.stat_matrix_parallel(df_sites_window).fillna(0)
		
		return(df_sites_ca)

	def RNA_FE(self, df_sites, RNA_bw, process_number):
		"""calculate RNA binding strength at the given site from big wiggle track"""
		
		df_sites_window = df_sites.copy()
		bwMatrix = Utility.BWMatrix(RNA_bw, point_number = 1, process_number = process_number, anno = False)
		df_sites_rna = bwMatrix.stat_matrix_parallel(df_sites_window).fillna(0)
		
		return(df_sites_rna)

def roc_analysis_balanced(df_sig_pos_FE, df_sig_neg_FE, sampling_count):
	"""perform ROC analysis for features at balanced positive sites and negative sites"""
	
	if df_sig_pos_FE.shape[0] >= df_sig_neg_FE.shape[0]:
		df_sig_merged_FE = pd.concat([df_sig_pos_FE.sample(n = df_sig_neg_FE.shape[0], replace = False), df_sig_neg_FE], axis = 0)
	else:
		df_sig_merged_FE = pd.concat([df_sig_pos_FE, df_sig_neg_FE.sample(n = df_sig_pos_FE.shape[0], replace = False)], axis = 0)

	X = df_sig_merged_FE.loc[:, np.setdiff1d(df_sig_merged_FE.columns, ["label"])]
	y = df_sig_merged_FE.loc[:, "label"]

	# calculate AUROC and plot ROC curve
	list_auc = []
	for feature in X.columns:
		x = X.loc[:, feature].values
		# fpr, tpr, thresholds = metrics.roc_curve(y, x)
		auc = round(metrics.roc_auc_score(y,x), 2)
		list_auc.append([feature, auc])
	df_auc = pd.DataFrame(list_auc)
	df_auc.columns = ["feature", "AUROC"]

	return(df_auc)


def PPI_frequency_cal(df_PPI_dataset, df_peakov_selected):
	"""calculate protein-protein interaction frequency"""

	PPI_frequencies = []
	for index, row in df_peakov_selected.iterrows():
		peakov_factor = row[row > 0].index.values
		peakov_factor_paired_count = len(peakov_factor) * (len(peakov_factor) - 1) * 0.5 # All potential pairs of N occupied factors: C(N,2)
		PPI_peakov_factors_paired_count = df_PPI_dataset.loc[(
				(
				df_PPI_dataset["SWISS-PROT Accessions Interactor A"].isin(peakov_factor)
				) & (
				df_PPI_dataset["SWISS-PROT Accessions Interactor B"].isin(peakov_factor)
				)
			), :].shape[0]
		PPI_frequency = PPI_peakov_factors_paired_count / peakov_factor_paired_count
		PPI_frequencies.append(PPI_frequency)

	return(PPI_frequencies)

def LLPS_frequency_cal(df_LLPS_status, df_peakov_selected):
	"""calculate LLPS frequency of occupied factors"""
	
	df_LLPS_frequencies = df_peakov_selected * df_LLPS_status.loc[df_peakov_selected.columns.values, "status"] # 0/1 matrix * 0/1 status array
	LLPS_frequencies = np.array(df_LLPS_frequencies.sum(axis = 1) / df_peakov_selected.sum(axis = 1))
		
	return(LLPS_frequencies)

def MLO_frequency_cal(df_MLO, df_peakov_selected):
	"""calculate MLO frequency of occupied factors"""

	MLO_protein_frequencies = []
	for index, row in df_peakov_selected.iterrows():
		peakov_uids = row[row > 0].index.values
		MLO_protein = []
		for index, row in df_MLO.iterrows():
			peakov_factor_inside_condensate = np.intersect1d(peakov_uids, row["Uniprot_IDs"].split(","))
			if(len(peakov_factor_inside_condensate)) >= 2:
				MLO_protein += peakov_factor_inside_condensate.tolist() # count factors co-occured in the same condensate with other co-occupied factor
				
		MLO_protein_frequency = len(np.unique(MLO_protein)) / len(peakov_uids)
		MLO_protein_frequencies.append(MLO_protein_frequency)

	return(MLO_protein_frequencies)

def disorder_frequency_cal(df_IDR_status, df_peakov_selected):
		
	df_disorder_frequencies = df_peakov_selected * df_IDR_status.loc[df_peakov_selected.columns.values, "status"]
	disorder_frequencies = np.array(df_disorder_frequencies.sum(axis = 1) / df_peakov_selected.sum(axis = 1))

	return(disorder_frequencies)

def RBD_frequency_cal(df_RBD_status, df_peakov_selected):
		
	df_RBD_frequencies = df_peakov_selected * df_RBD_status.loc[df_peakov_selected.columns.values, "status"]
	RBD_frequencies = np.array(df_RBD_frequencies.sum(axis = 1) / df_peakov_selected.sum(axis = 1))

	return(RBD_frequencies)

