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
		self.mean_auc_threshold = 0.7
		self.mean_shap_threshold = 1
		self.JI_TR_threshold = 1.0 / 3
		self.JI_TR_top6_threshold = 1.0 / 3

	def filter(self, df_filtered_topics, region):

		CLsig_features_list= []
		# filter condensate-like features
		for index, row in df_filtered_topics.iterrows():
			SigName = row["topic_name"]
			focus_TR = "_".join(SigName.split("_")[:-1])
			sig_idx = SigName.split("_")[-1]
			
			ROC_file = "ROC/{0}_{1}_{2}_AUROC.txt" . format(
				focus_TR,
				region,
				sig_idx
			)

			# XGB_file = "XGBoost/{0}_{1}_{2}_shap.txt" . format(
			# 	focus_TR,
			# 	region,
			# 	sig_idx
			# )

			df_ROC = pd.read_csv(ROC_file, header = 0, sep = "\t", index_col = "feature")
			# df_XGB = pd.read_csv(XGB_file, header = 0, sep = "\t", index_col = "feature")

			condensate_features = ["PPI", "RNA binding strength", "LLPS capacity", "MLO", "Disordered domain", "RNA binding domain"]
			condensate_features = np.intersect1d(df_ROC.index.values, condensate_features) # to avoid unavailable RNA binding strength
			
			roc_status = df_ROC.loc[condensate_features, "AUROC"].astype("float64")
			# xgb_status = df_XGB.loc[condensate_features, "mean_shap"].astype("float64")
			
			# condensate-like feature positive corresponds to signatures
			condensate_features_pos = roc_status[roc_status >= self.individual_auc_threshold].index.values
			if len(condensate_features_pos) >= 3:
				# top 3 condensate-like features positively respond to signature 
				condensate_features_pos_top3 = roc_status[condensate_features_pos].nlargest(3).index.values
				# if roc_status[condensate_features_pos_top3].mean() >= self.mean_auc_threshold and xgb_status[condensate_features_pos_top3].mean() >= self.mean_shap_threshold:
				# 	CLsig_features_list.append([SigName, ",".join(condensate_features_pos), len(condensate_features_pos), roc_status[condensate_features_pos].mean(), xgb_status[condensate_features_pos].mean()])
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
		# sort signatures based on topic coherence to keep the one with high topic coherence when redudant signatures are found
		# df_filtered_topics_CLsigs = df_filtered_topics_CLsigs.sort_values(by = "topic_coherence", ascending = False)
		# sort signatures based on count and level of enriched condensation features
		df_CLsig_features_sorted = df_CLsig_features.sort_values(by = ["qualified_CL_features_count", "mean_AUROC"], ascending = False)

		list_JI_sigs = []
		for sig_i, sig_j in combinations(df_CLsig_features_sorted.index.values, 2):
			TR_sig_i = df_filtered_topics_CLsigs.loc[sig_i, "component_word"].split(",")
			TR_sig_j = df_filtered_topics_CLsigs.loc[sig_j, "component_word"].split(",")
			# calculate JI between sig-i-pos TRs and sig-j-pos TRs
			JI_TR = len(np.intersect1d(TR_sig_i, TR_sig_j)) / len(np.union1d(TR_sig_i, TR_sig_j))
			
			# calculate JI between top 6 sig-i-pos TRs and sig-j-pos TRs
			JI_TR_top6 = len(np.intersect1d(TR_sig_i[:6], TR_sig_j[:6])) / len(np.union1d(TR_sig_i[:6], TR_sig_j[:6]))
			
			# calculate JI between sig-i-pos and sig-j-pos sites
			focus_TR_i = "_".join(sig_i.split("_")[:-1])
			sig_idx_i = sig_i.split("_")[-1]
			focus_TR_j = "_".join(sig_j.split("_")[:-1])
			sig_idx_j = sig_j.split("_")[-1]
			
			list_JI_sigs.append([sig_i, sig_j, JI_TR, JI_TR_top6])

		df_JI_sigs = pd.DataFrame(list_JI_sigs)
		df_JI_sigs.columns = ["primary", "secondary", "JI_TR", "JI_TR_top6"]

		df_JI_sigs_grouped = df_JI_sigs.groupby("primary")
		redundant_sigs = []
		for group, df_grpup in df_JI_sigs_grouped:
			if group not in redundant_sigs:
				redundant_sigs += list(df_grpup.loc[df_grpup["JI_TR_top6"] > self.JI_TR_top6_threshold, "secondary"].values)
		df_unique_CLsig = df_filtered_topics_CLsigs.loc[np.setdiff1d(df_filtered_topics_CLsigs.index.values, redundant_sigs), :]
		# df_unique_CLsig_featues = df_CLsig_features.loc[df_unique_CLsig.index.values, ["qualified_CL_features", "qualified_CL_features_count", "mean_AUROC", "mean_feature_importance"]]
		df_unique_CLsig_featues = df_CLsig_features.loc[df_unique_CLsig.index.values, ["qualified_CL_features", "qualified_CL_features_count", "mean_AUROC"]]

		df_unique_CLsig_anno = pd.concat([df_unique_CLsig, df_unique_CLsig_featues], axis = 1)
		df_unique_CLsig_anno = df_unique_CLsig_anno.sort_values(by = ["qualified_CL_features_count", "mean_AUROC"], ascending = False)
		# df_unique_CLsig_anno.columns = ['signature_name', 'component_CAP', 'qualified_CL_features', 'qualified_CL_features_count', 'mean_AUROC', 'mean_feature_importance']
		df_unique_CLsig_anno.columns = ['signature_name', 'component_CAP', 'qualified_CL_features', 'qualified_CL_features_count', 'mean_AUROC']
		
		df_unique_CLsig_anno.to_csv("CondSig/{0}_{1}_unique_condensate_like_signatures.txt" . format(self.args.name, region), header = True, index = False, sep = "\t")


class SigEval():
	"""Evaluation of condensate-like capacity for each signature was performed by checking 
	whether signature-positive sites showed higher condensate-like features than signature-negative sites"""

	def __init__(self, args):

		self.args = args

	def run(self, SigName, region):

		focus_TR = "_".join(SigName.split("_")[:-1])
		sig_idx = SigName.split("_")[-1]

		df_sig_pos_FE = pd.read_csv("Features/{0}_{1}_{2}_pos_FE.txt" . format(
				focus_TR,
				region,
				sig_idx), header = 0, sep = "\t")

		df_sig_neg_FE = pd.read_csv("Features/{0}_{1}_{2}_neg_FE.txt" . format(
				focus_TR,
				region,
				sig_idx), header = 0, sep = "\t")


		# balance the number of signature-positive and -negative sites
		# if round(df_sig_neg_FE.shape[0] / df_sig_pos_FE.shape[0]) > 1:
		# 	df_sig_pos_FE = pd.concat([df_sig_pos_FE] * round(df_sig_neg_FE.shape[0] / df_sig_pos_FE.shape[0]), ignore_index = True)
		# elif round(df_sig_pos_FE.shape[0] / df_sig_neg_FE.shape[0]) > 1:
		# 	df_sig_neg_FE = pd.concat([df_sig_neg_FE] * round(df_sig_pos_FE.shape[0] / df_sig_neg_FE.shape[0]), ignore_index = True)

		df_sig_pos_FE.loc[:, "label"] = True
		df_sig_neg_FE.loc[:, "label"] = False

		df_sig_merged_FE_raw = pd.concat([df_sig_pos_FE, df_sig_neg_FE], axis = 0)
		df_sig_merged_FE = df_sig_merged_FE_raw.loc[:, np.setdiff1d(df_sig_merged_FE_raw.columns, ["chrom", "start", "end"])]

		X = df_sig_merged_FE.loc[:, np.setdiff1d(df_sig_merged_FE.columns, ["label"])]
		y = df_sig_merged_FE.loc[:, "label"]

		self.roc_analysis(X, y, focus_TR, region, sig_idx)
		# self.XGB_analysis(X, y, focus_TR, region, sig_idx)

	def roc_analysis(self, X, y, focus_TR, region, sig_idx):
		"""perform ROC analysis for features in X to label y"""
		
		# calculate AUROC and plot ROC curve
		list_auc = []
		for feature in X.columns:
			x = X.loc[:, feature].values
			# fpr, tpr, thresholds = metrics.roc_curve(y, x)
			auc = round(metrics.roc_auc_score(y,x), 2)
			# ax.plot(fpr, tpr, label = "%s (AUC = %0.2f)"% (feature, auc))
			list_auc.append([feature, auc])
		# ax.plot([0, 1], [0, 1], "r--")
		# ax.set_xlim([0.0, 1.0])
		# ax.set_ylim([0.0, 1.05])
		# ax.set_xlabel("FPR")
		# ax.set_ylabel("TPR")
		# ax.set_title("ROC curve")
		# ax.legend(prop={'size':11}, bbox_to_anchor = (1.2, 0.5))

		df_auc = pd.DataFrame(list_auc)
		df_auc.columns = ["feature", "AUROC"]
		df_auc = df_auc.sort_values(by = "AUROC", ascending = True)
		df_auc.to_csv("ROC/{0}_{1}_{2}_AUROC.txt" . format(
				focus_TR,
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
				focus_TR,
				region,
				sig_idx))
		plt.close()
		logging.getLogger().setLevel(logging.INFO)

	def XGB_analysis(self, X, y, focus_TR, region, sig_idx):
		"""calculate feature importance"""

		xgb = XGBClassifier(max_depth = 5,
					n_estimators = 1000,
					learning_rate = 0.01,
					verbosity = 1,
					n_jobs = self.args.threads,
					scale_pos_weight = 1,
					min_child_weight = 1,
					random_state = 1007,
					subsample = 0.8,
					colsample_bytree = 0.8,
					importance_type = "gain")

		# build xgb classification model
		xgb.fit(X, y)

		# explain the model's predictions using SHAP
		explainer = shap.Explainer(xgb)
		shap_values = explainer(X)

		# calculate mean shap value for each feature
		df_shap = pd.DataFrame([X.columns.values, np.abs(shap_values.values).mean(0)]).T
		df_shap.columns = ["feature", "mean_shap"]
		df_shap = df_shap.sort_values(by = "mean_shap", ascending = True)
		df_shap.to_csv("XGBoost/{0}_{1}_{2}_shap.txt".format(focus_TR, region, sig_idx), header = True, sep = "\t", index = False)

		mycmap = cm.get_cmap("viridis")
		mean_shap_normalized = (df_shap.loc[:, "mean_shap"] - df_shap.loc[:, "mean_shap"].min()) / (df_shap.loc[:, "mean_shap"].max() - df_shap.loc[:, "mean_shap"].min())
		# info(mean_shap_normalized)
		# info(type(mean_shap_normalized))
		my_color = mycmap(mean_shap_normalized.astype("float"))
			
		fig, ax1 = plt.subplots(figsize = (4.5,4), nrows = 1, ncols = 1)
		ax1.barh(df_shap.loc[:, "feature"], df_shap.loc[:, "mean_shap"], color = my_color, height = 0.7)
		ax1.set_xlabel("SHAP feature importance\n(XGBoost classification model)")
		ax1.set_title("Feature importance")
		plt.subplots_adjust(bottom = 0.2, left = 0.45, right = 0.95)
		plt.savefig("XGBoost/{0}_{1}_{2}_FeatureImportance.pdf" . format(
			focus_TR,
			region,
			sig_idx))
		plt.close()


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
		

	def FE_preprocess(self):

		"""preprocessing input feature file"""
		
		# genearate motif presence matrix
		# info("1-1. Generate motif presence matrix ... ")
		
		# motif_path = self.args.motif
		# motif_factor_uid_file = "motif_factor_uids_edited.txt"
		# if not os.path.isfile("{0}/{1}".format(motif_path, motif_factor_uid_file)):
		# 	error("No annotation file in motif path: {0}/{1}" . format(motif_path, motif_factor_uid_file))
		# 	sys.exit(0)
		# df_motif_factor_uid = pd.read_csv("{0}/{1}".format(motif_path, motif_factor_uid_file), header = None, sep = "\t", encoding='latin-1').iloc[:, 0:4]
		# df_motif_factor_uid.columns = ["factor", "label", "uniprot_id", "uniprot_entry"]
		# df_motif_factor_uid.loc[:, "uniprot_id"] = df_motif_factor_uid.loc[:, "uniprot_id"].str.upper()
		
		# df_motif_factor_uid_dataset = df_motif_factor_uid.loc[df_motif_factor_uid["uniprot_id"].isin(self.df_dataset.loc[:, "uniprot_id"].values), :].copy() # only motifs in dataset annotation were selected 
		# df_motif_factor_uid_dataset.loc[:, "file"] = ["{0}/filtered_results/{1}".format(motif_path, motif_label) for motif_label in df_motif_factor_uid_dataset.loc[:, "label"].values]

		# motifov_bed = Utility.Motif_OverlapMatrix(df_motif_factor_uid_dataset, self.args).intervals_intersect_genomewide_bins()
		# df_motif_dataset = pd.read_csv(motifov_bed, header = 0, sep = "\t", index_col = "name").iloc[:, 3:]
		# df_motif_dataset[df_motif_dataset > 1] = 1
		
		# motifov_uid_list = []
		# uid_list = []
		# for uid in np.unique(df_motif_factor_uid_dataset.loc[:, "uniprot_id"].values):
		# 	motif_labels = df_motif_factor_uid_dataset.loc[df_motif_factor_uid_dataset["uniprot_id"] == uid, "label"] # a uniprot id may be responsible for multiple labels (multiple files for the same factor)
		# 	motifov_uid = df_motif_dataset.loc[:, motif_labels].max(axis = 1) # if any one motif file show co-binding, the factor was denoted to be motif occurrence
		# 	motifov_uid_list.append(motifov_uid)
		# 	uid_list.append(uid)
		# df_motif_dataset_labeled = pd.DataFrame(motifov_uid_list).T
		# df_motif_dataset_labeled.columns = uid_list
		# self.df_motif_dataset_labeled = df_motif_dataset_labeled

		# protein-protein interaction between proteins of dataset (to reduce the execution time of FE run)
		df_PPI_dataset = self.df_PPI_merged.loc[
						((
							self.df_PPI_merged["SWISS-PROT Accessions Interactor A"].isin(self.df_dataset.loc[:, "uniprot_id"].values)
							) & (
								self.df_PPI_merged["SWISS-PROT Accessions Interactor B"].isin(self.df_dataset.loc[:, "uniprot_id"].values)
								)), :]
		self.df_PPI_dataset = df_PPI_dataset

		# disordered domain content
		disorder_content_cutoff = 0.153 # determined by IDR percentage of LLPS proteins in human
		df_mobidb_lite = pd.read_csv(self.args.IDR, header = None, sep = "\t")
		df_mobidb_lite.columns = ["uniprot_id", "disorder_content"]
		df_mobidb_lite.loc[:, "uniprot_id"] = df_mobidb_lite.loc[:, "uniprot_id"].str.upper()
		df_mobidb_lite.loc[:, "status"] = 0
		df_mobidb_lite.loc[df_mobidb_lite["disorder_content"] >= disorder_content_cutoff, "status"] = 1
		df_mobidb_lite.index = df_mobidb_lite.loc[:, "uniprot_id"].values
		self.df_mobidb_lite = df_mobidb_lite

		# RNA binding domain content
		df_RBD = pd.read_csv(self.args.RBD, header = 0, sep = "\t")
		RBD_uids = df_RBD.loc[:, "uid"].dropna().values
		self.RBD_uids = RBD_uids


	def FE_run(self, SigName, region):

		focus_TR = "_".join(SigName.split("_")[:-1])
		sig_idx = SigName.split("_")[-1]

		pos_sites = pd.read_csv(
			"{0}/LearnSig/{1}/{2}/{1}_{2}_{3}_pos_sites.bed" . format(
				self.args.Sig_InputPath,
				focus_TR,
				region,
				sig_idx
				),
			header = None, sep = "\t"
			).iloc[:,3].values

		neg_sites = pd.read_csv(
			"{0}/LearnSig/{1}/{2}/{1}_{2}_{3}_neg_sites.bed" . format(
				self.args.Sig_InputPath,
				focus_TR,
				region,
				sig_idx
				),
			header = None, sep = "\t"
			).iloc[:, 3].values

		potential_combinatorial_TRs = pd.read_csv(
			"{0}/LearnSig/{1}/{1}_{2}_voca.txt" . format(
				self.args.Sig_InputPath,
				focus_TR,
				region
				),
			header = None, sep = "\t"
			).iloc[:, 1].values

		df_peaks_sig_pos = self.df_peakov.loc[pos_sites, ["chrom", "start",	"end"]].copy()
		df_peaks_sig_neg = self.df_peakov.loc[neg_sites, ["chrom", "start",	"end"]].copy()

		df_peakov_sig_pos = self.df_peakov.loc[pos_sites, potential_combinatorial_TRs].copy()
		df_peakov_sig_neg = self.df_peakov.loc[neg_sites, potential_combinatorial_TRs].copy()

		# convert factor label to uniprot ids
		df_peakov_sig_pos_uids = df_peakov_sig_pos.copy()
		df_peakov_sig_neg_uids = df_peakov_sig_neg.copy()
		df_peakov_sig_pos_uids.columns = self.df_dataset.loc[potential_combinatorial_TRs, "uniprot_id"].values
		df_peakov_sig_neg_uids.columns = self.df_dataset.loc[potential_combinatorial_TRs, "uniprot_id"].values

		
		# chromatin accessiblity
		# CA_pos = self.CA_FE(df_peaks_sig_pos, self.args.accessibility, self.args.threads)
		# CA_neg = self.CA_FE(df_peaks_sig_neg, self.args.accessibility, self.args.threads)

		# distance to TSS
		# cmd_pos = "bedtools closest -d -t first -a {0} -b {1} | cut -f 10" . format("{0}/LearnSig/{1}/{2}/{1}_{2}_{3}_pos_sites.bed" . format(self.args.Sig_InputPath,
		# 		focus_TR,
		# 		region,
		# 		sig_idx), self.tss_file)
		# cmd_neg = "bedtools closest -d -t first -a {0} -b {1} | cut -f 10" . format("{0}/LearnSig/{1}/{2}/{1}_{2}_{3}_neg_sites.bed" . format(self.args.Sig_InputPath,
		# 		focus_TR,
		# 		region,
		# 		sig_idx), self.tss_file)

		
		# p_pos = subprocess.Popen(
		#   cmd_pos, 
		#   stdout = subprocess.PIPE,
		#   stderr = subprocess.PIPE,
		# 	shell = "True")
		# (stdoutdata_pos, stderrdata_pos) = p_pos.communicate()
		# exit_code = p_pos.returncode
		# distance_pos = list(map(int, stdoutdata_pos.decode().split("\n")[:-1]))
		# add_sudocount = lambda x:x+1
		# distance_pos = list(map(add_sudocount, distance_pos))
		# distance_pos_log10 = -np.log10(distance_pos)

		# p_neg = subprocess.Popen(
		#   cmd_neg, 
		#   stdout = subprocess.PIPE,
		#   stderr = subprocess.PIPE,
		# 	shell = "True")
		# (stdoutdata_neg, stderrdata_neg) = p_neg.communicate()
		# exit_code = p_neg.returncode
		# distance_neg = list(map(int, stdoutdata_neg.decode().split("\n")[:-1]))
		# distance_neg = list(map(add_sudocount, distance_neg))
		# distance_neg_log10 = -np.log10(distance_neg)

		df_peakov_sig_pos_uids_splits = np.array_split(df_peakov_sig_pos_uids, self.args.threads)
		df_peakov_sig_neg_uids_splits = np.array_split(df_peakov_sig_neg_uids, self.args.threads)

		# motif presence
		# focused_motif_uid = np.intersect1d(df_peakov_sig_pos_uids.columns.values, self.df_motif_dataset_labeled.columns.values)
		# df_motif_sig_pos = self.df_motif_dataset_labeled.loc[pos_sites, focused_motif_uid].copy()
		# df_motif_sig_neg = self.df_motif_dataset_labeled.loc[neg_sites, focused_motif_uid].copy()
		
		# pool = multiprocessing.Pool(self.args.threads)
		# func_motif_pos = partial(motif_frequency_cal, df_motif_sig_pos)
		# motif_frequency_pos = np.concatenate(pool.map(func_motif_pos, df_peakov_sig_pos_uids_splits))
		# pool.close()
		# pool.join()

		# pool = multiprocessing.Pool(self.args.threads)
		# func_motif_neg = partial(motif_frequency_cal, df_motif_sig_neg)
		# motif_frequency_neg = np.concatenate(pool.map(func_motif_neg, df_peakov_sig_neg_uids_splits))
		# pool.close()
		# pool.join()

		# protein-protein interaction frequency
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
		func_LLPS = partial(LLPS_frequency_cal, self.df_LLPS)
		pool = multiprocessing.Pool(self.args.threads)
		LLPS_frequency_pos = np.concatenate(pool.map(func_LLPS, df_peakov_sig_pos_uids_splits))
		pool.close()
		pool.join()

		pool = multiprocessing.Pool(self.args.threads)
		LLPS_frequency_neg = np.concatenate(pool.map(func_LLPS, df_peakov_sig_neg_uids_splits))
		pool.close()
		pool.join()

		# MLO frequency
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
		func_disorder = partial(disorder_frequency_cal, self.df_mobidb_lite)
		pool = multiprocessing.Pool(self.args.threads)
		disorder_pos = np.concatenate(pool.map(func_disorder, df_peakov_sig_pos_uids_splits))
		pool.close()
		pool.join()
		
		pool = multiprocessing.Pool(self.args.threads)
		disorder_neg = np.concatenate(pool.map(func_disorder, df_peakov_sig_neg_uids_splits))
		pool.close()
		pool.join()

		# RNA-binding domain frequency
		func_RBD = partial(RBD_frequency_cal, self.RBD_uids)
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
				# pd.DataFrame(CA_pos), 
				# pd.DataFrame(motif_frequency_pos, index = df_peakov_sig_pos_uids.index.values), 
				# pd.DataFrame(distance_pos_log10, index = df_peakov_sig_pos_uids.index.values), 
				pd.DataFrame(PPI_frequency_pos, index = df_peakov_sig_pos_uids.index.values),
				pd.DataFrame(RNA_pos),
				pd.DataFrame(LLPS_frequency_pos, index = df_peakov_sig_pos_uids.index.values),
				pd.DataFrame(MLO_frequency_pos, index = df_peakov_sig_pos_uids.index.values),
				pd.DataFrame(disorder_pos, index = df_peakov_sig_pos_uids.index.values),
				pd.DataFrame(RBD_pos, index = df_peakov_sig_pos_uids.index.values)
				], axis = 1)

			df_sig_pos_FE.index = df_peaks_sig_pos.index.values
			# df_sig_pos_FE.columns = ["chrom", "start", "end", "Chromatin accessibility", "Motif presence", "Proximity to TSS", "PPI", "RNA binding strength", "LLPS capacity", "MLO", "Disordered domain", "RNA binding domain"]
			df_sig_pos_FE.columns = ["chrom", "start", "end", "PPI", "RNA binding strength", "LLPS capacity", "MLO", "Disordered domain", "RNA binding domain"]

			df_sig_neg_FE = pd.concat([df_peaks_sig_neg, 
				# pd.DataFrame(CA_neg), 
				# pd.DataFrame(motif_frequency_neg, index = df_peakov_sig_neg_uids.index.values), 
				# pd.DataFrame(distance_neg_log10, index = df_peakov_sig_neg_uids.index.values), 
				pd.DataFrame(PPI_frequency_neg, index = df_peakov_sig_neg_uids.index.values),
				pd.DataFrame(RNA_neg),
				pd.DataFrame(LLPS_frequency_neg, index = df_peakov_sig_neg_uids.index.values),
				pd.DataFrame(MLO_frequency_neg, index = df_peakov_sig_neg_uids.index.values),
				pd.DataFrame(disorder_neg, index = df_peakov_sig_neg_uids.index.values),
				pd.DataFrame(RBD_neg, index = df_peakov_sig_neg_uids.index.values)
				], axis = 1)

			df_sig_neg_FE.index = df_peaks_sig_neg.index.values
			df_sig_neg_FE.columns = ["chrom", "start", "end", "PPI", "RNA binding strength", "LLPS capacity", "MLO", "Disordered domain", "RNA binding domain"]

			df_sig_pos_FE.to_csv("Features/{0}_{1}_{2}_pos_FE.txt" . format(
				focus_TR,
				region,
				sig_idx), header = True, sep = "\t", index = False)

			df_sig_neg_FE.to_csv("Features/{0}_{1}_{2}_neg_FE.txt" . format(
				focus_TR,
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
			df_sig_pos_FE.columns = ["chrom", "start", "end", "PPI", "LLPS capacity", "MLO", "Disordered domain", "RNA binding domain"]

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
			df_sig_neg_FE.columns = ["chrom", "start", "end", "PPI", "LLPS capacity", "MLO", "Disordered domain", "RNA binding domain"]

			df_sig_pos_FE.to_csv("Features/{0}_{1}_{2}_pos_FE.txt" . format(
				focus_TR,
				region,
				sig_idx), header = True, sep = "\t", index = False)

			df_sig_neg_FE.to_csv("Features/{0}_{1}_{2}_neg_FE.txt" . format(
				focus_TR,
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

def motif_frequency_cal(df_motif_selected, df_peakov_selected):
	"""caluclate the ratio of occupancy events with motif presence"""

	motif_frequencies = []
	for index, row_peakov in df_peakov_selected.iterrows():
		peakov_factor = row_peakov[row_peakov > 0].index.values # factor with occupancy events
		if index in df_motif_selected.index.values:
			row_motif = df_motif_selected.loc[index, :]
			motif_factor = row_motif[row_motif > 0].index.values # factor with motif presence
			motif_peakov_number = len(np.intersect1d(motif_factor, peakov_factor)) # factor with both occupancy events and motif presence
			motif_frequency = motif_peakov_number / len(peakov_factor)
		else:
			motif_frequency = 0
		motif_frequencies.append(motif_frequency)

	return(motif_frequencies)

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

def LLPS_frequency_cal(df_LLPS, df_peakov_selected):
	"""calculate LLPS frequency of occupied factors"""
		
	LLPS_frequencies = []
	LLPS_protein_uids = df_LLPS.loc[:, "UniProt ID"].dropna().values
	for index, row in df_peakov_selected.iterrows():
		peakov_uids = row[row > 0].index.values
		peakov_uids_LLPS = np.intersect1d(peakov_uids, LLPS_protein_uids)
		LLPS_frequency = len(peakov_uids_LLPS) / len(peakov_uids)
		LLPS_frequencies.append(LLPS_frequency)
		
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

def disorder_frequency_cal(df_mobidb_lite, df_peakov_selected):
		
	disorder_frequencies = []
	for index, row in df_peakov_selected.iterrows():
		peakov_uids = row[row > 0].index.values
		disorder_frequency = np.sum(df_mobidb_lite.loc[peakov_uids, "status"]) / len(peakov_uids)
		disorder_frequencies.append(disorder_frequency)
	return(disorder_frequencies)

def RBD_frequency_cal(RBD_uids, df_peakov_selected):
		
	RBD_frequencies = []
	for index, row in df_peakov_selected.iterrows():
		peakov_uids = row[row > 0].index.values
		RBD_frequency = len(np.intersect1d(peakov_uids, RBD_uids)) / len(peakov_uids)
		RBD_frequencies.append(RBD_frequency)
	return(RBD_frequencies)

