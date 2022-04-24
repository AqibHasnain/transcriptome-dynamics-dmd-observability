import numpy as np
from sklearn.preprocessing import StandardScaler
from scipy.signal import savgol_filter as savgol 
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from copy import deepcopy

def get_groups_from_df(data,labels): 
	''' The df that is imported has the following structure
	each column is a datapoint with the progression 
	1rep1, 1rep2, 1rep3, 1Mrep1, 1Mrep2, 1Mrep3, 2rep1, ..., 12Mrep3
	This function builds the two matrices (one per group)
	1rep1, 1rep2, 1rep3, ..., 12rep1, 12rep2, 12rep3
	1Mrep1, 1Mrep3, 1Mrep3, ..., 12Mrep1, 12Mrep2, 12Mrep3 '''
	data_c = np.zeros([data.shape[0],int(data.shape[1]/2)]) # c for control
	data_t = np.zeros([data.shape[0],int(data.shape[1]/2)]) # t for treatment
	c = 0
	ct = 0
	for i in range(0,data.shape[1]):
		if 'M' not in labels[i]: # 'M' stands for malathion (the treatment of group 2) in the labels
			data_c[:,c] = data[:,i]
			c += 1
		elif 'M' in labels[i]:
			data_t[:,ct] = data[:,i]
			ct += 1
	return data_c, data_t

def get_reps(reps,ntimepts):
	''' get the columns inds of the replicates to keep.'''
	nreps = 3
	allinds = set(list(range(0,ntimepts*nreps)))
	if reps == [0]:
		keepers = list(allinds - set(list(range(1,ntimepts*nreps,3))) - set(list(range(2,ntimepts*nreps,3))))
	elif reps == [1]:
		keepers = list(allinds - set(list(range(0,ntimepts*nreps,3))) - set(list(range(2,ntimepts*nreps,3))))
	elif reps == [2]:
		keepers = list(allinds - set(list(range(0,ntimepts*nreps,3))) - set(list(range(1,ntimepts*nreps,3))))
	elif reps == [0,1]: # set subtract column inds for rep3 from allinds
		keepers = list(allinds - set(list(range(2,ntimepts*nreps,3)))) # column inds for rep1 and rep2
	elif reps == [0,2]: # set subtract column inds for rep2 from allinds
		keepers = list(allinds - set(list(range(1,ntimepts*nreps,3)))) # column inds for rep1 and rep3
	elif reps == [1,2]: # set subtract column inds for rep1 from allinds
		keepers = list(allinds - set(list(range(0,ntimepts*nreps,3))))
	elif reps == [0,1,2]:
		keepers = list(allinds)
	keepers.sort()
	return keepers

def put_groups_in_3D(data,nTraj,nT):
	'''Data from each trajectory (replicate) is placed in a new 2d array which is appended to one 3d array of 
	dimension n x m x r. n is number of genes, m is number of timepoints, r is number of replicates.'''
	X = np.zeros((data.shape[0],nT,nTraj))
	reps = list(range(0,nTraj))
	for i in reps:
		X[:,:,i] = data[:,get_reps([i],nT)]
	return X

def smooth_time_series(data,window_size=5,polyorder=2): 
	'''Using scipy's Savitsky-Golay filter to smooth the data. window_size must be an odd number 
	and greater than polyorder'''
	return savgol(data,window_size,polyorder,axis=1)

def recover_negatives(data):
	data[data < 0.0] = 0.0
	return data

def get_high_exp_genes(data_c, data_t,mean_cutoff=200): 
	mean_c_bool = np.mean(np.mean(data_c,axis=2),axis=1) > mean_cutoff
	mean_t_bool = np.mean(np.mean(data_t,axis=2),axis=1) > mean_cutoff
	keepers_c_high = [ii for ii, val in enumerate(mean_c_bool) if val]
	keepers_t_high = [ii for ii, val in enumerate(mean_t_bool) if val]
	keepers_high = list(set(keepers_c_high) & set(keepers_t_high))
	return keepers_high

def compute_fold_changes(data_c_keep, data_t_keep, pseudoTPM = 20):
	return (data_t_keep + pseudoTPM) / (data_c_keep + pseudoTPM) 

def standardize_time_series(data_fc,ntimepts,nreps):
	# performing standardization across time and replicates
	data_fc = data_fc.reshape(len(data_fc),data_fc.shape[1]*data_fc.shape[2],order='F')
	scaler = StandardScaler()
	scaler.fit(data_fc.T) 
	data_fc_norm = scaler.transform(data_fc.T).T
	data_fc_norm = data_fc_norm.reshape(len(data_fc),ntimepts,nreps,order='F')
	return data_fc_norm

def process_df(df,sampleLabels,nreps,reps,ntimepts,txIDs,sg_filter=True,standardize=True):
	data_c, data_t = get_groups_from_df(np.array(df),sampleLabels)
	data_c_orig, data_t_orig = put_groups_in_3D(data_c,nreps,ntimepts), put_groups_in_3D(data_t,nreps,ntimepts)
	data_c_orig, data_t_orig = deepcopy(data_c_orig)[:,2:-1,reps], deepcopy(data_t_orig)[:,2:-1,reps]
	if sg_filter: 
	   data_c_orig, data_t_orig = smooth_time_series(data_c_orig), smooth_time_series(data_t_orig) 
	   data_c_orig, data_t_orig = recover_negatives(data_c_orig), recover_negatives(data_t_orig)
	# downselect the timepoints not used in this study
	data_c, data_t = deepcopy(data_c_orig), deepcopy(data_t_orig) 
	keepers = get_high_exp_genes(data_c,data_t)
	data_c_keep, data_t_keep = data_c[keepers], data_t[keepers]
	txIDs_keep = [txIDs[x] for x in keepers]
	data_fc = compute_fold_changes(data_c_keep, data_t_keep)
	if standardize:
		data_fc_norm = standardize_time_series(data_fc,data_c.shape[1],data_c.shape[2])
	else: 
		data_fc_norm = deepcopy(data_fc)
	return data_c_orig, data_t_orig, data_c, data_t, keepers, data_c_keep, data_t_keep, txIDs_keep, data_fc, data_fc_norm

# def process_df(df,sampleLabels,nreps,reps,ntimepts,txIDs,sg_filter=True,standardize=True):
# 	data_c, data_t = get_groups_from_df(np.array(df),sampleLabels)
# 	data_c_orig, data_t_orig = put_groups_in_3D(data_c,nreps,ntimepts), put_groups_in_3D(data_t,nreps,ntimepts)
# 	if sg_filter: 
# 	   data_c_orig, data_t_orig = smooth_time_series(data_c_orig[:,:-1,reps]), smooth_time_series(data_t_orig[:,:-1,reps]) 
# 	   data_c_orig, data_t_orig = recover_negatives(data_c_orig), recover_negatives(data_t_orig)
# 	# downselect the timepoints not used in this study
# 	data_c, data_t = deepcopy(data_c_orig[:,2:]), deepcopy(data_t_orig[:,2:]) 
# 	keepers = get_high_exp_genes(data_c,data_t)
# 	data_c_keep, data_t_keep = data_c[keepers], data_t[keepers]
# 	txIDs_keep = [txIDs[x] for x in keepers]
# 	data_fc = compute_fold_changes(data_c_keep, data_t_keep)
# 	if standardize:
# 		data_fc_norm = standardize_time_series(data_fc,data_c.shape[1],data_c.shape[2])
# 	else: 
# 		data_fc_norm = deepcopy(data_fc)
# 	return data_c_orig, data_t_orig, data_c, data_t, keepers, data_c_keep, data_t_keep, txIDs_keep, data_fc, data_fc_norm

# get gene names and locus tags from transcript ID
def getRecords(fasta_path,genbank_path,transcriptIDs):
	# have to use cds_from_genome.fasta because this is the where the transcriptIDs came from (e.g. lcl|AM181176.4_cds_CAY53368.1_5775)
	fasta_records = list(SeqIO.parse(fasta_path,'fasta')) # full cds_from_genome fasta
	keep_fasta_records = [] # getting records of genes that we have used
	for tx in transcriptIDs:
		for rec in fasta_records:
			if rec.name == tx:
				keep_fasta_records.append(rec)
	# match locus_tags in keep_records (from fasta) with tags in genbank to easily parse rest of cds' description 
	# can grab locus_tags from the fasta description, other info need to grab from the genbank
	genes, locus_tags = [],[] 
	for rec in keep_fasta_records:
		rec_elems = [x.strip().strip(']') for x in rec.description.split(' [')]
		if 'gene=' in str(rec_elems): # sequence has gene name (e.g. gene=dnaA)
			genes.append(rec_elems[1][5:])
			locus_tags.append(rec_elems[2][10:])
		elif 'gene=' not in str(rec_elems): # sequence has no gene name, but has locus tag
			genes.append('N/A')
			locus_tags.append(rec_elems[1][10:])
	gb_records = next(SeqIO.parse(genbank_path,'genbank'))
	
	locations = []
	for tag in locus_tags:
			for feature in gb_records.features:
				if feature.type == 'CDS':
					if feature.qualifiers['locus_tag'][0] == tag: 
						if feature.strand == 1: # 5' -> 3'
							locations.append([feature.location.start.position,feature.location.end.position,1])
						elif feature.strand == -1: # 3' -> 5' (complementary)
							locations.append([feature.location.start.position,feature.location.end.position,-1])            
	
	return genes, locus_tags, locations

def txid2locustag(txid_list,txid,tag_list):
	return tag_list[txid_list.index(txid)]
















