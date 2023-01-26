### import modules
from argparse import ArgumentTypeError
import argparse
import os
from Bio import SeqIO
import re
import numpy as np
import random
import khmer
import warnings
import pandas as pd
import copy
import seaborn as sb
import matplotlib.pyplot as plt


from script.benchmark_1_pilot_ext_randstrobe_vs_regular_randstrobe import generate_string_list_from_genome, generate_mutated_strings, count_island_from_binary_vector

def is_kmer_open_syncmer_by_hash_value(kmer, s, shift=0, max_prime=9999999999971., rev_comp=False):
	"""
	Determines whether a kmer is an open syncmer with position shift
	:param kmer: input kmer
	:param s: submer length
	:param shift: location of the min submer (default 0, which is the start)
	:param rev_comp: use canonical submer
	:return: boolean T/F
	"""
	target_submer = kmer[shift:shift + s]
	if rev_comp:
		target_submer = min(target_submer, khmer.reverse_complement(target_submer))
	
	target_value = khmer.hash_no_rc_murmur3(target_submer) % max_prime
	
	for start_pos in list(range(len(kmer) - s + 1)):
		other_submer = kmer[start_pos:start_pos + s]
		if rev_comp:
			other_submer = min(other_submer, khmer.reverse_complement(other_submer))
		if target_value > khmer.hash_no_rc_murmur3(other_submer) % max_prime:
			return False
	
	return True

def write_or_extend_kmer_record_to_dict(input_dict, item_key, item_value):
	"""
	For a k-mer record, add to a dict if it's not there; or append its location if it already exists

	:param input_dict: the dict to be update
	:param item_key: the key of the item (k-mer) to add
	:param item_value: the value of the item (k-mer) to add
	:return: no direct return, but the input_dict will be updated
	"""
	if item_key in input_dict:
		input_dict[item_key].append(item_value)
	else:
		input_dict[item_key] = [item_value]

def get_syncmer_eval_matrix(seq_length: int, syncmer_length: int, overlap_pos_list: list, full_pos_list_length: int):
	"""
	Get eval matrix: m, sc, mc, E for syncmer matchs
	:param seq_length: length of the sequence
	:param syncmer_length: length of syncmer
	:param overlap_pos_list: list of start positions in original dict that overlap with tested dict
	:param full_pos_list_length: sketch size of all syncmers
	"""
	### m: % of matched substring
	m = 100.0 * len(overlap_pos_list) / full_pos_list_length
	
	### sc and mc are the same for syncmers
	sc_binary = np.array([0] * seq_length)
	for start_pos in overlap_pos_list:
		sc_binary[start_pos:start_pos + syncmer_length] = 1
	
	sc = 100.0 * sc_binary.sum() / seq_length
	mc = sc
	
	### island
	e_size = count_island_from_binary_vector(sc_binary)
	
	return [m, sc, mc, e_size]

def get_rs_ext_sync_eval_matrix(input_object: object, overlap_rs_mers_list: list):
	"""
	Get eval matrix: m, sc, mc, E for randstrobe-syncmers (abbv: rs_mer)
	:param input_object: original object to load
	:param overlap_rs_mers_list: list of overlapping rs_mers between original object and mutated objects
	"""
	if len(overlap_rs_mers_list) == 0:
		return [0,0,0,len(input_object.seq_string)]
	
	### calculate total rs_mer size: don't count dup
	full_rs_mer_size = len(input_object.rs_syncmer_dict)
	
	### m: % of matched substring
	m = 100.0 * len(overlap_rs_mers_list) / full_rs_mer_size
	
	### sc and mc: coverage of 2 syncmers in rs_mer
	seq_length = len(input_object.seq_string)
	rs_mer_length = len(overlap_rs_mers_list[0])
	one_side_ext_length = int(rs_mer_length/4)  # 2 syncmer * 2 sides
	
	sc_binary = np.array([0] * seq_length)
	mc_binary = np.array([0] * seq_length)
	
	for overlapped_rs_mer in overlap_rs_mers_list:
		sync_location_list = input_object.rs_syncmer_dict[overlapped_rs_mer]
		# this list may have multiple hits
		for item in sync_location_list:
			sync1_mid = item[0]
			sc_binary[sync1_mid - one_side_ext_length : sync1_mid + one_side_ext_length] = 1
			
			sync2_mid = item[1]
			sc_binary[sync2_mid - one_side_ext_length: sync2_mid + one_side_ext_length] = 1
			
			# mc is from start to end
			mc_binary[sync1_mid - one_side_ext_length : sync2_mid + one_side_ext_length] = 1
	
	sc = 100.0 * sc_binary.sum() / seq_length
	mc = 100.0 * mc_binary.sum() / seq_length
	
	### island
	e_size = count_island_from_binary_vector(sc_binary)
	
	return [m, sc, mc, e_size]


### kmer object for syncmer
class ext_syncmer(object):
	def __init__(self, seq_string, k_max=30, k_min=10, trunc_step=6, submer_length_dif=4):
		self.k_max = k_max
		self.k_min = k_min
		self.trunc_step = trunc_step
		self.submer_length_dif = submer_length_dif  # to ensure same syncmer sketch size during trunction
		self.seq_string = seq_string  # sequence to generate k-mers
		
		# only need a list for k_min syncmers
		self.k_min_sync_location=[]
		self.get_kmin_syncmer_location_list()
		
		# ext-syncmer can be directly retrived from input sequence
		self.current_ext_syncmer_length = 0
		self.ext_syncmer_dict = {}  # need location infor to get quality matrix
		
		# regular syncmer dict
		self.regular_syncmer_dict={}  # need location infor to get quality matrix
		
	def get_kmin_syncmer_location_list(self):
		"""
		Generate list of locations of k_min-based syncmers for future extension
		"""
		seq = self.seq_string
		k_value = self.k_min
		submer_length = k_value - self.submer_length_dif
		shift = 0
		location_list = self.k_min_sync_location
		
		for i in range(len(seq) - k_value + 1):
			kmer = seq[i:i + k_value]
			if is_kmer_open_syncmer_by_hash_value(kmer=kmer, s=submer_length, shift=shift):
				location_list.append(i)
				
	def get_ext_syncmer_dict(self, ext_ksize: int):
		if ext_ksize < self.k_min:
			raise Exception("Ext ksize should be no less than k_min: %s" %self.k_min)
		
		k_min_syncmer_location = self.k_min_sync_location
		seq = self.seq_string
		self.current_ext_syncmer_length = ext_ksize # update record
		# clear prior data when rebuild
		self.ext_syncmer_dict = {}
		out_dict = self.ext_syncmer_dict
		
		# save all syncmers to a dict with kmer:loc pairs
		for start_pos in k_min_syncmer_location:
			cur_kmer = seq[start_pos:start_pos + ext_ksize]
			write_or_extend_kmer_record_to_dict(input_dict=out_dict, item_key=cur_kmer, item_value=start_pos)
			
	def build_regular_syncmer_dict(self, k_value):
		"""
		Build regular syncmer dict, similar to kmin_sync_list function above
		"""
		seq = self.seq_string
		k_value = k_value
		submer_length = k_value - self.submer_length_dif
		shift = 0
		# clear prior data when run again
		self.regular_syncmer_dict = {}
		out_dict = self.regular_syncmer_dict
		
		for i in range(len(seq) - k_value + 1):
			kmer = seq[i:i + k_value]
			if is_kmer_open_syncmer_by_hash_value(kmer=kmer, s=submer_length, shift=shift):
				write_or_extend_kmer_record_to_dict(input_dict=out_dict, item_key=kmer, item_value=i)
			
		
		
### kmer object for rand-mid-shift-syncmers
# note:
# 1. maintain a syncmer list for strobe selection (sample from range of strobes, not in nt location)
# 2. also maintain a dict to recover match coverage
# 3. for location of syncmer, save mid point instead of starting point (for ext purpose)

class rand_sync(object):
	def __init__(self, seq_string, k_max=30, k_min=10, submer_length_dif=4, wmin=5, wmax=11):
		self.k_max = k_max
		self.k_min = k_min
		if submer_length_dif % 2 != 0:
			raise Exception("For symmetry, the submer length dif must be even")
		self.submer_length_dif = submer_length_dif  # to ensure same syncmer sketch size during trunction
		self.shift = int(self.submer_length_dif / 2)
		self.seq_string = seq_string  # sequence to generate k-mers
		self.wmin = wmin
		self.wmax = wmax
		
		# only need a list for k_min syncmers
		self.k_min_sync_location = []
		self.get_kmin_syncmer_location_list()
		
		# ext-syncmer
		self.ext_syncmer_list = []
		self.ext_syncmer_loc_list = [] # deal with rep syncmers to record location of every single syncmer
		
		# regular syncmer
		self.regular_syncmer_list = []
		self.regular_syncmer_loc_list = []
		
		# randstrobe-syncmer
		self.rs_syncmer_dict = {} # need location infor to get quality
	
	def get_kmin_syncmer_location_list(self):
		"""
		Generate list of locations of k_min-based syncmers for future extension
		"""
		seq = self.seq_string
		k_value = self.k_min
		submer_length = k_value - self.submer_length_dif
		shift = self.shift
		location_list = self.k_min_sync_location
		
		for i in range(len(seq) - k_value + 1):
			kmer = seq[i:i + k_value]
			if is_kmer_open_syncmer_by_hash_value(kmer=kmer, s=submer_length, shift=shift, rev_comp=True):
				location_list.append(int(i+ k_value/2))
	
	def build_regular_syncmer_list(self, k_value):
		"""
		Build regular syncmer list, and store location in an affilinary list
		"""
		seq = self.seq_string
		k_value = k_value
		submer_length = k_value - self.submer_length_dif
		shift = self.shift
		# clear prior data when run again
		self.regular_syncmer_list = []
		self.regular_syncmer_loc_list = []
		# output
		sync_list = self.regular_syncmer_list
		sync_loc_list = self.regular_syncmer_loc_list
		
		for i in range(len(seq) - k_value + 1):
			kmer = seq[i:i + k_value]
			if is_kmer_open_syncmer_by_hash_value(kmer=kmer, s=submer_length, shift=shift, rev_comp=True):
				sync_list.append(kmer)
				# save mid point as location to be in same manner as ext_rs_mer for countint seq coverage
				sync_loc_list.append(int(i+ k_value/2))
				
	def build_ext_syncmer_list(self, ext_ksize):
		if ext_ksize < self.k_min:
			raise Exception("Ext ksize should be no less than k_min: %s" % self.k_min)
		
		if (ext_ksize - self.k_min) % 2 != 0:
			raise Exception("Ext_ksize - kmin should be even to extend from both sides.")
		
		k_min_syncmer_location = self.k_min_sync_location
		seq = self.seq_string
		# clear prior data when rebuild
		self.ext_syncmer_list = []
		self.ext_syncmer_loc_list = []
		# output
		out_ext_list = self.ext_syncmer_list
		out_ext_loc_list = self.ext_syncmer_loc_list
		
		# note: extension is for 2 sides
		ext_length = int(ext_ksize / 2) # one side ext
		total_length = len(seq)
		
		for start_pos in k_min_syncmer_location:
			# use min/max to handle out_of_boundary situation
			cur_kmer = seq[ max(0, start_pos - ext_length) : min(start_pos + ext_length, total_length) ]
			out_ext_list.append(cur_kmer)
			# this is already mid point (calculated from previous function)
			out_ext_loc_list.append(start_pos)
			
	def build_rand_sync_based_on_syncmer_lists(self, sync_list, loca_list, prime=6229):
		# clear prior data when rebuild
		self.rs_syncmer_dict = {}
		# output
		out_dict = self.rs_syncmer_dict
		wmin = self.wmin
		wmax = self.wmax
		
		# loop through syncmer list
		for index in range(len(sync_list) - wmax):  # last strobe to have full size window width
			strobe1 = sync_list[index]
			
			# select 2nd strobe and record location
			list_strobe2_loc = list(range(index+wmin, index+wmax))
			
			current_hash = prime + 1
			current_pos = 0
			
			for pos in list_strobe2_loc:
				current_syncmer = sync_list[pos]
				new_hash = ( khmer.hash_no_rc_murmur3(strobe1) + khmer.hash_no_rc_murmur3(current_syncmer) ) % prime
				# find argmin
				if new_hash < current_hash:
					current_hash = new_hash
					current_pos = pos
					
			# store result
			strobe2 = sync_list[current_pos]
			rand_syncmer = strobe1 + strobe2
			# syncmer location is stored in loca_list
			loc_rand_syncmer = [loca_list[index], loca_list[current_pos]]
			
			# may encounter dup, so each rand_syncmer use a list of 2 location values
			write_or_extend_kmer_record_to_dict(input_dict=out_dict, item_key=rand_syncmer, item_value = loc_rand_syncmer)
				
		
			
		


### local tests
def local_variables():
	genome_file = os.path.abspath("./input/toy_data.fna")
	sample_size = 10
	k_max = 30
	k_list = list(range(k_max, 10, -4))
	
	p_range = "0.01,0.03,0.05,0.08,0.1"
	p_values = [float(x) for x in p_range.split(",")]
	seq_length = 10000
	wmin=5
	wmax=11
	
	


### run analysis
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Strobemer simulation",
	                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-g', '--genome', help="Path to input genomes")
	parser.add_argument('-k', '--maxk', type=int, help="Max k value to start", default="30")
	parser.add_argument('-n', '--sample_size', type=int, help="Number of technical replicates", default=10)
	parser.add_argument('-p', '--p_range', type=str, help="mutation rate range",
	                    default="0.01,0.03,0.05,0.08,0.1")
	parser.add_argument('-w', '--window_min', type=int, help="Min window size", default=10)
	parser.add_argument('-e', '--window_max', type=int, help="Max window size", default=20)
	parser.add_argument('-l', '--seq_length', type=int, help="Seq length for the test", default=10000)
	
	########################################## read parameters
	args = parser.parse_args()
	
	if not os.path.exists(args.genome):
		raise Exception("Input file %s does not exist." % args.genome)
	genome_file = os.path.abspath(args.genome)
	print("Input genome file is " + genome_file)
	
	sample_size = args.sample_size
	print("Permutation sample size is " + str(sample_size))
	
	k_max = args.maxk
	print("The max k value is ")
	print(k_max)
	k_list = list(range(k_max, 10, -4))
	
	p_range = args.p_range
	p_values = [float(x) for x in p_range.split(",")]
	print("The range of mutation rate p values are: ")
	print(p_values)
	
	seq_length = args.seq_length
	print("The sequence length to keep is %s" % seq_length)
	wmin = args.window_min
	print("Window size min is %s" % wmin)
	wmax = args.window_max
	print("Window size max is %s" % max)
	
	
	########################################## Step1: prepare genome
	not_unknown = re.compile('[nN]')
	# read genome files into a long string
	raw_seq = generate_string_list_from_genome(genome_file)
	cleaned_seq_list = []
	for fragment in raw_seq:
		fragment = fragment.upper()
		frag_split_onlyACTG = not_unknown.split(fragment)
		cleaned_seq_list.extend([x for x in frag_split_onlyACTG if len(x) > 0])
	
	### for this test, we only need 1 string not the multiple records per chr, so manually concatenate it
	out_seq = "".join(cleaned_seq_list)[:seq_length]
	print("The seq length is " + str(len(out_seq)))
	
	
	
	
	### benchmark 2: test syncmer truncation
	original_obj = ext_syncmer(seq_string=out_seq, k_max=k_max, k_min=10, submer_length_dif=4)
	print("Current ext-syncmer size is " + str(original_obj.current_ext_syncmer_length))
	
	for mutation_rate_p in p_values:
		print("Simuation with mutation rate %s" % mutation_rate_p)
		# get a list of sequences
		mut_obj_list = []
		for i in range(sample_size):
			temp_seq = generate_mutated_strings(out_seq, p=mutation_rate_p)
			mut_obj_list.append(ext_syncmer(seq_string=temp_seq, k_max=k_max, k_min=10, submer_length_dif=4))
			
			
		### compare regular syncmer
		out_dict = {}
		
		for k in k_list:
			out_eval_list = [0,0,0,0]
			original_obj.build_regular_syncmer_dict(k_value=k)
			for mut_obj in mut_obj_list:
				mut_obj.build_regular_syncmer_dict(k_value=k)
				overlap_keys = original_obj.regular_syncmer_dict.keys() & mut_obj.regular_syncmer_dict.keys()
				overlap_list = []
				for kmer in overlap_keys:
					# because every loc is a list (to deal with dup), so use extend
					overlap_list.extend(original_obj.regular_syncmer_dict[kmer])
					
				eval_matrix = get_syncmer_eval_matrix(seq_length=len(original_obj.seq_string),
				                                      syncmer_length=k,
				                                      overlap_pos_list=overlap_list,
				                                      full_pos_list_length=len(original_obj.regular_syncmer_dict))
				
				out_eval_list = [sum(x) for x in zip(out_eval_list, eval_matrix)]
			
			# output for this k value
			out_eval_list = [x/sample_size for x in out_eval_list]
			out_dict[str(k)] = out_eval_list
		
		# save to file
		out_df = pd.DataFrame(out_dict).transpose()
		out_df.columns = ['m', 'sc', 'mc', 'E']
		filename = "_".join(["regular_syncmer_mut_rate", str(mutation_rate_p), ".csv"])
		out_df.to_csv(filename, header=True)
			
		
		### compare ext-syncmer
		out_dict = {}
		
		for k in k_list:
			out_eval_list = [0, 0, 0, 0]
			original_obj.get_ext_syncmer_dict(ext_ksize=k)
			for mut_obj in mut_obj_list:
				mut_obj.get_ext_syncmer_dict(ext_ksize=k)
				overlap_keys = original_obj.ext_syncmer_dict.keys() & mut_obj.ext_syncmer_dict.keys()
				overlap_list = []
				for kmer in overlap_keys:
					# because every loc is a list (to deal with dup), so use extend
					overlap_list.extend(original_obj.ext_syncmer_dict[kmer])
				
				eval_matrix = get_syncmer_eval_matrix(seq_length=len(original_obj.seq_string),
				                                      syncmer_length=k,
				                                      overlap_pos_list=overlap_list,
				                                      full_pos_list_length=len(original_obj.ext_syncmer_dict))
				
				out_eval_list = [sum(x) for x in zip(out_eval_list, eval_matrix)]
				
			# output for this k value
			out_eval_list = [x / sample_size for x in out_eval_list]
			out_dict[str(k)] = out_eval_list
			
		# save to file
		out_df = pd.DataFrame(out_dict).transpose()
		out_df.columns = ['m', 'sc', 'mc', 'E']
		filename = "_".join(["ext_syncmer_mut_rate", str(mutation_rate_p), ".csv"])
		out_df.to_csv(filename, header=True)




	### benchmark 3: ext-rand-sync vs regular
	original_obj = rand_sync(seq_string=out_seq, k_max=k_max, k_min=10, submer_length_dif=4, wmin=5, wmax=11)
	
	for mutation_rate_p in p_values:
		print("Simuation with mutation rate %s" % mutation_rate_p)
		# get a list of sequences
		mut_obj_list = []
		for i in range(sample_size):
			temp_seq = generate_mutated_strings(out_seq, p=mutation_rate_p)
			mut_obj_list.append(rand_sync(seq_string=temp_seq, k_max=k_max, k_min=10, submer_length_dif=4, wmin=5, wmax=11))
			
		### compare regular rs-mer
		out_dict = {}
		
		for k in k_list:
			out_eval_list = [0, 0, 0, 0]
			# build rs-mer dict
			original_obj.build_regular_syncmer_list(k_value=k)
			original_obj.build_rand_sync_based_on_syncmer_lists(sync_list=original_obj.regular_syncmer_list,
			                                                    loca_list=original_obj.regular_syncmer_loc_list)
			
			for mut_obj in mut_obj_list:
				# build rs-mer dict
				mut_obj.build_regular_syncmer_list(k_value=k)
				mut_obj.build_rand_sync_based_on_syncmer_lists(sync_list=mut_obj.regular_syncmer_list,
							                                    loca_list=mut_obj.regular_syncmer_loc_list)
				
				# get overlaps
				overlap_list = list(original_obj.rs_syncmer_dict.keys() & mut_obj.rs_syncmer_dict.keys())

				eval_matrix = get_rs_ext_sync_eval_matrix(input_object=original_obj, overlap_rs_mers_list=overlap_list)
				
				out_eval_list = [sum(x) for x in zip(out_eval_list, eval_matrix)]
			
			# output for this k value
			out_eval_list = [x / sample_size for x in out_eval_list]
			out_dict[str(k)] = out_eval_list
		
		# save to file
		out_df = pd.DataFrame(out_dict).transpose()
		out_df.columns = ['m', 'sc', 'mc', 'E']
		filename = "_".join(["regular_rs_mer_mut_rate", str(mutation_rate_p), ".csv"])
		out_df.to_csv(filename, header=True)
		
		
		### compare ext rs-mer
		out_dict = {}
		
		for k in k_list:
			out_eval_list = [0, 0, 0, 0]
			# build rs-mer dict
			original_obj.build_ext_syncmer_list(ext_ksize=k)
			original_obj.build_rand_sync_based_on_syncmer_lists(sync_list=original_obj.ext_syncmer_list,
			                                                    loca_list=original_obj.ext_syncmer_loc_list)
			
			for mut_obj in mut_obj_list:
				# build rs-mer dict
				mut_obj.build_ext_syncmer_list(ext_ksize=k)
				mut_obj.build_rand_sync_based_on_syncmer_lists(sync_list=mut_obj.ext_syncmer_list,
				                                               loca_list=mut_obj.ext_syncmer_loc_list)
				
				# get overlaps
				overlap_list = list(original_obj.rs_syncmer_dict.keys() & mut_obj.rs_syncmer_dict.keys())
				
				eval_matrix = get_rs_ext_sync_eval_matrix(input_object=original_obj, overlap_rs_mers_list=overlap_list)
				
				out_eval_list = [sum(x) for x in zip(out_eval_list, eval_matrix)]
			
			# output for this k value
			out_eval_list = [x / sample_size for x in out_eval_list]
			out_dict[str(k)] = out_eval_list
		
		# save to file
		out_df = pd.DataFrame(out_dict).transpose()
		out_df.columns = ['m', 'sc', 'mc', 'E']
		filename = "_".join(["ext_rs_mer_mut_rate", str(mutation_rate_p), ".csv"])
		out_df.to_csv(filename, header=True)

	
	
	### benchmark 5: how k-length affect ext-rand-sync match status
	original_obj = rand_sync(seq_string=out_seq, k_max=k_max, k_min=8, submer_length_dif=4, wmin=5, wmax=11)
	
	# p:m, and p:sc pairs to plot
	plot_match_dict = {}
	plot_sc_dict = {}
	plot_index_name = list(range(8,31,2))
	
	
	for mutation_rate_p in p_values:
		print("Simuation with mutation rate %s" % mutation_rate_p)
		# get a list of sequences
		mut_obj_list = []
		for i in range(sample_size):
			temp_seq = generate_mutated_strings(out_seq, p=mutation_rate_p)
			mut_obj_list.append(
				rand_sync(seq_string=temp_seq, k_max=k_max, k_min=8, submer_length_dif=4, wmin=5, wmax=11))
		
		### use ext-rs-mer for demo
		out_dict = {}
		
		for k in range(8,31,2):
			out_eval_list = [0, 0, 0, 0]
			# build rs-mer dict
			original_obj.build_ext_syncmer_list(ext_ksize=k)
			original_obj.build_rand_sync_based_on_syncmer_lists(sync_list=original_obj.ext_syncmer_list,
			                                                    loca_list=original_obj.ext_syncmer_loc_list)
			
			for mut_obj in mut_obj_list:
				# build rs-mer dict
				mut_obj.build_ext_syncmer_list(ext_ksize=k)
				mut_obj.build_rand_sync_based_on_syncmer_lists(sync_list=mut_obj.ext_syncmer_list,
				                                               loca_list=mut_obj.ext_syncmer_loc_list)
				
				# get overlaps
				overlap_list = list(original_obj.rs_syncmer_dict.keys() & mut_obj.rs_syncmer_dict.keys())
				
				eval_matrix = get_rs_ext_sync_eval_matrix(input_object=original_obj, overlap_rs_mers_list=overlap_list)
				
				out_eval_list = [sum(x) for x in zip(out_eval_list, eval_matrix)]
			
			# output for this k value
			out_eval_list = [x / sample_size for x in out_eval_list]
			out_dict[str(k)] = out_eval_list
		
		# save to file
		out_df = pd.DataFrame(out_dict).transpose()
		out_df.columns = ['m', 'sc', 'mc', 'E']
		filename = "_".join(["ext_rs_mer_mut_rate", str(mutation_rate_p), ".csv"])
		out_df.to_csv(filename, header=True)
		
		plot_match_dict[str(mutation_rate_p)] = list(out_df['m'])
		plot_sc_dict[str(mutation_rate_p)] = list(out_df['sc'])
		
	# save merged dict
	df1 = pd.DataFrame(plot_match_dict, index=plot_index_name)
	df1.to_csv("match_status_by_mutp_and_single_syncmer_length_k.csv", header=True)
	df2 = pd.DataFrame(plot_sc_dict, index=plot_index_name)
	df2.to_csv("seq_coverage_by_mutp_and_single_syncmer_length_k.csv", header=True)
	
	# make plot
	def generate_multiple_line_plot(input_df, input_title, out_filename):
		"""
		Generate multiple line plots for m or sc VS mutp and syncmer length k (full k-mer = 2*k)
		"""
		temp_df = input_df.copy()
		temp_df['sync_len'] = temp_df.index
		plot_df = temp_df.melt('sync_len', var_name="mut_rate", value_name="percent")
		
		fig, axs = plt.subplots(1, 1)
		sb.pointplot(x="sync_len", y="percent", hue='mut_rate', data=plot_df, ax=axs)
		axs.set_title(input_title)
		axs.set_ylim([0, 100])
		fig.savefig(out_filename, dpi=300)
		plt.close(fig)
	
	
	generate_multiple_line_plot(input_df=df1, input_title="Matched k-mers by mut_rate and syncmer length", out_filename="compare_m.png")
	generate_multiple_line_plot(input_df=df2, input_title="Seq coverage by mut_rate and syncmer length", out_filename="compare_sc.png")

	
	
	