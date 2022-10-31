# Helper class for creating extension-based strobemers and syncmers

import os
from Bio import SeqIO
import re
import numpy as np
import random
import khmer
import pandas as pd
import copy
import marisa_trie as mt
import pickle

# Basic functions
def generate_string_dict_from_genome(genome_file: str):
	"""
	Generate a dict of seq_id:sequence pairs from a fasta file
	
	:param genome_file: a fasta file for input genome
	:return: a dict of seq_id:sequence pairs
	"""
	seq_id = []
	seq_list = []
	fasta_records = SeqIO.parse(genome_file, 'fasta')
	for fasta in fasta_records:
		seq_list.append(str(fasta.seq).upper())
		seq_id.append(fasta.id.split(' ')[0])
	
	if len(seq_list) == 0:
		raise ValueError('No sequences')
	
	out_dict = dict(zip(seq_id, seq_list))
	return out_dict

def mutate_strings_from_a_string_dict(input_dict: dict, p: float):
	"""
	Mutate a string dict with mutation rate p
	
	:param input_dict: a string dict generate from previous helper function
	:param p: mutation rate p in [0, 1]. IID in every loci
	:return: a mutated string dict
	"""
	string_dict = copy.deepcopy(input_dict)
	
	for seq_id in string_dict.keys():
		input_string = string_dict[seq_id]
		
		# simple mutation event by Unif(0,1)
		string_length = len(input_string)
		random_unif1 = list(np.random.uniform(size=string_length))
		mutation_id = [x < p for x in random_unif1]  # use collections.Counter(mutation_id) to check number of elements
		
		# generate mutate string
		out_string = ""
		counter = 0
		while counter < string_length:
			if not mutation_id[counter]:
				# no mutation
				out_string += input_string[counter]
			else:
				# choose from insertion / deletion / substitution of 1/3 prob
				mut_type = random.choice(['ins', 'del', 'sub'])
				if mut_type == 'ins':
					# inserted letter
					out_string += random.choice('ACGT')
					# original letter
					out_string += input_string[counter]
				elif mut_type == 'del':
					# just don't add anything
					continue
				else:
					# substitution
					dest_space = "ACGT".replace(input_string[counter], "")
					out_string += random.choice(dest_space)
			counter += 1
		
		# replace the string
		string_dict[seq_id] = out_string
		
		return string_dict

def load_pkl_obj(pkl_file):
	"""
	A shortcut to load pickle file storing "ext_strobemer_obj"
	"""
	with open(pkl_file, "rb") as fp:
		out_obj = pickle.load(fp)
	return out_obj


# create the main class for extension-based strobemers
class ext_strobemer_obj:
	def __init__(self, k: int, n: int, l: int, smin: int=5, wmin: int=25, wmax: int=50,
	             syncmer_lenth: int=20, sync_submer_lenth: int=5,
	             prime: int=6229, label: str="no_label", filename: str="no_filename"):
		"""
		Init a ext_strobemer_obj
		
		:param k: length of full k-mer (k = n*l), it's also the max k value for extension method
		:param n: number of strobes
		:param l: length of single strobe
		:param smin: minimal strobe length for the selection of strobe to build extension-based strobemers
		:param wmin: min window size for strobe selection
		:param wmax: max window size for strobe selection
		:param syncmer_lenth: length of the syncmer
		:param sync_submer_lenth: length of substr s to select syncmer from long sequences
		:param prime: a prime number for strobe selection (don't need a large one given that the window is small)
		:param label: keyword for input file / output file
		:param filename: store genome file if necessary
		"""
		# load parameters
		self.ksize = k
		self.smin = smin
		if smin > l:
			raise Exception("smin should NOT be larger than strobe size")
		self.num_strobe = n
		self.len_strobe = l
		if k != n * l:
			raise Exception("Strobemer length NOT EQUAL to k-mer length, plaese double check parameters")
		self.wmin = wmin
		self.wmax = wmax
		self.syncmer_length = syncmer_lenth
		self.sync_submer_lenth = sync_submer_lenth
		self.prime = prime
		self.label = label
		self.filename = filename
		
		# init empty properties
		self.seq_dict = {}
		self.kmer_dict = {}
		self.kmer_tst = None # for TST object
		self.regu_randstrobe_dict = {}
		self.ext_randstrobe_dict = {}
		self.ext_randstrobe_tst = None # for TST object
		self.syncmer_dict = {}
		self.syncmer_tst = None # for TST object
		self.aux_syncmer = {}
	
	def print_info(self):
		"""
		Print out the information for this object
		"""
		print("The information of this ext_randstrobe object: \n\n")
		print("K-mer length: %s" % self.ksize)
		print("Strobe number: %s" % self.num_strobe)
		print("Strobe length: %s" % self.len_strobe)
		print("Extension-strobe length: %s" % self.smin)
		print("Min window size: %s" % self.wmin)
		print("Max window size: %s" % self.wmax)
		print("Suncymer length is %s" % self.syncmer_length)
		print("Syncmer substr length is %s" % self.sync_submer_lenth)
		print("Prime number for randstrobe: %s" % self.prime)
		print("File label: %s" % self.label)
		print("File name: %s \n" % self.filename)
		
		# also check the dict inside if any
		if self.seq_dict:
			print("Sequences have been loaded with name and length:")
			print(list(self.seq_dict.keys()))
			print([len(x) for x in list(self.seq_dict.values())])
			print("")
		
		if self.kmer_dict and self.kmer_tst:
			print("K-mer dict: %s unique records" % len(self.kmer_dict))
			print("In addition, a KTST has been built for prefix lookup \n")
			
		if self.regu_randstrobe_dict:
			print("Regular randstrobe dict: %s unique records \n" % len(self.regu_randstrobe_dict))
			
		if self.ext_randstrobe_dict and self.ext_randstrobe_tst:
			print("Extension-based randstrobe dict: %s unique records" % len(self.ext_randstrobe_dict))
			print("In addition, a TST has been built for ext_randstrobes \n")
			
		if self.syncmer_dict and self.syncmer_tst:
			print("Syncmer dict: %s unique records" %len(self.syncmer_dict))
			print("In addition, a TST has been built for syncmers")
	
	def load_seq_dict(self, seq_dict):
		"""
		Load a seq dict to this object
		"""
		self.seq_dict = copy.deepcopy(seq_dict)
	
	@staticmethod
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
			
	def build_kmer_dict_with_start_pos(self):
		"""
		Build a regular k-mer dict for this object
		"""
		kmer_dict = self.kmer_dict
		ksize = self.ksize
		
		for seq_id in self.seq_dict:
			seq = self.seq_dict[seq_id]
			for i in range(len(seq) - ksize + 1):
				kmer = seq[i:i + ksize]
				self.write_or_extend_kmer_record_to_dict(kmer_dict, kmer, item_value="_".join([seq_id, "i", str(i)]))
				
		# construct KTST for kmer_dict
		self.kmer_tst = mt.Trie(kmer_dict.keys())
				
	def build_both_randstrobe_dict(self):
		"""
		Build a regular randstrobe dict AND a ext-based randstrobe for this object
		"""
		randstrobe_dict = self.regu_randstrobe_dict
		ext_rs_dict = self.ext_randstrobe_dict
		n = self.num_strobe
		l = self.len_strobe
		wmin = self.wmin
		wmax = self.wmax
		prime = self.prime
		seq_dict = self.seq_dict
		smin = self.smin
		
		for seq_id in seq_dict:
			seq = seq_dict[seq_id]
			for i in range(len(seq) - n * l + 1):
				# max window size: only activate when approaching the tail of seq, o.w. is w_max
				wu = min(wmax, int((len(seq) - i - l + 1) / (n - 1)) )  # in case of wu -wl < 1, use int
				# min window size: make sure no overlap strobes: lenth l < window min
				wl = max(wmin - (wmax - wu), l)
				# special cases near the end of sequences (may generate empty records in the j loop)
				if wu < wl:
					break
				elif wu == wl:
					# start from here, all remaining strobemers are regular kmer (may have wu - wl < 1 -> several kmers at end)
					fix_seq = seq[i:i + n * l]
					fix_pos = [str(i)]
					for j in range(1, n):
						fix_pos.append(str(i + j * l))
					# add to randstrobe_dict
					fix_pos.insert(0, seq_id+"_i")
					self.write_or_extend_kmer_record_to_dict(randstrobe_dict, fix_seq, item_value="_".join(fix_pos))
					continue
					
				# ranstrobe: 1st strobe
				rand_seq = last_strobe = seq[i:i + l]  # string instance is immutable
				rand_pos = [str(i)]
				# randstrobe loop, range(2,n+1) is for window size purpose, while in L182 we don't need it (so use 1,n)
				for j in range(2, n + 1):
					# this will generate empty list when wu <= wl
					select_window = list(range(i + wl + (j - 2) * wu, i + (j - 1) * wu, 1))
					current_pos = 0
					current_hash = prime + 10
					# find argmin in current window
					for pos in select_window:
						new_hash = (khmer.hash_no_rc_murmur3(seq[pos:pos + l]) + khmer.hash_no_rc_murmur3(last_strobe)) % prime
						if new_hash < current_hash:
							current_hash = new_hash
							current_pos = pos
					# update last strobe
					last_strobe = seq[current_pos:current_pos + l]
					# record
					rand_pos.append(str(current_pos))
					rand_seq += seq[current_pos:current_pos + l]
				# add to the dict
				rand_pos.insert(0, seq_id+"_i")
				self.write_or_extend_kmer_record_to_dict(randstrobe_dict, rand_seq, item_value="_".join(rand_pos))
				del rand_seq, rand_pos, last_strobe
				
				# ext_randstrobe: 1st strobe same
				rand_seq = last_strobe = seq[i:i + l]
				rand_pos = [str(i)]
				# ext_rs loop
				for j in range(2, n + 1):
					# this will generate empty list when wu <= wl
					select_window = list(range(i + wl + (j - 2) * wu, i + (j - 1) * wu, 1))
					current_pos = 0
					current_hash = prime + 10
					# find argmin via shorter strobe of length smin
					for pos in select_window:
						new_hash = (khmer.hash_no_rc_murmur3(seq[pos:pos + smin]) + khmer.hash_no_rc_murmur3(last_strobe[:smin])) % prime
						if new_hash < current_hash:
							current_hash = new_hash
							current_pos = pos
					# update last strobe
					last_strobe = seq[current_pos:current_pos + l]
					# record
					rand_pos.append(str(current_pos))
					rand_seq += seq[current_pos:current_pos + l]
				# add to the dict
				rand_pos.insert(0, seq_id + "_i")
				self.write_or_extend_kmer_record_to_dict(ext_rs_dict, rand_seq, item_value="_".join(rand_pos))
				
		# also build a TST for ext_strobemer
		self.ext_randstrobe_tst = mt.Trie(ext_rs_dict.keys())
	
	@staticmethod
	def is_kmer_an_open_syncmer(input_string: str, submer_length: int):
		"""
		Determine if a given kmer (string) is an OPEN syncmer at given submer_length via lexicographic order.
		:param input_string: a kmer
		:param submer_length: the length substring in this kmer for syncmer selection
		"""
		temp_start = input_string[:submer_length]
		for i in range(len(input_string) - submer_length + 1):
			if temp_start > input_string[i:i+submer_length]:
				return False
		# o.w., the start is minimum, so this is an open syncmer
		return True
	
	def build_open_syncmer(self):
		"""
		Build a syncmer dict for this object, and also maintain a auxiliary dict to store chr:binary list of syncmer start location.
		This auxiliary dictionary is for faster randstrobe-syncmer index.
		"""
		sync_len = self.syncmer_length
		sub_len = self.sync_submer_lenth
		sync_dict = self.syncmer_dict
		aux_infor = self.aux_syncmer
		
		for seq_id in self.seq_dict:
			seq = self.seq_dict[seq_id]
			# use a binary list to save the syncmer start location
			sync_binary_record = [0] * len(seq)
			# get substring, i.e. k-mers
			for i in range(len(seq) - sync_len + 1):
				kmer = seq[i:i + sync_len]
				if self.is_kmer_an_open_syncmer(input_string=kmer, submer_length=sub_len):
					sync_binary_record[i] = 1
					self.write_or_extend_kmer_record_to_dict(sync_dict, kmer, item_value="_".join([seq_id, "i", str(i)]))
			# update auxiliary information (this is int list and will be re-initiated in next round loop)
			aux_infor[seq_id] = sync_binary_record
			
		# make TST
		self.syncmer_tst = mt.Trie(sync_dict.keys())
		
	def export_to_pkl(self):
		"""
        Explort to a pkl file.
        """
		export_file_name = "".join(['ext_strobemer_obj_', self.label, '.pkl'])
		with open(export_file_name, 'wb') as outp:
			pickle.dump(self, outp, pickle.HIGHEST_PROTOCOL)
			
	def kmer_query(self, kmer: str, kmer_type: str):
		"""
		Query (exact match) if a given kmer is in a specific kmer dict
		
		:param kmer: query kmer to search for exact match
		:param kmer_type: which kmer dict to search, including ['kmer', 'rs', 'ers', 'sync']
		:return: genomic location of the kmer if found, o.w. return empty list
		"""
		if kmer_type == 'kmer':
			try:
				return self.kmer_dict[kmer]
			except KeyError:
				return []
		elif kmer_type == 'rs':
			try:
				return self.regu_randstrobe_dict[kmer]
			except KeyError:
				return []
		elif kmer_type == 'ers':
			try:
				return self.ext_randstrobe_dict[kmer]
			except KeyError:
				return []
		elif kmer_type == 'sync':
			try:
				return self.syncmer_dict[kmer]
			except KeyError:
				return []
			
	def prefix_query(self, kmer: str, kmer_type: str):
		"""
		Find the k-mer that contains a given prefix from the pre-built TST
		
		:param kmer: query prefix
		:param kmer_type: which TST to search, including [kmer, ers, sync]
		:return: a list of kmers containing the given prefix, or empty list if not found
		"""
		if kmer_type == 'kmer':
			return self.kmer_tst.keys(kmer)
		elif kmer_type == 'ers':
			return self.ext_randstrobe_tst.keys(kmer)
		elif kmer_type == 'sync':
			return self.syncmer_tst.keys(kmer)

		
		
		
		
