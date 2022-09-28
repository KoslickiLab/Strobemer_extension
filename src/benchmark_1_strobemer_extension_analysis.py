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


### basic functions
def generate_string_list_from_genome(input_genome):
	"""
	Transfer input genome to a string
	"""
	# store fasta to a list (in case multiple records line)
	fasta_sequences = SeqIO.parse(input_genome, 'fasta')
	seq_list = []
	for fasta in fasta_sequences:
		seq_list.append(str(fasta.seq))
	# merge every record to a string
	if len(seq_list) == 0:
		raise Exception("couldn't read fasta from input file, please double check!")
	else:
		print("There are " + str(len(seq_list)) + " record(s).")
		return seq_list

def prepare_input_string(genome_file, seq_length):
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
	
	return out_seq

def generate_mutated_strings(input_string, p):
	"""
	Generate a mutated string by simple mutation model
	"""
	# simple mutation event by Unif(0,1)
	string_length = len(input_string)
	# determine the location to be mutated
	random_unif1 = list(np.random.uniform(size=string_length))
	mutation_id = [x < p for x in random_unif1]  # use collections.Counter(mutation_id) to check number of elements
	# mutate string
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
				out_string += "-"
			else:
				dest_space = "ACGT".replace(input_string[counter], "")
				out_string += random.choice(dest_space)
		counter += 1
	# return
	out_string = out_string.replace("-", "")  # trail the deletion letter
	return out_string

def build_kmer_dict_with_start_pos(seq, ksize: int):
	"""
	Build a dict with key=kmer, value=kmer start pos (include multiple hits)
	:param seq:
	:param k:
	:return: a dict
	"""
	temp = dict()
	for i in range(len(seq) - ksize + 1):
		kmer = seq[i:i + ksize]
		if kmer in temp:
			temp[kmer].append(i)
		else:
			temp[kmer] = [i]
	return temp

def dup_count_for_kmer_dict(input_dict):
	"""
	Count dup kmers in a given kmer dict: kmer:[pos]
	"""
	for item in ['kmer', 'spaced_sparse', 'spaced_dense', 'minstrobe', 'randstrobe', 'hybridstrobe']:
		temp_dict = input_dict[item]
		all_value_len = [len(x) for x in temp_dict.values()]
		multi_hit = [x > 1 for x in all_value_len]
		print("Dup count of %s \nTotal distinct: %s\nDuplicated items: %s \n\n" % (
			item, len(all_value_len), sum(multi_hit)))

def count_island_from_binary_vector(input_array):
	"""
	Count the island size from a binary vector in kmer match coverage
	:param input_array:
	:return: a dict with pos:length pairs for islands
	"""
	zero_dict = dict()
	zero_dict['start'] = []
	zero_dict['length'] = []
	i = 0
	n = len(input_array)
	input_array = list(input_array)
	# append 1 to terminate the loop without error
	input_array.append(1)
	count = 0
	
	while i < n:
		if input_array[i] != 0:
			i += 1
		else:
			zero_dict['start'].append(i)
			# start to count length
			while input_array[i] == 0:
				count += 1
				i += 1
			zero_dict['length'].append(count)
			# reset count
			count = 0
	
	# calculate E size
	temp = np.array(zero_dict['length'])
	out_E = np.sum(temp ** 2) / n
	
	return out_E

def truncate_kmer_dict(kmer_dict, kmer_truncate, single_strobe_truncate):
	"""
	Truncate strobmers from a dict object
	:return:
	"""
	dict_input = copy.deepcopy(kmer_dict)  # use deepcopy to avoid unwanted change in input df
	
	if kmer_truncate != single_strobe_truncate * dict_input['strobemer_label'][0]:
		raise Exception("K-mer truncation length doesn't match strobemer, please double check!!!")
	
	out_dict = dict()
	out_dict['length_str'] = dict_input['length_str']
	
	# new strobe size
	temp_label = dict_input['strobemer_label']
	old_strobe_size = temp_label[1]
	new_strobe_size = temp_label[1] - single_strobe_truncate
	temp_label[1] = new_strobe_size
	out_dict['strobemer_label'] = temp_label
	
	# new ksize
	new_ksize = dict_input['ksize'] - kmer_truncate
	out_dict['ksize'] = new_ksize
	
	# k-mers can be directly truncated
	for item in ['kmer', 'spaced_sparse', 'spaced_dense']:
		temp_dict = dict()
		for single_kmer in dict_input[item]:
			truncated_kmer = single_kmer[:new_ksize]
			if truncated_kmer in temp_dict:  # dup found
				temp_dict[truncated_kmer].extend(dict_input[item][single_kmer])
			else:
				temp_dict[truncated_kmer] = dict_input[item][single_kmer]
		# add to output
		out_dict[item] = temp_dict
		del temp_dict
	
	# single strobe level truncation
	for item in ['minstrobe', 'randstrobe', 'hybridstrobe']:
		temp_dict = dict()
		for single_kmer in dict_input[item]:
			truncated_strobe = []
			for i in range(temp_label[0]):
				truncated_strobe.append(single_kmer[i * old_strobe_size: i * old_strobe_size + new_strobe_size])
			truncated_kmer = "".join(truncated_strobe)
			if truncated_kmer in temp_dict:  # dup found
				temp_dict[truncated_kmer].extend(dict_input[item][single_kmer])
			else:
				temp_dict[truncated_kmer] = dict_input[item][single_kmer]
		# add to output
		out_dict[item] = temp_dict
		del temp_dict
	
	return out_dict


### strobemer function
def construct_strobemer(seq, n, l, wmin, wmax, prime=6229, x=3):
	"""
	Generate minstrobe for given seq
	:param seq:
	:param l: strobe length
	:param n: strobe number
	:param wmin: window gap min
	:param wmax: window gap max
	:param prime: a prime # to mod, don't need a large one
	#:param: x: partition parameter, use 3
	:return: dict with strobemer:[ [m1...mn], ...  ]
	"""
	minstrobe_dict = dict()
	randstrobe_dict = dict()
	hybridstrobe_dict = dict()
	# loop through index
	for i in range(len(seq) - n * l + 1):
		# max window size: only activate when approaching the tail of seq, o.w. is w_max
		wu = int(min(wmax, (len(seq) - i - l + 1) / (n - 1)))  # in case of wu -wl < 1, use int
		# min window size: make sure no overlap strobes: lenth l < window min
		wl = max(wmin - (wmax - wu), l)
		# manually stop (may generate empty records in the j loop)
		if wu < wl:
			break
		elif wu == wl:
			# start from here, all remaining strobemers are regular kmer (may have wu - wl < 1 -> several kmers at end)
			fix_seq = seq[i:i + n * l]
			fix_pos = [i]
			for j in range(1, n):
				fix_pos.append(i + j * l)
			# add to min_dict
			if fix_seq in minstrobe_dict:
				minstrobe_dict[fix_seq].append(fix_pos)
			else:
				minstrobe_dict[fix_seq] = [fix_pos]
			# also add to rand_dict
			if fix_seq in randstrobe_dict:
				randstrobe_dict[fix_seq].append(fix_pos)
			else:
				randstrobe_dict[fix_seq] = [fix_pos]
			# end this loop here
			continue
		
		# minstrobe: 1st strobe
		m1 = seq[i:i + l]
		min_seq = m1
		min_pos = [i]
		# minstrobe loop
		for j in range(2, n + 1):
			# this will generate empty list when wu <= wl
			select_window = list(range(i + wl + (j - 2) * wu, i + (j - 1) * wu, 1))
			current_pos = 0
			current_hash = prime + 10  # init loop
			# find argmin in current window
			for pos in select_window:
				new_hash = khmer.hash_no_rc_murmur3(seq[pos:pos + l]) % prime
				if new_hash < current_hash:
					current_hash = new_hash
					current_pos = pos
			# record
			min_pos.append(current_pos)
			min_seq += seq[current_pos:current_pos + l]
		# add to the dict
		if min_seq in minstrobe_dict:
			# add strobe starting indices [x,x,x] into a list, so use append
			minstrobe_dict[min_seq].append(min_pos)
		else:
			minstrobe_dict[min_seq] = [min_pos]
		
		# ranstrobe: 1st strobe
		rand_seq = m1
		rand_pos = [i]
		last_strobe = m1
		# randstrobe loop
		for j in range(2, n + 1):
			# this will generate empty list when wu <= wl
			select_window = list(range(i + wl + (j - 2) * wu, i + (j - 1) * wu, 1))
			current_pos = 0
			current_hash = prime + 10  # init loop
			# find argmin in current window
			for pos in select_window:
				new_hash = (khmer.hash_no_rc_murmur3(seq[pos:pos + l]) + khmer.hash_no_rc_murmur3(last_strobe)) % prime
				if new_hash < current_hash:
					current_hash = new_hash
					current_pos = pos
			# update last strobe
			last_strobe = seq[current_pos:current_pos + l]
			# record
			rand_pos.append(current_pos)
			rand_seq += seq[current_pos:current_pos + l]
		# add to the dict
		if rand_seq in randstrobe_dict:
			randstrobe_dict[rand_seq].append(rand_pos)
		else:
			randstrobe_dict[rand_seq] = [rand_pos]
		
		# hybridstrobe: 1st strobe
		hy_seq = m1
		hy_pos = [i]
		last_strobe = m1
		# hybridstrobe loop
		for j in range(2, n + 1):
			# this will generate empty list when wu <= wl
			select_window = list(range(i + wl + (j - 2) * wu, i + (j - 1) * wu, 1))
			# additional step to subset selection window by x
			slice_index = khmer.hash_no_rc_murmur3(last_strobe) % x
			if len(select_window) < x:
				# no need to slice with few elements
				sub_window = select_window
			else:
				sub_window = list(np.array_split(select_window, x)[slice_index])
			# now we get the r-th subwindow to play
			current_pos = 0
			current_hash = prime + 10  # init loop
			# find argmin in current window
			for pos in sub_window:
				new_hash = (khmer.hash_no_rc_murmur3(seq[pos:pos + l]) + khmer.hash_no_rc_murmur3(last_strobe)) % prime
				if new_hash < current_hash:
					current_hash = new_hash
					current_pos = pos
			# update last strobe
			last_strobe = seq[current_pos:current_pos + l]
			# record
			hy_pos.append(current_pos)
			hy_seq += seq[current_pos:current_pos + l]
		# add to the dict
		if hy_seq in hybridstrobe_dict:
			hybridstrobe_dict[hy_seq].append(hy_pos)
		else:
			hybridstrobe_dict[hy_seq] = [hy_pos]
	
	# return results
	return minstrobe_dict, randstrobe_dict, hybridstrobe_dict

def construct_strobemer_from_kmin(seq, kmin, n, l, wmin, wmax, prime=6229, x=3):
	"""
	Generate minstrobe for given seq
	:param seq:
	:param kmin: shortest substring for strobe selection
	:param l: strobe length
	:param n: strobe number
	:param wmin: window gap min
	:param wmax: window gap max
	:param prime: a prime # to mod, don't need a large one
	#:param: x: partition parameter, use 3
	:return: dict with strobemer:[ [m1...mn], ...  ]
	"""
	minstrobe_dict = dict()
	randstrobe_dict = dict()
	hybridstrobe_dict = dict()
	# loop through index
	for i in range(len(seq) - n * l + 1):
		# max window size: only activate when approaching the tail of seq, o.w. is w_max
		wu = int(min(wmax, (len(seq) - i - l + 1) / (n - 1)))  # in case of wu -wl < 1, use int
		# min window size: make sure no overlap strobes: lenth l < window min
		wl = max(wmin - (wmax - wu), l)
		# manually stop (may generate empty records in the j loop)
		if wu < wl:
			break
		elif wu == wl:
			# start from here, all remaining strobemers are regular kmer (may have wu - wl < 1 -> several kmers at end)
			fix_seq = seq[i:i + n * l]
			fix_pos = [i]
			for j in range(1, n):
				fix_pos.append(i + j * l)
			# add to min_dict
			if fix_seq in minstrobe_dict:
				minstrobe_dict[fix_seq].append(fix_pos)
			else:
				minstrobe_dict[fix_seq] = [fix_pos]
			# also add to rand_dict
			if fix_seq in randstrobe_dict:
				randstrobe_dict[fix_seq].append(fix_pos)
			else:
				randstrobe_dict[fix_seq] = [fix_pos]
			# end this loop here
			continue
		
		# minstrobe: 1st strobe
		m1 = seq[i:i + l]
		min_seq = m1
		min_pos = [i]
		# minstrobe loop
		for j in range(2, n + 1):
			# this will generate empty list when wu <= wl
			select_window = list(range(i + wl + (j - 2) * wu, i + (j - 1) * wu, 1))
			current_pos = 0
			current_hash = prime + 10  # init loop
			# find argmin in current window for kmin-mer and extend to length l !!!
			for pos in select_window:
				new_hash = khmer.hash_no_rc_murmur3(seq[pos:pos + kmin]) % prime
				if new_hash < current_hash:
					current_hash = new_hash
					current_pos = pos
			# record
			min_pos.append(current_pos)
			min_seq += seq[current_pos:current_pos + l]
		# add to the dict
		if min_seq in minstrobe_dict:
			# add strobe starting indices [x,x,x] into a list, so use append
			minstrobe_dict[min_seq].append(min_pos)
		else:
			minstrobe_dict[min_seq] = [min_pos]
		
		# ranstrobe: 1st strobe
		rand_seq = m1
		rand_pos = [i]
		last_strobe = m1
		# randstrobe loop
		for j in range(2, n + 1):
			# this will generate empty list when wu <= wl
			select_window = list(range(i + wl + (j - 2) * wu, i + (j - 1) * wu, 1))
			current_pos = 0
			current_hash = prime + 10  # init loop
			# find argmin in current window
			for pos in select_window:
				new_hash = (khmer.hash_no_rc_murmur3(seq[pos:pos + kmin]) + khmer.hash_no_rc_murmur3(
					last_strobe[:kmin])) % prime
				if new_hash < current_hash:
					current_hash = new_hash
					current_pos = pos
			# update last strobe
			last_strobe = seq[current_pos:current_pos + l]
			# record
			rand_pos.append(current_pos)
			rand_seq += seq[current_pos:current_pos + l]
		# add to the dict
		if rand_seq in randstrobe_dict:
			randstrobe_dict[rand_seq].append(rand_pos)
		else:
			randstrobe_dict[rand_seq] = [rand_pos]
		
		# hybridstrobe: 1st strobe
		hy_seq = m1
		hy_pos = [i]
		last_strobe = m1
		# hybridstrobe loop
		for j in range(2, n + 1):
			# this will generate empty list when wu <= wl
			select_window = list(range(i + wl + (j - 2) * wu, i + (j - 1) * wu, 1))
			# additional step to subset selection window by x
			slice_index = khmer.hash_no_rc_murmur3(last_strobe[:kmin]) % x
			if len(select_window) < x:
				# no need to slice with few elements
				sub_window = select_window
			else:
				sub_window = list(np.array_split(select_window, x)[slice_index])
			# now we get the r-th subwindow to play
			current_pos = 0
			current_hash = prime + 10  # init loop
			# find argmin in current window
			for pos in sub_window:
				new_hash = (khmer.hash_no_rc_murmur3(seq[pos:pos + kmin]) + khmer.hash_no_rc_murmur3(
					last_strobe[:kmin])) % prime
				if new_hash < current_hash:
					current_hash = new_hash
					current_pos = pos
			# update last strobe
			last_strobe = seq[current_pos:current_pos + l]
			# record
			hy_pos.append(current_pos)
			hy_seq += seq[current_pos:current_pos + l]
		# add to the dict
		if hy_seq in hybridstrobe_dict:
			hybridstrobe_dict[hy_seq].append(hy_pos)
		else:
			hybridstrobe_dict[hy_seq] = [hy_pos]
	
	# return results
	return minstrobe_dict, randstrobe_dict, hybridstrobe_dict

def unit_test_2_ways_strobemer(test_seq):
	# when there is no extension: 2 ways get exact same results
	(min1, rand1, hy1) = construct_strobemer(test_seq, n=2, l=10, wmin=25, wmax=50)
	# kmin = 10 = l, this is the marginal case where extension strobemer is same to regular one
	(min2, rand2, hy2) = construct_strobemer_from_kmin(test_seq, kmin=10, n=2, l=10, wmin=25, wmax=50)
	assert min1 == min2
	assert rand1 == rand2
	assert hy1 == hy2
	
	# regular l=15, extension from 10 to 15: totally different kmer sets will be obtained
	(min1, rand1, hy1) = construct_strobemer(test_seq, n=2, l=15, wmin=25, wmax=50)
	(min2, rand2, hy2) = construct_strobemer_from_kmin(test_seq, kmin=10, n=2, l=15, wmin=25, wmax=50)
	assert min1 != min2
	# this extension one will share 10-mer prefix with regular strobemer at l=10
	(min3, rand3, hy3) = construct_strobemer(test_seq, n=2, l=10, wmin=25, wmax=50)
	
	count = 0
	for kmers in min2.keys():
		temp_10_prefix = "".join([kmers[:10], kmers[15:25]])
		if temp_10_prefix not in min3:
			count += 1
	print("Number of unmatched keys for MinStrobe is %s" %count)
	
	count = 0
	for kmers in rand2.keys():
		temp_10_prefix = "".join([kmers[:10], kmers[15:25]])
		if temp_10_prefix not in rand3:
			count += 1
	print("Number of unmatched keys for RandStrobe is %s" % count)
	
	count = 0
	for kmers in hy2.keys():
		temp_10_prefix = "".join([kmers[:10], kmers[15:25]])
		if temp_10_prefix not in hy3:
			count += 1
	print("Number of unmatched keys for RandStrobe is %s" % count)


### workflow functions
def s1_generate_all_types_of_kmers(input_string, kmin=5, k_value=30, strobe_size=10, wmin=25, wmax=50):
	"""
	Generate regular k-mer, spaced k-mer, and 3 types of strobmer
	:return: a dict with keys ['kmer', 'spaced_sparse', 'spaced_dense', 'minstrobe', 'randstrobe', 'hybridstrobe'], values are dicts of kmers with starting point
	e.g. dict['kmer'] stores "AAA":[1,3,5] (if dup 3 times, o.w. is [1]; strobmers are [ [1,5,9] ] for 3 strobes
	"""
	out_dict = dict()
	out_dict['ksize'] = k_value
	out_dict['length_str'] = len(input_string)
	### regular kmer
	out_dict['kmer'] = build_kmer_dict_with_start_pos(input_string, k_value)
	
	### spaced_kmer sparse: L=3k, 1/3 occupied
	temp_dict = build_kmer_dict_with_start_pos(input_string, k_value * 3)
	# this is a 3k-mer, then we modify it to a kmer with fixed pattern (uniform)
	for k3mer in list(temp_dict):
		# spaced_kmer sparse, keep every 3rd letter
		new_kmer = k3mer[::3]
		# refresh keys, this may introduce dup! So need to use "EXTEND", but in regular kmer we don't need it
		if new_kmer in temp_dict:
			temp_dict[new_kmer].extend(temp_dict.pop(k3mer))
		else:
			temp_dict[new_kmer] = temp_dict.pop(k3mer)
	out_dict['spaced_sparse'] = temp_dict
	del temp_dict, k3mer
	
	### spaced_kmer dense: L=1.5k, 2/3 occupied
	temp_dict = build_kmer_dict_with_start_pos(input_string, int(k_value * 1.5))
	for k1_5mer in list(temp_dict):
		# spaced_kmer dense
		kmer_2_list = list(k1_5mer)
		del kmer_2_list[2::3]  # delete every 3rd letter
		new_kmer = "".join(kmer_2_list)
		# refresh keys, this may introduce dup!
		if new_kmer in temp_dict:
			temp_dict[new_kmer].extend(temp_dict.pop(k1_5mer))
		else:
			temp_dict[new_kmer] = temp_dict.pop(k1_5mer)
	out_dict['spaced_dense'] = temp_dict
	del temp_dict, k1_5mer
	
	strobe_num = k_value // strobe_size
	if k_value / strobe_size != int(k_value / strobe_size):
		warnings.warn(
			"Notice: k_value can't be divided by strobe size, the length may NOT match between kmer and strobemer")
	
	### minstrobe, randstrobe and hybridstrobe
	out_dict['strobemer_label'] = [strobe_num, strobe_size, wmin, wmax]
	out_dict['minstrobe'], out_dict['randstrobe'], out_dict['hybridstrobe'] = \
		construct_strobemer_from_kmin(input_string, kmin=kmin, n=strobe_num, l=strobe_size, wmin=wmin, wmax=wmax, prime=6229, x=3)
	
	### output
	return out_dict


def s2_get_eval_matrix_between_kmer_dict_from_2_strings(dict_input, dict_mutate):
	"""
	Get py-version SIM-R table for 1 pair, may repeat for many times
	:param dict1: kmer dict for str1
	:param dict2: kmer dict for mutated str1
	:return: a list [category, parameter, m, sc, mc, E]
	"""
	out_dict = dict()  # store the 5*4 statis for each category
	ksize = dict_input['ksize']
	str_length = dict_input['length_str']
	strobe_size = dict_input['strobemer_label'][1]
	# kmers (start pos) and strobes ([start1, start2, ...] has different data structure, so need 2 loops
	
	for item in ['kmer', 'spaced_sparse', 'spaced_dense']:
		out_dict[item] = []
		### overlapped substrings
		overlap_list = list(dict_input[item].keys() & dict_mutate[item].keys())
		
		### m: % of matched substring
		out_dict[item].append(100 * len(overlap_list) / len(dict_input[item].keys()))
		
		### sc: seq coverage of input string
		sc_binary = np.array([0] * str_length)
		for shared in overlap_list:
			for start_pos in dict_input[item][shared]:  # this will count all dups, use ....[0] if only want 1st hit
				if item == 'kmer':
					sc_binary[start_pos:start_pos + ksize] = 1
				elif item == 'spaced_sparse':  # keep every 3rd letter, k3mer[::3]
					sc_binary[range(start_pos, start_pos + 3 * ksize, 3)] = 1
				else:  # dense: delete every 3rd letter, del kmer_2_list[2::3]
					# keep 1st and 2nd for every 3 letters in the 1.5*K region
					sc_binary[range(start_pos + 0, int(start_pos + 1.5 * ksize), 3)] = 1
					sc_binary[range(start_pos + 1, int(start_pos + 1.5 * ksize), 3)] = 1
		# coverage in original string
		out_dict[item].append(100 * sc_binary.sum() / str_length)
		
		### mc: match overage, subsring full cover in original seq
		sc_binary = np.array([0] * str_length)
		for shared in overlap_list:
			for start_pos in dict_input[item][shared]:  # this will count all dups, use ....[0] if only want 1st hit
				if item == 'kmer':  # kmer sc = mc
					sc_binary[start_pos:start_pos + ksize] = 1
				elif item == 'spaced_sparse':  # full 3k region
					sc_binary[start_pos: start_pos + 3 * ksize] = 1
				else:  # dense: 1.5k region
					sc_binary[start_pos: int(start_pos + 1.5 * ksize)] = 1
		# match coverage
		out_dict[item].append(100 * sc_binary.sum() / str_length)
		
		### E size
		e_size = count_island_from_binary_vector(sc_binary)
		out_dict[item].append(round(e_size))
	
	for item in ['minstrobe', 'randstrobe', 'hybridstrobe']:
		out_dict[item] = []
		### overlapped substrings
		overlap_list = list(dict_input[item].keys() & dict_mutate[item].keys())
		
		### m: % of matched substring
		out_dict[item].append(100 * len(overlap_list) / len(dict_input[item].keys()))
		
		### sc: seq coverage of input string
		sc_binary = np.array([0] * str_length)
		for shared in overlap_list:
			for strobe_start_pos in dict_input[item][
				shared]:  # this will count all dups, use ....[0] if only want 1st hit
				for single_strobe in strobe_start_pos:
					sc_binary[single_strobe: single_strobe + strobe_size] = 1
		# coverage in original string
		out_dict[item].append(100 * sc_binary.sum() / str_length)
		
		### mc: match coverage
		sc_binary = np.array([0] * str_length)
		for shared in overlap_list:
			# this will count all dups, use ....[0] if only want 1st hit
			for strobe_start_pos in dict_input[item][shared]:
				# we only need the 1st and last strobe pos for strobemer spanning
				sc_binary[strobe_start_pos[0]: strobe_start_pos[-1] + strobe_size] = 1
		# match coverage
		out_dict[item].append(100 * sc_binary.sum() / str_length)
		
		### E size will use the mc coverage
		e_size = count_island_from_binary_vector(sc_binary)
		out_dict[item].append(round(e_size))
	
	return out_dict


def s3_workflow_with_npk_from_input_dict_and_mutate_str(raw_dict: dict, seq_list: list, strobe_num: int, p: float,
                                                        k: int):
	"""
	Get eval matrix for fixed n,p,k based on a kmer dict (from raw str) and replicates of mutated strings
	:param raw_dict: a pre-built kmer dict from input string, can be a full kmer or TRUNCATED kmer
	:param seq_list: mutated string in a list
	:param strobe_num: number of strobes
	:param p: mutation rate
	:param k: k value
	:return: a dict with average values on each item
	"""
	out_list = []  # store those output dict to calculate mean
	
	# confirm strobe size is integer (to match the k)
	strobe_size = k / strobe_num
	if strobe_size != int(strobe_size):
		raise Exception("Please carefully pick k and strobe num s.t. there is no decimal!")
	
	# mutate str for n times
	for mutated_str in seq_list:
		dict_mutate = s1_generate_all_types_of_kmers(input_string=mutated_str, k_value=k, strobe_size=int(strobe_size),
		                                             wmin=wmin, wmax=wmax)
		temp_eval = s2_get_eval_matrix_between_kmer_dict_from_2_strings(raw_dict, dict_mutate)
		out_list.append(temp_eval.copy())
		del temp_eval
	
	# sumup results
	temp_df = pd.DataFrame(out_list)
	temp_dict = dict()
	for item in temp_df.keys():
		temp_array = np.array(list(temp_df[item]))
		out_mean = list(np.average(temp_array, axis=0))
		temp_dict[item] = [round(x, 2) for x in out_mean]
		del out_mean
	
	# add labels
	temp_dict['kmer'].insert(0, "k=" + str(k))
	temp_dict['spaced_sparse'].insert(0, 'sparse')
	temp_dict['spaced_dense'].insert(0, 'dense')
	for item in ['minstrobe', 'randstrobe', 'hybridstrobe']:
		temp_dict[item].insert(0, raw_dict['strobemer_label'])
	
	out_df = pd.DataFrame(temp_dict).transpose()
	out_df.columns = ['label', 'm', 'sc', 'mc', 'E']
	
	return out_df


### test code
def local_test_variable():
	"""
	setup var for local test
	"""
	print("Loading local variables")
	genome_file = '/Users/shaopeng/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/PSU_academic/Koslicki_lab/research/strobemer_truncation/input_genomes/GCA_000014345.1_ASM1434v1_genomic.fna'
	sample_size = 3
	k_max = 48
	k_min = 24
	p_values = [0.01, 0.05, 0.1]
	seq_length = 10000
	wmin=25
	wmax=50
	k_list=[48, 42, 36, 30, 24]


### run analysis
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Strobemer extension analysis, note this is NOT a regular strobemer.",
	                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-g', '--genome', help="Path to input genomes")
	parser.add_argument('-k', '--maxk', type=int, help="Max k value to start", default="48")
	parser.add_argument('-m', '--mink', type=int, help="Minimul k value to start", default="24")
	parser.add_argument('-n', '--sample_size', type=int, help="Number of technical replicates", default=10)
	parser.add_argument('-p', '--p_range', type=str, help="mutation rate range",
	                    default="0.01,0.05,0.1")
	parser.add_argument('-w', '--window_min', type=int, help="Min window size", default=25)
	parser.add_argument('-e', '--window_max', type=int, help="Max window size", default=50)
	parser.add_argument('-l', '--seq_length', type=int, help="Seq length for the test", default=10000)
	
	########################################## read parameters
	args = parser.parse_args()
	
	print("Note!\nThis strobemer code contains genomic location for analysis purpose and only take single string as input.")
	print("This will be different from the formal version strobemer as this contains extra information.")
	
	if not os.path.exists(args.genome):
		raise Exception("Input file %s does not exist." % args.genome)
	genome_file = os.path.abspath(args.genome)
	print("Input genome file is " + genome_file)
	
	sample_size = args.sample_size
	print("Permutation sample size is " + str(sample_size))
	
	k_max = args.maxk
	print("The max k value is ")
	print(k_max)
	k_min = args.mink
	print("The minimul k value is ")
	print(k_min)
	k_list = list(range(k_max, k_min-1, -6))
	
	p_range = args.p_range
	p_values = [float(x) for x in p_range.split(",")]
	print("The range of mutation rate p values are: ")
	print(p_values)
	
	seq_length = args.seq_length
	print("The sequence length to keep is %s" % seq_length)
	wmin = args.window_min
	print("Window size min is %s" % wmin)
	wmax = args.window_max
	print("Window size max is %s" % wmax)

	########################################## input string
	out_seq = prepare_input_string(genome_file=genome_file, seq_length=seq_length)
	
	# unit test
	unit_test_2_ways_strobemer(out_seq[:4000])
	
	########################################## Step2: loop through p and many k values
	for mutation_rate_p in p_values:
		print("Simuation with mutation rate %s" % mutation_rate_p)
		# get a list of sequences
		mutated_seq_list = []
		for i in range(sample_size):
			mutated_seq_list.append(generate_mutated_strings(out_seq, p=mutation_rate_p))
		
		# truncate dict_raw
		for strobe_num in [2, 3]:
			# max kmer dict
			strobe_size = int(k_list[0] / strobe_num)
			temp_raw_dict = s1_generate_all_types_of_kmers(input_string=out_seq, k_value=k_list[0],
			                                               strobe_size=strobe_size, wmin=wmin, wmax=wmax)
			for k in k_list[1:]:
				truncated_raw_dict = truncate_kmer_dict(temp_raw_dict, kmer_truncate=k_list[0] - k,
				                                        single_strobe_truncate=int((k_list[0] - k) / strobe_num))
				# run analysis
				temp_df = s3_workflow_with_npk_from_input_dict_and_mutate_str(raw_dict=truncated_raw_dict,
				                                                              seq_list=mutated_seq_list,
				                                                              strobe_num=strobe_num,
				                                                              p=mutation_rate_p,
				                                                              k=k)
				# save results
				filename = "_".join(
					["Truncated_matrix", "mutation-" + str(mutation_rate_p), "strobe_num-" + str(strobe_num),
					 "k-" + str(k), "rep-" + str(sample_size)]) + ".csv"
				temp_df.to_csv(filename, header=True)
		
		# make regular run
		for k in k_list:
			for strobe_num in [2, 3]:
				strobe_size = int(k / strobe_num)
				temp_raw_dict = s1_generate_all_types_of_kmers(input_string=out_seq, k_value=k, strobe_size=strobe_size,
				                                               wmin=wmin, wmax=wmax)
				# run analysis
				temp_df = s3_workflow_with_npk_from_input_dict_and_mutate_str(raw_dict=temp_raw_dict,
				                                                              seq_list=mutated_seq_list,
				                                                              strobe_num=strobe_num, p=mutation_rate_p,
				                                                              k=k)
				# save results
				filename = "_".join(
					["Regular_matrix", "mutation-" + str(mutation_rate_p), "strobe_num-" + str(strobe_num),
					 "k-" + str(k), "rep-" + str(sample_size)]) + ".csv"
				temp_df.to_csv(filename, header=True)
