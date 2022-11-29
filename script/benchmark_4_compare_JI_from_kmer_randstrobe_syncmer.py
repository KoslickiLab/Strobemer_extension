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
import matplotlib.pyplot as plt
import seaborn as sb

from script.benchmark_1_pilot_ext_randstrobe_vs_regular_randstrobe import generate_string_list_from_genome, generate_mutated_strings

def get_kmer_set_from_str(seq, ksize):
	"""
	Generate a kmer set from input string
	"""
	temp = set()
	for i in range(len(seq) - ksize + 1):
		temp.add(seq[i:i + ksize])
	return temp

def get_randstrobe_set_from_str(seq, n, l, wmin, wmax, prime):
	"""
	Get a randstrobe set from input string
	"""
	temp = set()
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
			temp.add(fix_seq)
			continue

		# ranstrobe: 1st strobe
		rand_seq = seq[i:i + l]
		last_strobe = seq[i:i + l]
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
			# add this strobe to randstrobe
			rand_seq += seq[current_pos:current_pos + l]
			
		# End of current cycle, add the current randstrobe to the set
		temp.add(rand_seq)
	
	return temp
	
def is_kmer_an_open_syncmer_both_strand_for_submer(input_string: str, submer_length: int):
	"""
	Determine if a given kmer (string) is an OPEN syncmer at given submer_length via lexicographic order.
	Consider submer for both directions (same as Strobe-align)
	:param input_string: a kmer
	:param submer_length: the length substring in this kmer for syncmer selection
	"""
	temp_start = input_string[:submer_length]
	temp_start = min(temp_start, khmer.reverse_complement(temp_start))
	for i in range(len(input_string) - submer_length + 1):
		another_submer = input_string[i:i+submer_length]
		another_submer = min(another_submer, khmer.reverse_complement(another_submer))
		if temp_start > another_submer:
			return False
	# o.w., the start is minimum, so this is an open syncmer
	return True

def get_syncmer_set_from_kmer_set(kmer_set, submer_length: int):
	"""
	Get a syncmer set from input string
	"""
	temp = set()
	for kmer in kmer_set:
		if is_kmer_an_open_syncmer_both_strand_for_submer(kmer, submer_length):
			temp.add(kmer)
	
	return temp
	


### main object to store 3 types of k-mers and calculate JI
class simple_kmer_list_for_seq:
	def __init__(self, seq: str, k: int=20, submer_length:int=15, n: int=2, l: int=10, wmin: int=25, wmax: int=50, prime: int=6229):
		"""
		Create a simple kmer object to store 3 types of k-mers and calculate JI
		"""
		self.seq = seq # input sequence (a single string)
		self.k = k  # kmer length
		self.n = n  # strobe number
		self.l = l  # strobe length
		self.wmin = wmin # strobe window min
		self.wmax = wmax # strobe window max
		self.prime = prime # prime for randstrobe
		self.submer_length = submer_length # submer for syncmer
		### kmer lists:
		self.kmer_set = get_kmer_set_from_str(seq, ksize=k)
		self.randstrobe_set = get_randstrobe_set_from_str(seq, n=n, l=l, wmin=wmin, wmax=wmax, prime=prime)
		self.syncmer_set = get_syncmer_set_from_kmer_set(self.kmer_set, submer_length=submer_length)
		
	def calculate_ji_with_another_obj(self, other):
		kmer_ji = 1.0 * len(self.kmer_set.intersection(other.kmer_set)) / len(self.kmer_set.union(other.kmer_set))
		randstrobe_ji = 1.0 * len(self.randstrobe_set.intersection(other.randstrobe_set)) / len(self.randstrobe_set.union(other.randstrobe_set))
		syncmer_ji = 1.0 * len(self.syncmer_set.intersection(other.syncmer_set)) / len(self.syncmer_set.union(other.syncmer_set))
		
		results = [round(kmer_ji, 3), round(randstrobe_ji,3), round(syncmer_ji, 3)]
		return results
	

### local variables
def local_test_variable():
	"""
	setup var for local test
	"""
	print("Loading local variables")
	genome_file = '/Users/shaopeng/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/PSU_academic/Koslicki_lab/research/strobemer_truncation/input_genomes/GCA_000014345.1_ASM1434v1_genomic.fna'
	sample_size = 10
	k_value = 20
	p_values = [0.01, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2, 0.3]
	seq_length = 30000
	wmin=25
	wmax=50


### run analysis
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Strobemer simulation",
	                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-g', '--genome', help="Path to input genomes")
	parser.add_argument('-k', '--kvalue', type=int, help="K value to use", default="20")
	parser.add_argument('-n', '--sample_size', type=int, help="Number of technical replicates", default=10)
	parser.add_argument('-p', '--p_range', type=str, help="mutation rate range",
	                    default="0.01,0.03,0.05,0.08,0.1,0.15,0.2,0.3")
	parser.add_argument('-w', '--window_min', type=int, help="Min window size", default=25)
	parser.add_argument('-e', '--window_max', type=int, help="Max window size", default=50)
	parser.add_argument('-l', '--seq_length', type=int, help="Seq length for the test", default=30000)
	
	########################################## read parameters
	args = parser.parse_args()
	
	if not os.path.exists(args.genome):
		raise Exception("Input file %s does not exist." % args.genome)
	genome_file = os.path.abspath(args.genome)
	print("Input genome file is " + genome_file)
	
	sample_size = args.sample_size
	print("Permutation sample size is " + str(sample_size))
	
	k_value = args.kvalue
	print("The k value is ")
	print(k_value)
	
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
	
	
	########################################## Step2: loop through mutation rate p and rep 10 times for JI calculation
	original_obj = simple_kmer_list_for_seq(seq=out_seq)
	out_dict = {}
	
	for mutation_rate_p in p_values:
		out_ji = [0, 0, 0] # init
		for i in range(sample_size):
			mutated_seq = generate_mutated_strings(input_string=out_seq, p=mutation_rate_p)
			mutated_obj = simple_kmer_list_for_seq(seq=mutated_seq)
			ji_out = original_obj.calculate_ji_with_another_obj(mutated_obj)
			out_ji = [sum(x) for x in zip(out_ji, ji_out)]
		# get mean for sample_size times
		out_dict[str(mutation_rate_p)] = [round(x/sample_size,3) for x in out_ji]
		
	final_table = pd.DataFrame(out_dict)
	final_table.index = ['kmer_JI', 'randstrobe_JI', 'syncmer_JI']
	final_table.to_csv("Compare_JI_of_3_types_of_kmers.csv", header=True)
	
	# make plot
	out_df = final_table.transpose()
	out_df['mut_p'] = out_df.index
	plot_df = out_df.melt('mut_p', var_name="cols", value_name="JI values")
	
	fig, axs = plt.subplots(1, 1)
	sb.pointplot(x="mut_p", y="JI values", hue='cols', data=plot_df, ax=axs)
	axs.set_title("Comparison of JIs for 3 types of kmers")
	fig.savefig("_".join(["Compare_JI_3_types_of_kmers.png"]), dpi=300)
	plt.close(fig)
	
	
	
