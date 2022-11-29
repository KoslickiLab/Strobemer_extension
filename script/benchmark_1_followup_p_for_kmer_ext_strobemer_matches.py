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
import sys
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Strobemer/kmer match vs p.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-g', '--genome', help="Path to input genomes")
	parser.add_argument('-k', '--maxk', type=int, help="Max k value to start", default=30)
	parser.add_argument('-m', '--mink', type=int, help="Minimul k value to start", default=30)
	parser.add_argument('-n', '--sample_size', type=int, help="Number of technical replicates", default=10)
	parser.add_argument('-p', '--p_range', type=str, help="mutation rate range",
	                    default="0.01,0.05,0.08,0.11,0.14,0.17,0.2,0.23,0.26,0.30")
	parser.add_argument('-w', '--window_min', type=int, help="Min window size", default=25)
	parser.add_argument('-e', '--window_max', type=int, help="Max window size", default=50)
	parser.add_argument('-l', '--seq_length', type=int, help="Seq length for the test", default=10000)
	
	########################################## read parameters
	args = parser.parse_args()
	# load functions from benchmark_1 file
	src_folder = os.path.realpath(os.path.dirname(__file__))
	sys.path.insert(1, src_folder)
	print("Loading object functions form %s" % src_folder)
	from benchmark_1_strobemer_extension_analysis import *
	
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
	k_list = list(range(k_max, k_min - 1, -6))
	print(k_list)
	
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
	
	
	########################################## Step2: loop through p and many k values
	# make regular run
	strobe_num = 2
	
	for mutation_rate_p in p_values:
		print("Simuation with mutation rate %s" % mutation_rate_p)
		# get a list of sequences
		mutated_seq_list = []
		for i in range(sample_size):
			mutated_seq_list.append(generate_mutated_strings(out_seq, p=mutation_rate_p))
		

		for k in k_list:
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

	### summary
	for k in k_list:
		kmer_match = []
		kmer_mc = []
		strobe_match = []
		strobe_mc = []
		
		for mutation_rate_p in p_values:
			filename = "_".join(
				["Regular_matrix", "mutation-" + str(mutation_rate_p), "strobe_num-" + str(strobe_num),
				 "k-" + str(k), "rep-" + str(sample_size)]) + ".csv"
			temp_df = pd.read_csv(filename, header=0, index_col=0)
			kmer_match.append(temp_df.loc['kmer', 'm'])
			kmer_mc.append(temp_df.loc['kmer','mc'])
			strobe_match.append(temp_df.loc['randstrobe','m'])
			strobe_mc.append(temp_df.loc['randstrobe','mc'])
			
		out_dict = {'p':p_values, 'kmer_match':kmer_match, 'strobe_match':strobe_match, 'kmer_mc':kmer_mc, 'strobe_mc':strobe_mc}
		out_df = pd.DataFrame(out_dict)
		out_df.to_csv("_".join(["Summary_k", str(k), "comparison.csv"]), header=True, index=False)
		# make plot
		plot_df = out_df.melt('p', var_name="cols", value_name="percent")
		
		fig, axs = plt.subplots(1, 1)
		sb.pointplot(x="p", y="percent", hue='cols', data=plot_df, ax=axs)
		axs.set_title("Comparison of kmer/randstrobe")
		fig.savefig("_".join(["Summary_k", str(k), "comparison.png"]), dpi=300)
		plt.close(fig)
		