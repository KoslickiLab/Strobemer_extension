from src.extension_strobemer_obj import *
import pytest

print("Pipe start")

def test_object_performance():
	input_file = "test_data/toy_data.fna"
	raw_seq_dict = generate_string_dict_from_genome(input_file)
	temp_obj = ext_strobemer_obj(k=30, n=2, l=15, smin=10, label="test", filename="none")
	temp_obj.load_seq_dict(raw_seq_dict)
	temp_obj.build_kmer_dict_with_start_pos()
	temp_obj.build_both_randstrobe_dict()
	temp_obj.build_open_syncmer()
	temp_obj.print_info()
	#temp_obj.export_to_pkl()
	
	# check 2 version of randstrobe
	test_ext = ext_strobemer_obj(k=30, n=2, l=15, smin=15, label="test")
	test_ext.load_seq_dict(raw_seq_dict)
	test_ext.build_kmer_dict_with_start_pos()
	test_ext.build_both_randstrobe_dict()
	# smin = l means a regular strobemer
	assert test_ext.ext_randstrobe_dict == temp_obj.regu_randstrobe_dict
	assert test_ext.ext_randstrobe_dict == test_ext.regu_randstrobe_dict
	del test_ext
	
	# trie practice
	some_manual_prefix = 'CTAATGGAGAAACTCAT'
	assert len(temp_obj.kmer_tst.keys(some_manual_prefix)) > 0
	# find kmer in trie that is a prefix of a given seq
	some_unique_kmer = 'GTCTTGCAATAATGGCAAAACTAAATGTAC'
	assert temp_obj.kmer_tst.prefixes(some_unique_kmer + 'aaaaa') == [some_unique_kmer]
	
	# kmer and prefix lookup
	temp_obj.kmer_query('GTCTTGCAATAATGGCAAAACTAAATGTAC', kmer_type='kmer')
	temp_obj.kmer_query('GTCTTGCAATAATGGCAAAACTAAATGTAC', kmer_type='rs')
	temp_obj.kmer_query('GTCTTGCAATAATGGCAAAACTAAATGTAC', kmer_type='ers')
	temp_obj.prefix_query('GTCTTA', kmer_type='kmer')
	temp_obj.prefix_query('GTCTTA', kmer_type='ers')
	temp_obj.prefix_query('GTCTTGCAATAATGGCAAAACTAAATGTAC', kmer_type='kmer')