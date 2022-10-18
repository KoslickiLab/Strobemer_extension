# Apply the truncation idea to Strobemer

### Current tasks:

- [ ] Init a genome object to incorporate the following features: k-mer, randstrobe, extension-randstrobe, syncmer, and extension-syncmer
- [ ] Write a document for the object usage



### Brainstorm

1. Conclusion from benchmark tests:
   1. **Strobemer == Extension_strobemer**: evaluation matrices (match, sequence coverage etc.) are similar in strobemer and extension-strobemer
   2. **Strobemer always have better sequence coverage** than k-mer; and its match ratio will outperform k-mer when mutation rate p is high (>0.2)

2. Future directions:
   1. Apply the extension idea to syncmer and integrate it into Strobealign
   2. Apply the extension idea to randstrobe and use the adjustable sequence match to compare sequences (show superority than k-mer)
   3. Relate extension-randstrobe-based CI to k-mer based CI. Similar to the Mash Distance plot, try if we can find high correlation between the CIs. If so, some theory works might be added. 



### Contents

- [Reproducibility](#repro)

  - [Usage (TBD)](#usage)

  - [K-mer object document](#obj)

- [Algorithm overview: k-min extension for strobemer](#alg)

  

</br>



### Reproducibility <a name="repro"></a>

---

#### Create environment

```
conda create -n ext_strobemer -c bioconda -c conda-forge -c conda --file ./src/requirements.txt python=3.7
conda activate ext_strobemer
```

</br>

#### Usage <a name="usage"></a>

Build scripts to perform specific analysis.

</br>

#### K-mer object document <a name="obj"></a>

Basic functions of the k-mer object for extension-based randstrobe analysis. To accommodate the purpose of sequence alignment, the sequence name (e.g. “>chr1”) and k-mer locations are both stored as values in the k-mer dict. 

Parameters (also attributes)

| Keyword  | Default value | Information                                         |
| -------- | ------------- | --------------------------------------------------- |
| k        |               | k-mer length                                        |
| n        |               | number of strobes                                   |
| l        |               | length of a single strobe                           |
| smin     | 5             | length of sub-strobe for extension-based randstrobe |
| wmin     | 25            | min window size for strobemer                       |
| wmax     | 50            | max window size for strobemer                       |
| prime    | 6229          | prime number to mod for randstrobe selection        |
| label    | "no_label"    | keyword for output files                            |
| filename | "no_filename" | store the corresponding genome file if needed       |

</br>

Functions:

| Function name                  | Purpose                                                      |
| ------------------------------ | ------------------------------------------------------------ |
| print_info                     | Print parameters and status of k-mer dict of this object     |
| load_seq_dict                  | Load a dict of [seq_id : sequence] pairs into this object    |
| build_kmer_dict_with_start_pos | Build k-mer dict of [k-mer : seq_id + location] pairs and save k-mers into a ternary search tree |
| build_both_randstrobe_dict     | 1. Build regular randstrobe dict of [k-mer : seq_id + location] pairs <br />2. Build extension-based randstrobe dict of [k-mer : seq_id + location] pairs<br />3. Save k-mers from extension-based randstrobe into a ternary search tree for prefix lookup |
| export_to_pkl                  | Save this object into a pickle object                        |
| kmer_query                     | Check if a k-mer exists in this object                       |
| prefix_query                   | Check if a prefix exists in this object                      |

</br>

[Additional attributes](https://github.com/KoslickiLab/Strobemer_extension/blob/902eab72676788b25ebadf2f387a65d8e9a23b2f/src/extension_strobemer_obj.py#L93)

| Name                 | Contents                                              |
| -------------------- | ----------------------------------------------------- |
| seq_dict             | a dict for [seq_id : sequence] pairs                  |
| kmer_dict            | a dict for [k-mer : seq_id + location] pairs          |
| kmer_tst             | a ternary search tree for all k-mers                  |
| regu_randstrobe_dict | a dict for [randstrobe : seq_id + location] pairs     |
| ext_randstrobe_dict  | a dict for [ext_randstrobe : seq_id + location] pairs |
| ext_randstrobe_tst   | a ternary search tree for all ext_randstrobes         |
| (To be added)        |                                                       |



</br>

Examples

```
from src.extension_strobemer_obj import *

# use the toy data
input_file = "test/test_data/toy_data.fna"

# build a dict of [seq_id : sequence] pair from fasta file. It's necessary for multi-contig genomes.
# the reason why I leave this step outside the object is for comparison of mutated genomes
raw_seq_dict = generate_string_dict_from_genome(input_file)

# build a extension_strobemer object with k=30 and 2 strobes (each of length 15). By default, window size wmin=25, wmax=50
temp_obj = ext_strobemer_obj(k=30, n=2, l=15, smin=10, label="test", filename="none")

# load sequences into the object
temp_obj.load_seq_dict(raw_seq_dict)

# build k-mer dict
temp_obj.build_kmer_dict_with_start_pos()

# build randstrobe and extension-randstrobe dict
temp_obj.build_both_randstrobe_dict()

# print the stats of this object
temp_obj.print_info()


### Now we can perform some simple searching of k-mers or prefixes
# 'ers' for extension-randstrobe
# 'rs' for randstrobe
# 'kmer' for regular k-mer

some_manual_prefix = 'CTAATGGAGAAACTCAT'
# now only support prefix search in [kmer, ers]
temp_obj.prefix_query(some_manual_prefix, kmer_type='ers') #this is empty
temp_obj.prefix_query(some_manual_prefix, kmer_type='kmer')

temp_obj.kmer_query('GTCTTGCAATAATGGCAAAACTAAATGTAC', kmer_type='kmer')
temp_obj.kmer_query('GTCTTGCAATAATGGCAAAACTAAATGTAC', kmer_type='rs')
temp_obj.kmer_query('GTCTTGCAATAATGGCAAAACTAAATGTAC', kmer_type='ers')
```







</br>

## Algorithm overview <a name="alg"></a>
A figure illustration:
<img src="https://github.com/KoslickiLab/Strobemer_extension/blob/main/figure/extension_example.png" alt="darkmatter_plan" style="zoom:33%;" />

1. Select $k_{min} \text{ and } k_{max} \leq W_{min}$, the extension/truncation only happens in this region

2. Build strobemers from reference genomes on $k=k_{max}$:

   a) generate strobemer with parameter $(n, k_{min}, W_{min}, W_{max})$

   b) extend each strobe from $k_{min}$ to $k_{max}$

   c) note that the 2 steps above can be done at the same time as 2a) actually select strobe starting point

3. Build KTST via strobe concatenation

   a) concatenate all $k_{min}-mers$ from each strobe

   b) for all remaining letters from each strobe, we iteratively append 1 letter from each strobe to the sequence in a)

   c) step b) ensures that truncation will get a prefix match for strobemers

   d) then we can truncate this KTST by step of n instead of 1 (n is number of strobes) to get shorter strobemers

4. For arbitrary query data with interested $k-mer$ length $k=k' \in [k_{min}, k_{max}]$, we can build its strobemer in the same way

   a) generate strobemer with parameter $(n, k_{min}, W_{min}, W_{max})$

   b) extend each strobe from $k_{min}$ to $k'$

   c) this ensures that the sampling protocol are the same in query and ref though being random
   

</br>

