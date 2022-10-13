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

- [Algorithm overview: k-min extension for strobemer](#alg)

  

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

