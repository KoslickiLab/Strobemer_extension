# Strobemer extension
Apply the idea of CMash to Strobemer and benchmark the performance.
- [Algorithm overview: k-min extension for strobemer](#alg)
- [Py code implementation](#code)
- [Benchmark 1: evaluation matrix for extended strobemers vs regular strobemers](#eval)
- [Benchmark 2: influence of k value on different mutation rate](#optk)

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

## (TBD) Py code details <a name="code"></a>

</br>

## (TBD) Benchmark 1: extension performance  <a name="eval"></a>

</br>

## (TBD) Benchmark 2: influence of k value with different mutation rate <a name="optk"></a>
</br>
