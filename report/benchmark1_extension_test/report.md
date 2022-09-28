## K-min extension is a valid method for strobemer truncation

### Here are 2 EXCEL files:

1. "**old_result_standard_strobemer.xlsx**": strobe-level truncation results for standard strobemer. It turned out that we CAN NOT directly truncate strobemers because of the **re-shuffling of shorter strobes**. More details can be found [here](https://github.com/ShaopengLiu1/Koslicki_lab_metagenomic_analysis/blob/main/20220905_strobemer_truncation_simulation/report_strobemer_truncation.md).
2. "**kmin_extension_strobemers.xlsx**": results for k-min extension-based strobemers. By fixing a k-min to select strobes, the sampling pattern remains the same (though still being random) for truncated data. And I **found this to be feasible**!

</br>

### Results

Conclusion:

1. Comparing the "Regular p=0.01" matrices (EXCEL 1) with "Extend regular p=0.01" matrices (EXCEL 2) for all k values, we find that $k_{min}$ extension way performs similarly as standard strobmer. **It's a valid method to use**. 
2. However, due to ["tail disturbance"](#tail), the extension way has poor coverage in the tail. This could be a reason why it has smaller match ratio than the standard strobemer. 
3. In EXCEL 2, comparing "Extend regular" VS "Extend trunc", they are almost non-distinguishable, indicating that **the truncation idea can be used**. 
4. In EXCEL 2, we can also find that strobe number of 3 outperforms 2 for when mutation rate p=0.05. n=3 might be used for future analysis.

</br>

Parameters:

| Name            | Value                |
| --------------- | -------------------- |
| k range         | [24, 30, 36, 42, 48] |
| kmin            | 24                   |
| mutation rate p | [0.01, 0.05, 0.1]    |

Method annotation:

| Name           | Method                                                       |
| -------------- | ------------------------------------------------------------ |
| Regular        | Standard strobemer implemented in Py                         |
| Trunc          | Strobe-wise truncation based on standard strobemer           |
| Extend regular | Strobemers built from the k-min extension way                |
| Extend trunc   | Strobe-wise truncation based on k-min extension way of strobmer |











</br>



### Code detail:

Please [check here](https://github.com/KoslickiLab/Strobemer_extension/blob/main/src/benchmark_1_strobemer_extension_analysis.py).

Key points:

1. use $k_{min}-mer$ to determine strobe start points and extend to the real length k. [Check here](https://github.com/KoslickiLab/Strobemer_extension/blob/98bd4398b44a9261bed0962f24795c63f9e9ad66/src/benchmark_1_strobemer_extension_analysis.py#L281) for standard strobemer VS [here](https://github.com/KoslickiLab/Strobemer_extension/blob/98bd4398b44a9261bed0962f24795c63f9e9ad66/src/benchmark_1_strobemer_extension_analysis.py#L453) for this $k_{min}$ extension.
2. this is a validation step for the $k_{min}$ extension method above. In this test, I validated 2 things:
   1. when $k_{min} = k$ (i.e. no extension), the new method will generate exactly same strobemers as a standard strobmer method does
   2. when $k_{min} < k$, the extended strobemers do share expected prefixs
3. "**Tail disturbance**": when the k-mer starting point approaches the tail of the input string, we don't have enough space to deploy a strobemer. The solution ([check here](https://github.com/KoslickiLab/Strobemer_extension/blob/98bd4398b44a9261bed0962f24795c63f9e9ad66/src/benchmark_1_strobemer_extension_analysis.py#L353)) from the paper is to shrink the window size until you get a k-mer (i.e. window width =1). However, the shrink pattern depends on k, which is different in extended strobemers. Therefore, extended strobemers could be different due to window shift. This would **affect at most M k-mers in the range of the last full-size strobemer** (so I call it "tail disturbance"). In reality, I observed 1~2 dozens of mismatches, this could explain the slight difference between 2nd results and 1st results.  <a name="tail"></a>
