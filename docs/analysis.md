
# AA

Result files are parsed and transformed to `tidy` format:
- each row presents a junction and pair-wise comparison

at a first moment, we don't filter any results based on 
the probability of obtaining the observed results of a test,
giving a null hypothesis is true (p-value) or effect size change (deltaPSI).
This is done in a second step with files with 'annotated'


Difference in interpretation among the tools:

Analysis design

For Junctionseq model include all conditions at once
Majiq pairwise, Leafcutter currently pairwise, althought it also supports more complex design 

As result, SJ from JunctionSeq have single p-value, while multiple   

- 



Arbitrary cutoffs:
- Majiq: xx == 1 & deltapsi 
- Leafcutter: padj < 0.05
- JunctionSeq:


