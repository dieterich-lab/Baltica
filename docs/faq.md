### What you mean with `baltica provides a interface between users and DJU methods`:
There are many specificities to the DJU methods, and while running one method is not too complicated, figuring out how to run multiple methods is time-demanding. Baltica aims to facilitate this task so that methods results can be produced and compared. 

### Snakemake `Error: Directory cannot be locked.`:
This error happens when there is an error or failure during the workflow execution, and Snakemake's process does not have the opportunity to unlock the directory.
Use `baltica <workflow> <config> --unlock` to resolve it. See more [here](https://snakemake.readthedocs.io/en/stable/project_info/faq.html#how-does-snakemake-lock-the-working-directory)

### ERROR lines in majiq_gtf_to_gff:
Most of the errors in this rules are not fatal. See [here](https://manpages.debian.org/unstable/gbrowse/gtf2gff3.1p.en.html) for the diagnostics. 


### rmats empty or mostly empty outputs 
This error can be either issue with:
- [The read length parameter](https://github.com/Xinglab/rmats-turbo/issues/95). To resolve it, increase the read length parameter or use `--variable-read-length` in Baltica configuration.
    - Or an error with [the stack size limit](https://github.com/Xinglab/rmats-turbo/issues/91). To resolve it increase the stack size in bash with `ulimit -c unlimited`.

### Junctionseq bpapply error:

This error occurs in the `junctionseq_analysis` rule,  which uses multiple threads with BiocParallel. One can overcome this issue by setting threads to 1 on the said rule.

### Junctionseq analysis thread error

```Exception in thread "main" java.lang.ArrayIndexOutOfBoundsException: Index xxx out of bounds for length xxx
```
It is complaining the maximum read length is longer than the read length input. First, check the maximum read length in the quality control report and then increase the `read_len` parameter on Baltica config.

<!-- ### How is the integration, annotation and AS assignment process occur?

Integrations proceeds as the following:
1. Find the overlapping genomic ranges allowing for a maximum difference of 2 nucleotides on the start and stop to allow differences in the genomic coordinate system
1. Detect groups of introns that fit the criteria in 1
1. 

Integration, then  annotation.

Uses introns in the *de novo* annotation to annotate the introns coming from multiple sources. Only introns matching this novel annotation are annotated.
 -->
### How does Baltica compute the score for each DJU method? 

The different DJU methods are pretty different in many aspects, including how they compute the final test statistic, and we use the following rule to compute the score for the Baltica table (higher is better):  
- **majiq** score = 1 - non-changing-threshold (probability of |&Delta;&Psi;| > 0.2, by default)  
- **leafcutter** = 1 - p.adjust  
- **junctionseq** = 1 - padjust  
- **rmats** = 1 - FDR  


### I see a message: `/bin/bash: /root/.bashrc: Permission denied`. What is wrong?
This error message is benign and, in our experience, does not affect workflow execution. 