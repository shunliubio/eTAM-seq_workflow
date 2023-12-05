# Workflow of eTAM-seq data processing

## Dependencies:

- **Cutadapt**, tested with 1.18
- **Samtools**, tested with 1.14
- **HISAT-3N**, tested with 2.2.1-3n
- **pileup2var**, tested with 1.1.0
- **R**, tested with 4.1.0. Required package: eTAMseq

`eTAMseq` can be downloaded from the repository and installed using the following command in R:

```R
install.packages("eTAMseq_x.x.x.tar.gz", repos = NULL, type = "source")
```

## Citation:

Xiao, Y., Liu S., Ge R., et al. Transcriptome-Wide Profiling and Quantification of *N*<sup>6</sup>-methyladenosine by Enzyme-Assisted Adenosine Deamination. Nat Biotechnol. 2023 Jul;41(7):9931003. [PMID: 36593412](https://pubmed.ncbi.nlm.nih.gov/36593412)
