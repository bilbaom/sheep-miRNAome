# sheep-miRNAome
This repository contains the scripts used in the manuscript: "The sheep miRNAome: Characterization and distribution of miRNAs in 21 tissues" https://doi.org/10.1016/j.gene.2022.146998

## Scripts
* [mirdeep2-core-command.sh](/mirdeep2-core-command.sh): Code used to run the preprocessing of the samples, genome mapping and the core miRDeep2 algorithm.

* [mirdeep2-quantifier.sh](/mirdeep2-quantifier.sh): Code used to run the miRDeep2 quantifier algorithm.

* [Statistic_analysis.R](/Statistic_analysis.R): Code used for the analysis of miRNA expression and tissue specificity.

* [novel_mirnas.R](/novel_mirnas.R): Code used to give correct "3p" and "5p" names to the filtered miRNAs and to prepare the mature and pre-miRNA fasta files for quantification.  

* [mirna_blast.py](/mirna_blast.py): Code used for sequence conservation analysis of novel miRNAs, selection of unique miRNA sequences for quantification and search of clusters in genome.

* [expression_plots.py](/expression_plots.py): Code used to plot miRNA expression by conservation status, tissue and specificity.
