# sheep-miRNAome
This repository contains the scripts used in the manuscript: [under revision]

## Scripts
* [mirdeep2-core-command.sh](/mirdeep2-core-command.sh): Code used to run the preprocessing of the samples, genome mapping and the core miRDeep2 algorithm.

* [mirdeep2-quantifier.sh](/mirdeep2-quantifier.sh): Code used to run the miRDeep2 quantifier algorithm.

* [Statistic analysis.R](/code.R): Code used for the analysis of miRNA expression and tissue specificity.

* [novel_mirnas.R](/code2.R): Code used to rename correctly the filtered miRNAs and to prepare the mature and pre-miRNA fastas for quantification.  

* [mirna_blast.py](/mirna_blast.py): Code used for sequence conservation analysis of novel miRNAs, selection of unique miRNA sequences for quantification and search of clusters in genome.

* [expression_plots.py](/expression_plots.py): Code used to plot miRNA expression by conservation status, tissue and specificity.
