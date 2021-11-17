#!/bin/bash

# The quantification of known and novel miRNAs in all the samples.

quantifier.pl -p /Qallpremirna.fa -m /Qallmirna.fa -r reads/all_samples_reads.fa -y now -P -j

