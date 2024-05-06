#!/bin/bash -ue
salmon quant --threads 4 --libType=u -i salmon_index -1 SRR28420796_1.fastq -2 SRR28420796_2.fastq -o SRR28420796
