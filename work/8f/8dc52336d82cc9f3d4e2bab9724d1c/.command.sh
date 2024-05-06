#!/bin/bash -ue
salmon quant --threads 4 --libType=u -i salmon_index -1 SRR28420797_1.fastq -2 SRR28420797_2.fastq -o SRR28420797
