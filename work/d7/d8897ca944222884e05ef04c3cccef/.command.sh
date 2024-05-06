#!/bin/bash -ue
mkdir fastqc_SRR28420797_logs
fastqc -o fastqc_SRR28420797_logs -f fastq -q SRR28420797_1.fastq SRR28420797_2.fastq
