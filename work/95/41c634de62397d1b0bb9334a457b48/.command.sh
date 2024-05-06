#!/bin/bash -ue
mkdir fastqc_SRR28420795_logs
fastqc -o fastqc_SRR28420795_logs -f fastq -q SRR28420795_1.fastq SRR28420795_2.fastq
