#!/bin/bash -ue
mkdir fastqc_SRR28420798_logs
fastqc -o fastqc_SRR28420798_logs -f fastq -q SRR28420798_1.fastq SRR28420798_2.fastq
