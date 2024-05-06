#!/bin/bash -ue
mkdir fastqc_SRR28420796_logs
fastqc -o fastqc_SRR28420796_logs -f fastq -q SRR28420796_1.fastq SRR28420796_2.fastq
