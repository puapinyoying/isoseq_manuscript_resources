#!/bin/bash

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.primary_assembly.annotation.gtf.gz
gunzip gencode.vM10.primary_assembly.annotation.gtf.gz
ln -s gencode.vM10.primary_assembly.annotation.gtf mm10.gtf