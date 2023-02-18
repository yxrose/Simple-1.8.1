#!/bin/bash
#pipeline for mapping EMS mutants
#variables file
trait=HD
#define the path to Java java-1.8.0 version
java='/public/home/zmhan/anaconda3/bin/java'

#input files
mut_files=fastq/${trait}*mut*
wt_files=fastq/${trait}*wt*

#output names
mutation=recessive #change to dominant if the mutation is dominant

line=${trait}  ##if you prefer, change EMS to the name of your line.  Letters and underscores only.
mut=${line}_mut
wt=${line}_wt



my_species=Oryza_sativa_Japonica #paste your species name here to replace Arabidopsis_thaliana

SNPEFF="~/software/SNPEFF/snpEff/snpEff.jar"
annName="Japonica7.0"



