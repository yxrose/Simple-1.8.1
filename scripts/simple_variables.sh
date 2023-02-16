
#!/bin/bash
#pipeline for mapping EMS mutants
#variables file
trait=FT
#define the path to Java java-1.8.0 version
java='~/anaconda3/bin/java'

#input files
mut_files=fastq/${trait}*mut*
wt_files=fastq/${trait}*wt*

#output names
mutation=recessive #change to dominant if the mutation is dominant

line=${trait}  ##if you prefer, change EMS to the name of your line.  Letters a\
nd underscores only.
mut=${line}_mut
wt=${line}_wt

