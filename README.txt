Always use this GitHb README.txt version (and not the Plant Physiology one).

# SIMPLE — A SIMPLE pipeline for mapping point mutations

#####   SYSTEM REQUIREMENTS  #####

# Runs on mac (10.11.6) with X11 and Linux (tested for centOS 6.7)

# Requires Java JKD 1.8 (8u131; current url: https://www.oracle.com/java/technologies/javase/javase8-archive-downloads.html)
# Requires Curl
# Requires bwa
# Requires samtools
# Requires snpEff

# Requires R and the following packages: ggplot2, reshape2 and ggrepel
# requires Command Line Tools in Mac OS X (for instructions: https://apple.stackexchange.com/questions/254380/why-am-i-getting-an-invalid-active-developer-path-when-attempting-to-use-git-a)

# internet connection

##### RUNNING SIMPLE #####

1. Download and unpack the Simple package from the following link: https://github.com/yxrose/Simple-1.8.1.git (A fork and revised version) by pressing the green link: "Clone or download"

2. Place the Simple folder in your home directory.

3. Rename your fastq files as follow. For the mutant and WT bulks, the names should have the following formats, trait_X.mut.RY.fastq and trait_X.wt.RY.fastq, respectively (note the dots); trait represent any trait like FT (flowering time), X is the line name (but you can also omit the "X." part of the name) and "Y" is either 1 or 2. E.g., if you sent a pair-end mutant bulk DNA and single-end WT bulk DNA, you could name the three files as follow: FT_lineX.mut.R1.fastq, FT_lineX.mut.R2.fastq and FT_lineX.wt.R1.fastq (the wild_cards for recognize files is : FT*mut*, and FT*wt*). Use only letters and underscores.

4. Place the renamed fastq files in the fastq folder located in the Simple folder.

5. Open the folder "scripts" inside Simple; open the data_base.txt file. Locate your species in the first column and copy it.
Open the file simple_variables.sh inside the folder "scripts" with a text editor and paste the species name you've just copied to replace "Arabidopsis_thaliana" as the species name (e.g., this line should look like: my_species=Arabidopsis_thaliana or my_species=Oryza_sativa_Japonica). Save the file.

6. In simple_variables.sh, assign a trait name that correspond to your renamed fastq file. For this instance, it's "FT". trait=FT.

7. In the same file, simple_variables.sh choose the path to your Java1.8 executable in line 6. Note:This should be a absolute pathway start from root.

8. If the mutation you are mapping is dominant, change the mutation from recessive to dominant in the simple_variables.sh file.

9. Download your reference genome (both fasta and gff file ) to refs/ directory, and rename your genome file to be {species name}.fa, for example: Oryza_sativa_Japonica.fa.

10. Use snpEff software to make a local annotation database based on your {species name}.fa file. Look the manual of snpEff for detail.Then open simple_variables.sh, assign you snpEff pathway to variable SNPEFF and assign you annotated genome database name to variable annDB.

11. Open the Terminal application.

12. Type: cd ~/Simple. Press return.

13. Type: chmod +x ./scripts/simple.sh. Press return.

14. Type: ./scripts/simple.sh. Press return.

15. The last command will execute the program.

16. The script will run for a few hours up to a couple of days, depending on the size of your fastq files and the size of the genome you are working with. You will know it finished once the prompt shows the following colorful text: “Simple is done”.

17. Output files (EMS might be the name of the prefix you chose above for your output file names. For example, FT.).
	a. EMS.allSNPs.txt: All the SNPs that Simple found compared to the background genome, usually caused by the mutagenesis process.
	b. EMS.candidates.txt: The genes that are likely to be the causal mutation; notice that there can be more than one gene due to linkage. Another case is that there are no genes. In this scenario, browse the EMS.allSNPs.txt file for genes with high score in the ratio column, that affect the protein (e.g., missense, stop_gained) and that have reasonable coverage (columns mut.ref, mut.alt, wt.ref and wt.alt).
	c. EMS.Rplot_allele.pdf: Representation of the allelic distribution in the WT and mutant bulks. This file is useless and you can ignore it.
	d. EMS.Rplot.loess.1.pdf and EMS.Rplot.loess.3.pdf: Graphical representation of the results with allelic ratio column above 0.1 and 0.3, respectively. Your gene is likely to show in one of the highest peaks with a y-axis value of ~0.66 and if it appears in the EMS.candidates.txt file, it will also be marked in these plots.

#######################Example#########################

Let's say that you want to perform a forward genetic screen and find  mutations that promote wings development in dragons. You collect 5000 dragon eggs, mutagenize them and identify an M2 mutant line with some individual dragons that have no wings. Next, You extract two DNA bulks, one from the wingless dragons and one from their flying siblings and send each DNA prep for PE (paired ends) WGS (you can also have a SE but make sure you have approximately 40x coverage). Once you receive your results, unzip the files, rename them wing_have.wt.R1.fq, wing_have.wt.R2.fq, wing_no.mut.R1.fq, wing_no.mut.R2.fq (fq and fastq suffix are both OK) and place these files in the fastq folder. Open the file data_base.txt in the folder scripts and copy the species name, i.e., Dragon_fly to replace Arabidopsis_thaliana in the simple_variables.sh file (located in the same folder). Type the path to the Java1.8 executable. Since you noticed that in this M2 mutant line approximately 0.75 of the dragons have no wings, you assume the mutation is dominant and change the string "recessive" to dominant in the simple_variables.sh file (line 13). Follow the instructions 9-14 above to complete the process. It is that simple.



















