# ADV-PCC-Final-Project
Final project for AS.410.712.82.FA23

* ABOUT *

This web interface tool is designed to analyze the ratio of nonsynonymous (dN) to synonymous (dS) mutations in protein-coding genes for pairwise analysis. The user will upload their own pairwise sequence alignment file of their genomes of choice.

Source code can be obtained here:
http://bfx3.aap.jhu.edu/mnguye87/final_dnds_tool.tar.gz

Demo output with pre-uploaded files can be seen in the mySQL MariaDB [mnguye87_chado]. The table final_project5, show the results from the test of PWA files from /test_pwa_files/

Utilizes tools published at:

https://www.ncbi.nlm.nih.gov/gene/

	* FASTA files for pairwise sequence alignment can be found here

https://www.ebi.ac.uk/Tools/psa/

	* These pairwise sequence alignment tools are recommended to be used
	* The resulting aligned files should be downloaded in Pearson/FASTA format

* DETAILED USAGE *

(1) Choose a gene
(2) Download DNA/nucleotide sequence for two different species from NCBI GeneBank in FASTA format
(3) Perform global/local pairwise sequence alignment (PSA)
(4) Upload your aligned file in Pearson/FASTA format
