#!/usr/local/bin/perl
use strict;
use warnings ;
#
# Institution: The University of Memphis
# Project: Mouse Cerebellum
# Written By: Evan Savage
# Edit Date: 05/20/14
#
#
# This program takes the TF motifs/matrix found in proliferation query that matched with SOM 
# and runs them against our known TFs to produce matches for possible genes.
#
# The purpose of this program is to allow for a more robust GCAT/GeneIndexer analysis of genes
# associated with proliferation during development of the mouse cerebellum.
#
#
# Initializing some variables
#
my @temp = () ;
my @temp2 = () ;
my @sig_motifs = ([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[])  ;
my $sig_motif = () ;
my $match_motif = () ;
my $match_genes = () ;
my $i = () ;
#
# Creating files for editing/pulling info and matching. EDIT TF_lIST for correct geneindexer search.
#
my $infile = "C:/Users/Evan/Desktop/Bioinf Master's Program/Graduate Assistantship Materials and Reports/Cerebellum Study/tf_sig_match_literal.csv" ;
my $tf_list = "C:/Users/Evan/Desktop/Bioinf Master's Program/Graduate Assistantship Materials and Reports/Cerebellum Study/motif2TFmapping.csv" ;
my $output =  "C:/Users/Evan/Desktop/Bioinf Master's Program/Graduate Assistantship Materials and Reports/Cerebellum Study/tf_match_output_literal.txt" ;
#
# Opening the files
#
open (SIG_MOTIF, $infile) or die "Unable to open $infile: $!" ;
#
#!!!!!!!!!!!!!!!!!!!!!!!Pushing the motifs into an array, WE NEED TO ALTER THE NUMBER OF AoA for each search (LITERAL VS CONCEPTUAL !!!!!!!!!!
#
while ($sig_motif = <SIG_MOTIF> ) {
	chomp $sig_motif ;
	@temp = split(",", $sig_motif) ;
	push (@{$sig_motifs[0]}, $temp[0]) ;
	push (@{$sig_motifs[1]}, $temp[1]) ;
	push (@{$sig_motifs[2]}, $temp[7]) ;
	push (@{$sig_motifs[3]}, $temp[8]) ;
	push (@{$sig_motifs[4]}, $temp[9]) ;
	push (@{$sig_motifs[5]}, $temp[10]) ;
	push (@{$sig_motifs[6]}, $temp[11]) ;
	push (@{$sig_motifs[7]}, $temp[12]) ;
	push (@{$sig_motifs[8]}, $temp[13]) ;
	push (@{$sig_motifs[9]}, $temp[14]) ;
	push (@{$sig_motifs[10]}, $temp[15]) ;
	push (@{$sig_motifs[11]}, $temp[16]) ;
	push (@{$sig_motifs[12]}, $temp[17]) ;
	push (@{$sig_motifs[13]}, $temp[18]) ;
	push (@{$sig_motifs[14]}, $temp[19]) ;
	push (@{$sig_motifs[15]}, $temp[20]) ;
	push (@{$sig_motifs[16]}, $temp[21]) ;
	push (@{$sig_motifs[17]}, $temp[22]) ;
	push (@{$sig_motifs[18]}, $temp[23]) ;
	push (@{$sig_motifs[19]}, $temp[24]) ;
	push (@{$sig_motifs[20]}, $temp[25]) ;
	push (@{$sig_motifs[21]}, $temp[26]) ;
	push (@{$sig_motifs[22]}, $temp[27]) ;
	push (@{$sig_motifs[23]}, $temp[28]) ;
	push (@{$sig_motifs[24]}, $temp[29]) ;
	push (@{$sig_motifs[25]}, $temp[30]) ;
	push (@{$sig_motifs[26]}, $temp[31]) ;
	push (@{$sig_motifs[27]}, $temp[32]) ;
	push (@{$sig_motifs[28]}, $temp[33]) ;
	push (@{$sig_motifs[29]}, $temp[34]) ;
	push (@{$sig_motifs[30]}, $temp[35]) ;
	push (@{$sig_motifs[31]}, $temp[36]) ;
	push (@{$sig_motifs[32]}, $temp[37]) ;
	push (@{$sig_motifs[33]}, $temp[38]) ;
	push (@{$sig_motifs[34]}, $temp[39]) ;
	push (@{$sig_motifs[35]}, $temp[40]) ;
	push (@{$sig_motifs[36]}, $temp[41]) ;
	push (@{$sig_motifs[37]}, $temp[42]) ;
	push (@{$sig_motifs[38]}, $temp[43]) ;
	

}
close SIG_MOTIF ;
#
# Opening output file
#
open (OUTPUT, ">", $output) or die "Unable to open > $output: $!" ;
print OUTPUT "Motif,TF,P-value,hits in ENSG00000163739 CXCL1,hits in ENSG00000118513 MYB,hits in ENSG00000137720 C11orf1,hits in ENSG00000101057 MYBL2,hits in ENSG00000123473 STIL,hits in ENSG00000173207 CKS1B,hits in ENSG00000166508 MCM7,hits in ENSG00000166803 KIAA0101,hits in ENSG00000163808 KIF15,hits in ENSG00000110484 SCGB2A2,hits in ENSG00000213347 MXD3,hits in ENSG00000165304 MELK,hits in ENSG00000091128 LAMB4,hits in ENSG00000115738 ID2,hits in ENSG00000105810 CDK6,hits in ENSG00000144354 CDCA7,hits in ENSG00000112238 PRDM13,hits in ENSG00000125319 C17orf53,hits in ENSG00000177426 TGIF1,hits in ENSG00000057657 PRDM1,hits in ENSG00000137804 NUSAP1,hits in ENSG00000042493 CAPG,hits in ENSG00000111145 ELK3,hits in ENSG00000111206 FOXM1,hits in ENSG00000107105 ELAVL2,hits in ENSG00000019549 SNAI2,hits in ENSG00000106511 MEOX2,hits in ENSG00000103257 SLC7A5,hits in ENSG00000185811 IKZF1,hits in ENSG00000134825 C11orf10,hits in ENSG00000165891 E2F7,hits in ENSG00000164611 PTTG2,hits in ENSG00000178252 WDR6,hits in ENSG00000179348 GATA2,hits in ENSG00000139263 LRIG3,hits in ENSG00000121957 GPSM2,hits in ENSG00000112312 GMNN\n" ;
#
# Opening the TF_list to iterate through for our TFs
#
open (TF_LIST, $tf_list) or die "Unable to open $tf_list: $!:" ;
#
# Iterating through each line to find out motif and corresponding TFs
#
while (my $line = <TF_LIST>) {
	chomp $line ;
	@temp2 = split(",",$line) ;
	for $i (0..$#{$sig_motifs[0]}) {
		if ($sig_motifs[0][$i] =~ /$temp2[0]/i) {
			print "Found match at TF: $temp2[0]\n" ;
			print OUTPUT "$temp2[0],$temp2[1],$sig_motifs[1][$i],$sig_motifs[2][$i],$sig_motifs[3][$i],$sig_motifs[4][$i],$sig_motifs[5][$i],$sig_motifs[6][$i],$sig_motifs[7][$i],$sig_motifs[8][$i],$sig_motifs[9][$i],$sig_motifs[10][$i],$sig_motifs[11][$i],$sig_motifs[12][$i],$sig_motifs[13][$i],$sig_motifs[14][$i],$sig_motifs[15][$i],$sig_motifs[16][$i],$sig_motifs[17][$i],$sig_motifs[18][$i],$sig_motifs[19][$i],$sig_motifs[20][$i],$sig_motifs[21][$i],$sig_motifs[22][$i],$sig_motifs[23][$i],$sig_motifs[24][$i],$sig_motifs[25][$i],$sig_motifs[26][$i],$sig_motifs[27][$i],$sig_motifs[28][$i],$sig_motifs[29][$i],$sig_motifs[30][$i],$sig_motifs[31][$i],$sig_motifs[32][$i],$sig_motifs[33][$i],$sig_motifs[34][$i],$sig_motifs[35][$i],$sig_motifs[36][$i],$sig_motifs[37][$i],$sig_motifs[38][$i]\n" ;
		}
	}
}
close TF_LIST ;
close OUTPUT ;
	