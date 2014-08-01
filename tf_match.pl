#!/usr/local/bin/perl
use strict;
use warnings ;
use Data::Dumper;

#
# Institution: The University of Memphis
# Project: Mouse Cerebellum
# Written By: Ethan Willis
# Edit Date: 7/31/2014 
#
#
# This program takes the TF motifs/matrix found in proliferation query that matched with SOM 
# and runs them against our known TFs to produce matches for possible genes.
#
# The purpose of this program is to allow for a more robust GCAT/GeneIndexer analysis of genes
# associated with proliferation during development of the mouse cerebellum.
#

# Creating files for editing/pulling info and matching. EDIT TF_lIST for correct geneindexer search.
my $infile_motif = "data/input/SOM_1_1_TF_motifs_sig.csv";
my $infile_TFMap = "data/input/motif2TFmapping.csv";
my $output = "data/output/tfmatch.markdown";


# Read our tf-motif mapping and motif signature files into matrices and associative arrays.
my %TF_map = &readTFMap($infile_TFMap);
my %motif_sigs = &readMotifSigFile($infile_motif);

# Find the transcription factors that our motifs map to..
my %tf_motif_matches = &findTFMotifMatches(\%TF_map, {%motif_sigs});

# Save results
&outputTFMotifMatches({%tf_motif_matches}, $output);


# Give a tf-motif mapping and a list of motif sigs, this subroutine finds all tf matches for each motif.
sub findTFMotifMatches {
	# Define our local variables and map appropriate ones to passed parameters.
	my(%matches);
	my %TF_map = %{$_[0]};
	my %motif_sigs = %{$_[1]};
	# Lookup each motif in our tf-map and store tuples of motifs and tf strings.
	foreach my $motif_name (keys %motif_sigs) {
		# if we have a tf mapped to this motif..
		if( exists($TF_map{$motif_name}) && $TF_map{$motif_name} ne '-') {
			# Store this mapping in our matches.
			$matches{$motif_name} = $TF_map{$motif_name};
		}
	}
	return %matches;
}

# Outputs a map of motif-tf matches.
sub outputTFMotifMatches {
	our %matches; 
	local *matches = %{$_[0]};
 	my $output;
	$output = $_[1];	
	open(MATCHES, ">", $output) 
		or die "Unable to open output file: $output";
	foreach my $motif (keys %tf_motif_matches) {
		print MATCHES $motif . "," .  $tf_motif_matches{$motif} . "\n";
	}
}

# TODO-{NEEDS WORK}  Subroutine to print out some debugging text
sub printDebug {
	print Dumper(\%TF_map);
	foreach my $key (keys %motif_sigs) {
		print $key . "=>";
		print Dumper( $motif_sigs{$key});
	}
}
# Given a list of motif sigs and a tf-motif mapping this subroutine finds all tf matches for each motif.

# Reads in the Motif signature file. 
sub readMotifSigFile {
	# Read in parameters passed to this function
	my($motif_sig_file, %motif_sigs, @cur_sig, $motif_name);
	($motif_sig_file) = @_;

	open(MOTIF_SIGS, "<", $motif_sig_file)
		or die "Unable to open Motif Signature file: $motif_sig_file";

	# Read and process our motif sig file line by line
	my $dropHeader = <MOTIF_SIGS>;	
	while(<MOTIF_SIGS>) {
		# Get our line
		chomp;

		# Build an array out of the current line split around the delimeter ','
		@cur_sig = split(',', $_); 
		 
		# Splice our current signature into two arrays. One array containing the motif name, the other 
		# containing all the data for the motif.
		$motif_name = splice(@cur_sig, 0, 1);
		# Store copy of array in hash 
		$motif_sigs{$motif_name} = [@cur_sig];	
	}	
	close(MOTIF_SIGS);

	return %motif_sigs;	
}

# Reads in a TF Map and returns the associative array.
# parameters for this function 
# TF_file - String - The location of the TF mapping file.
sub readTFMap {
	# Read in parameters passed to this function.
	my($TF_file, %TF_map, $motif_name, $tf_name);
	($TF_file) = @_;
	
	# Initialize variables
	%TF_map = ();

	open(TF_MOTIF_MAP, "<", $TF_file)
		or die "Unable to open TF-Motif Mapping file: $TF_file";

	# Read and process our tf-motif map file line by line
	my $dropHeader = <TF_MOTIF_MAP>;
	while(<TF_MOTIF_MAP>) {
		chomp;
		# split our current line, $_, by ',' delimiters.
		# store each element into the respective variables.
		($motif_name, $tf_name) = split(',', $_);
		$TF_map{$motif_name} = $tf_name;
	}	 

	close(TF_MOTIF_MAP);

	return %TF_map;
}


	
