#!/usr/bin/env perl 
package SequenceFunctions;
use strict;
use warnings;
use Getopt::Long;

(my $SCRIPT_NAME = $0) =~ s/.*\///;
my $version = "1.0";
my $creation_date = "1 Nov 2016";
my $last_updated = "11 Jan 2017";

my $usage = "\n$SCRIPT_NAME v$version\n" .
	    "\nDESCRIPTION:\n" .
	    "This module is designed to perform operations on nucleotide sequences but some operations can also work on other strings such as quality score strings. This module can be accessed through function calls or run from the commandline as well.\n" .
	    "\nUSAGE:\n" .
	    "perl $SCRIPT_NAME <options>\n" .
	    "\nRequired options:\n" .
	    "-revcomp|rc\tReverse complements the sequence. Only supported for nucleotide sequences. If the nucleotide sequence contains U's these will be complemented to A but A's will always be complemented to T not U.\n" .
	    "-rev|rv\t\tReverse complements the sequence or string.\n" .
	    "-sequence|s\tNucleotide sequence (required)\n" .
	    "\nGeneral options:\n" .
	    "-help|h\t\t\tDisplays this message\n" .
            "\nAUTHORS:\n" .
	    "Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	    "Chad Brown\n" .
	    "\nCREATED:\n$creation_date\n" .
	    "\nLAST UPDATED:\n$last_updated\n" .
	    "\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
	    "\n";

sub runCommandline{
	my %opts = (sequence=>"");
	GetOptions(\%opts, qw(help|h revcomp|rc rev|rv sequence|s=s));
	if($opts{help}){
		print $usage;
		exit;
	}elsif($opts{sequence} eq ""){
		print STDERR "Please enter a sequence\n$usage";
		exit 1;
	}else{
		my $operation = 0;
		if ($opts{revcomp}){
			$operation = 1;
			print revcomp($opts{sequence}) . "\n";
		}
		if ($opts{rev}){
			$operation = 1;
			print rev($opts{sequence}) , "\n";
		}
		if($operation == 0){
			print STDERR "Please select one or more operations to perform.\n$usage";
			exit 1;
		}
	}
}

sub revcomp{
	my $read = shift;
	$read = shift if ($read eq "SequenceFunctions");
	my $revcomp ="";
	foreach my $base (split //, $read){
		$base=~tr/ATUCGatucg/TAAGCtaagc/;
		$revcomp = $base . $revcomp;
	}
	return $revcomp;
}

sub rev{
	my $read = shift;
	$read = shift if ($read eq "SequenceFunctions");
	my $rev ="";
	foreach my $base (split //, $read){
		$rev = $base . $rev;
	}
	return $rev;
}

__PACKAGE__->runCommandline unless caller;
