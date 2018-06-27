#!/usr/bin/env perl
package TaxonomyDivisions;
use strict;
use warnings;
use Getopt::Long;

(my $SCRIPT_NAME = $0) =~ s/.*\///;
my $version = "1.0";
my $creation_date = "1 May 2016";
my $last_updated = "18 Jan 2018";

my $usage="\n$SCRIPT_NAME v$version\n".
	  "\nDESCRIPTIONs:\n" .
	  "This module creates an HLA taxonomy division file used as input to the EA-modified Kraken and specifies how reads classified to each taxonomy is to be distributed into filtered fastq files. This module can be run from commandline or directly through function calls from another perl script.\n\n" . 
	  "\nUSAGE:\n" .
	  "perl $SCRIPT_NAME <options>\n" .
	  "\nRequired options:\n" .
	  "-names|n\tTaxonomy names file (names.dmp) (required)\n" .
	  "-nodes|t\tTaxonomy nodes file (nodes.dmp) (required)\n" .
	  "\nGeneral options:\n" .
	  "-output|o\tName of the division output file (default:'divisions.txt')\n" .
	  "-log|l\t\tLog file(default:STDERR)\n" .
          "-help|h\t\tDisplays this message\n" .
          "\nAUTHORS:\n" .
	  "Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	  "Chad Brown\n" .
	  "\nCREATED:\n$creation_date\n" .
	  "\nLAST UPDATED:\n$last_updated\n" .
	  "\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
	  "\n";

sub runCommandline{
	my %opts = (names=>"", nodes=>"", output=>"divisions.txt");
	GetOptions(\%opts, qw(names|n=s nodes|t=s output|o=s log|l=s help|h));
	if($opts{help}){
		print $usage;
		exit;
	}elsif($opts{names} eq ""){
		print STDERR "Please specify a taxonomy names file. $usage";
		exit 1;
	}elsif($opts{nodes} eq ""){
		print STDERR "Please specify a taxonomy nodes file. $usage";
		exit 1;
	}else{
		createTaxonomyDivisions($opts{names}, $opts{nodes}, $opts{output}, $opts{logfile});	
	}			
}

sub createTaxonomyDivisions{
	my $names_file = shift;
	$names_file = shift if ($names_file eq "TaxonomyDivisions");
	my $nodes_file = shift;
	my $output_file = shift || "divisions.txt";
	my $logfile = shift;

        my $log;
        if($logfile){
                open($log, ">$logfile");
        }else{
                $log = *STDERR;
        }

	if($names_file eq ""){
		print $log "Please specify a names.dmp file. Exiting.\n";
		exit 1;
	}elsif(! -e $names_file && ! -s $names_file){
		print $log "Names file $names_file does not exist or is empty. Please specify an existing and non-empty names.dmp file. Exiting.\n";
		exit 1;
	}

	if($nodes_file eq ""){
		print $log "Please specify a nodes.dmp file. Exiting.\n";
		exit 1;
	}elsif(! -e $nodes_file && ! -s $nodes_file){
		print $log "Nodes file $nodes_file does not exist or is empty. Please specify an existing and non-empty nodes.dmp file. Exiting.\n";
		exit 1;
	}

	my %names=();
	$names{0}="root";
	open(OUT, ">$output_file") || (print $log "Cannot open file $output_file for writing. Exiting.\n" && exit 1);

	open(IN, $names_file) || (print $log "Cannot open file $names_file for reading. Exiting.\n" && exit 1);
	while(<IN>){
		chomp;
		my @cols = split /\t/, $_;
		$names{$cols[0]}=$cols[2];
	}
	close(IN);


	my %classes = (root=>-3,hlaroot=>-2,superkingdom=>-1,kingdom=>0,phylum=>1,family=>2,genus=>3,species=>4,unknown=>5);
	my %designation = ();
	my %parents = ();
	open(IN, $nodes_file) || (print $log "Cannot open file $nodes_file for reading. Exiting.\n" && exit 1);
	while(<IN>){
		chomp;
		my @cols = split /\t/, $_;
		my $class= $cols[4] || "unknown";
		##Stores the designation level
		$designation{$cols[0]}=$classes{$class};
		##Stores the parent
		$parents{$cols[0]}=$cols[2];
	}
	close (IN);
	$designation{0}=$classes{root};

	my %divisions;
	foreach my $node(sort keys %parents){
		my $parent = $parents{$node};
		if($parent ne ""){
			##Travel up the tree until the ancestor at the kingdom level is found
			while($designation{$parent}>0){
				$parent = $parents{$parent};
			}
			##If below parent level then place into kingdom division
			if($designation{$node}>0){
				$divisions{$names{$parent}}{$node}=1 if($designation{$parent} ne "");
			}elsif($designation{$node} == 0){ ##The node is at the kingdom
				$divisions{$names{$node}}{$node}=1;
			}
		}else{ 
			##While this should be rare. Place any orphaned kingdom nodes in their own devision
			$divisions{$names{$node}}{$node}=1 if($designation{$node} == 0);	
		}
	}
	 
	foreach my $division (sort keys %divisions){
		my @nodes = sort {$a<=>$b} keys %{$divisions{$division}};
		print OUT "$division\t$nodes[0]";
		for(my $i=1; $i<=$#nodes; $i++){
			print OUT ",$nodes[$i]";
		}
		print OUT "\n";
	}

	close(OUT);
}
__PACKAGE__->runCommandline() unless caller;
