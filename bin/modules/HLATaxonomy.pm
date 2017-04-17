#!/usr/bin/env perl
package HLATaxonomy;
use strict;
use warnings;
use Getopt::Long;

(my $SCRIPT_NAME = $0) =~ s/.*\///;
my $version = "1.0";
my $creation_date = "1 May 2016";
my $last_updated = "11 Jan 2017";

my $usage="\n$SCRIPT_NAME v$version\n".
	  "\nDESCRIPTIONs:\n" .
	  "This module creates an HLA taxonomy based on allele names in the supplied IMGT reference. This module can be run from commandline or directly through function calls from another perl script.\n\n" . 
	  "\nUSAGE:\n" .
	  "perl $SCRIPT_NAME <options>\n" .
	  "\nRequired options:\n" .
	  "-reference_fa|r\t\tFasta file containing the nucleotide sequences of alleles in the hla reference (required)\n" .
	  "\nGeneral options:\n" .
	  "-output_dir|o\tName of taxonomy output directory (default:'.')\n" .
	  "-log|l\t\tLog file(default:STDERR)\n" .
          "-help|h\t\tDisplays this message\n" .
          "\nAUTHORS:\n" .
	  "Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	  "Chad Brown:chad.brown\@q2labsolutions.com\n" .
	  "\nCREATED:\n$creation_date\n" .
	  "\nLAST UPDATED:\n$last_updated\n" .
	  "\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
	  "\n";

sub runCommandline{
	my %opts = (fasta=>"", output_dir=>".", log=>"");
	GetOptions(\%opts, qw(help|h reference_fasta|r=s output_dir|o=s log|l=s));
	if($opts{help}){
		print $usage;
		exit;
	}elsif($opts{reference_fasta} eq ""){
		print STDERR "Please include a reference fasta. $usage";
		exit 1
	}
	createHLATaxonomy($opts{reference_fasta}, $opts{output_dir}, $opts{log});
}

sub createHLATaxonomy{
	my $fa = shift;
	$fa = shift if ($fa eq "HLATaxonomy");
	my $output_dir = shift || ".";
	my $logfile = shift;

        my $log;
	my $log_is_file = 0;
        if($logfile && $logfile ne ""){
                open($log, ">>$logfile");
		$log_is_file = 1;
        }else{
                $log = *STDERR;
        }
		
	if($fa eq "" || !-e $fa || -z $fa ){
		print $log "Please enter a fastq file that exists and is not empty. Exiting\n";
		exit 1;
	}
	
	if($output_dir eq "."){
		print $log "WARNING: No output directory specified. Using . instead.\n"; 
	}
	
	if(-e $output_dir){
		if(-d $output_dir){
			if(-e "$output_dir/taxonomy"){
				if(-d "$output_dir/taxonomy"){
				}else{
					print $log "\n$usage\n$output_dir/taxonomy exists but is not a directory. Exiting.\n";
					exit 1;
				}
			}else{
				mkdir "$output_dir/taxonomy";
			}	
			if(-e "$output_dir/data"){
				if(-d "$output_dir/data"){
					if(-e "$output_dir/data/reference"){
						if(-d "$output_dir/data/reference"){	
						}else{
							print $log "\n$usage\n$output_dir/data/reference exists but is not a directory. Exiting.\n";
							exit 1;
						}	
					}else{
						mkdir "$output_dir/data/reference";
					}
				}else{
					print $log "\n$usage\n$output_dir/data exists but is not a directory. Exiting\n";
					exit 1;
				}
			}else{
				mkdir "$output_dir/data";
				mkdir "$output_dir/data/reference";
			}	
			
		}else{
			print $log "\n$usage\n$output_dir exists but is not a directory. Exiting.\n";
			
		}
	}else{
		mkdir "$output_dir";	
		mkdir "$output_dir/taxonomy";
		mkdir "$output_dir/data";
		mkdir "$output_dir/data/reference";
	}
	##This assigns known HLA genes to a class
	my %proteins = (A=>"ClassI", B=>"ClassI", C=>"ClassI", DMA=>"ClassII", DMB=>"ClassII", DOA=>"ClassII", DOB=>"ClassII", DPA1=>"ClassII",DPA2=>"ClassII",DPA2=>"ClassII", DPB1=>"ClassII", DPB2=>"ClassII", DQA1=>"ClassII", DQB1=>"ClassII", DQB2=>"ClassII", DQB3=>"ClassII", DRA=>"ClassII", DRB1=>"ClassII", DRB2=>"ClassII", DRB3=>"ClassII", DRB4=>"ClassII", DRB5=>"ClassII", DRB6=>"ClassII", DRB7=>"ClassII", DRB8=>"ClassII", DRB9=>"ClassII", E=>"ClassI", F=>"ClassI", G=>"ClassI", H=>"ClassI", HFE=>"Non-HLA", J=>"ClassI", K=>"ClassI", L=>"ClassI", MICA=>"Non-HLA", MICB=>"Non-HLA", TAP1=>"Non-HLA", TAP2=>"Non-HLA", V=>"ClassI", N=>"ClassI", P=>"ClassI",S=>"ClassI",T=>"ClassI",U=>"ClassI",W=>"ClassI",X=>"ClassI",Y=>"ClassI",Z=>"ClassI", PSMB9=>"Non-HLA", PSMB8=>"Non-HLA",MICC=>"Non-HLA",MICD=>"Non-HLA",MICE=>"Non-HLA");
	
	##NCBI taxonomy file formats (from NCBI website)
	#---------names.dmp------------------------------------------------------------------------------
	#Taxonomy names file has these fields:
	#tax_id					-- the id of node associated with this name
	#name_txt				-- name itself
	#unique name				-- the unique variant of this name if name not unique
	#name class				-- (synonym, common name, ...)
	
	#--------nodes.dmp--------------------------------------------------------------------------------
	#This file represents taxonomy nodes. The description for each node includes the following fields:
	#tax_id					-- node id in GenBank taxonomy database
 	#parent tax_id				-- parent node id in GenBank taxonomy database
 	#rank					-- rank of this node (superkingdom, kingdom, ...) 
 	#embl code				-- locus-name prefix; not unique
 	#division id				-- see division.dmp file
 	#inherited div flag  (1 or 0)		-- 1 if node inherits division from parent
 	#genetic code id				-- see gencode.dmp file
 	#inherited GC  flag  (1 or 0)		-- 1 if node inherits genetic code from parent
 	#mitochondrial genetic code id		-- see gencode.dmp file
 	#inherited MGC flag  (1 or 0)		-- 1 if node inherits mitochondrial gencode from parent
 	#GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
 	#hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
 	#comments				-- free-text comments and citations
	#-------------------------------------------------------------------------------------------------


	##Establish the root of the tree
	open(NAMES, ">$output_dir/taxonomy/names.dmp");
	open(NODES, ">$output_dir/taxonomy/nodes.dmp");
	my %nodes;
	my $taxonomy_start=1;
	print NAMES  "$taxonomy_start\t|\troot\t|\tHLA:root\t|\tscientific name\n";
	print NODES  "$taxonomy_start\t|\t0\t|\troot\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n";
	$taxonomy_start++;
	print NAMES  "$taxonomy_start\t|\tDistractome\t|\tHLA:distractome\t|\tscientific name\n";
	print NODES  "$taxonomy_start\t|\t1\t|\tsuperkingdom\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n";
	$taxonomy_start++;
	my $hla_node = $taxonomy_start;
	print NAMES  "$taxonomy_start\t|\tHLA alleles\t|\tHLA:alleles\t|\tscientific name\n";
	print NODES  "$taxonomy_start\t|\t1\t|\thlaroot\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n";
	$taxonomy_start++;
	print NAMES  "$taxonomy_start\t|\tClassI\t|\tHLA:ClassI\t|\tscientific name\n";
	print NODES  "$taxonomy_start\t|\t$hla_node\t|\tsuperkingdom\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n";
	$nodes{"ClassI"}=$taxonomy_start;
	$taxonomy_start++;
	print NAMES  "$taxonomy_start\t|\tClassII\t|\tHLA:ClassII\t|\tscientific name\n";
	print NODES  "$taxonomy_start\t|\t$hla_node\t|\tsuperkingdom\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n";
	$nodes{"ClassII"}=$taxonomy_start;
	$taxonomy_start++;
	print NAMES  "$taxonomy_start\t|\tNon-HLA\t|\tHLA:Non-HLA\t|\tscientific name\n";
	print NODES  "$taxonomy_start\t|\t$hla_node\t|\tsuperkingdom\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n";
	$nodes{"Non-HLA"}=$taxonomy_start;
	$taxonomy_start++;
	
	
	my %tree;
	my %leaves;
	open(FA, ">$output_dir/data/reference/hla.ref.forKraken.fa");
	open(IN, $fa);
	while(<IN>){
		chomp;
		if(substr($_,0,1) eq ">"){
			my $name = substr($_,5); #Removes >HLA: or >IPD: from the fasta header

			#Adding kraken:taxid is necessary for kraken to work
			print FA ">kraken:taxid|$taxonomy_start|HLA:$name\n";
			print NAMES  "$taxonomy_start\t|\t$name\t|\tHLA:$name\t|\tscientific name\n";
			my @parts = split / /, $name;
			my @aparts = split /\*/,$parts[1];
			my @pparts = split /:/, $aparts[1];
			my $class = $proteins{$aparts[0]} || "Non-HLA";
			
			#$parts[0] is the accession
			#$parts[1] is the allele name
			#$aparts[1] is the gene name
			#@pparts contains the digits of precision 

			if($#pparts <0){ #Gene level allele
				$tree{$class}{$aparts[0]}=$parts[1];
			}elsif($#pparts ==0){ #2-digit allele
				$tree{$class}{$aparts[0]}{$pparts[0]}=$parts[1];
			}elsif($#pparts ==1){ #4-digit allele
				$tree{$class}{$aparts[0]}{$pparts[0]}{$pparts[1]}=$parts[1];
			}elsif($#pparts ==2){ #6-digit allele
				$tree{$class}{$aparts[0]}{$pparts[0]}{$pparts[1]}{$pparts[2]}=$parts[1];
			}elsif($#pparts ==3){ #8-digit allele
				$tree{$class}{$aparts[0]}{$pparts[0]}{$pparts[1]}{$pparts[2]}{$pparts[3]}=$parts[1];
			}			
			$leaves{$parts[1]}=$taxonomy_start;
			$taxonomy_start++;
		}else{
			print FA "$_\n";
		}
	}
	close(IN);
	close(FA);

	#This iterates through the tree to create the Taxonomy
	#$protein = Gene: kingdom level
	#$d1 = Allele grop: phylum level  2-digit
	#$d2 = HLA proteind: family level 4-digit
	#$d3 = synonomous substitutions: genus level 6-digit
	#$d4 = non-coding changes: species level 8-digit
	
	foreach my $class (keys %tree){
		my %phash = %{$tree{$class}};
		foreach my $protein (sort keys %phash){
			if(ref($phash{$protein}) eq "HASH"){ 
				print NAMES  "$taxonomy_start\t|\t$protein\t|\tHLA:Protein $protein\t|\tscientific name\n";
				print NODES  "$taxonomy_start\t|\t$nodes{$class}\t|\tkingdom\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n";
				$nodes{$protein}=$taxonomy_start;
				my $ptax = $taxonomy_start;
				$taxonomy_start++;
				
				my %d1hash = %{$phash{$protein}};
				foreach my $d1 (sort keys %d1hash){
					if(ref($d1hash{$d1}) eq "HASH"){
						print NAMES  "$taxonomy_start\t|\t$protein:$d1\t|\tHLA:Protein $protein:$d1\t|\tscientific name\n";
						print NODES  "$taxonomy_start\t|\t$ptax\t|\tphylum\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n";
						$nodes{$d1}=$taxonomy_start;
						my $d1tax = $taxonomy_start;
						$taxonomy_start++;
						my %d2hash = %{$d1hash{$d1}};
						foreach my $d2 (sort keys %d2hash){
							if(ref($d2hash{$d2}) eq "HASH"){
								print NAMES  "$taxonomy_start\t|\t$protein:$d1:$d2\t|\tHLA:Protein $protein:$d1:$d2\t|\tscientific name\n";
								print NODES  "$taxonomy_start\t|\t$d1tax\t|\tfamily\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n";
								$nodes{$d2}=$taxonomy_start;
								my $d2tax = $taxonomy_start;
								$taxonomy_start++;
								
								my %d3hash = %{$d2hash{$d2}};
								foreach my $d3 (sort keys %d3hash){
									if(ref($d3hash{$d3}) eq "HASH"){
										print NAMES  "$taxonomy_start\t|\t$protein:$d1:$d2:$d3\t|\tHLA:Protein $protein:$d1:$d2:$d3\t|\tscientific name\n";
										print NODES  "$taxonomy_start\t|\t$d2tax\t|\tgenus\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n";
										$nodes{$d3}=$taxonomy_start;
										my $d3tax = $taxonomy_start;
										$taxonomy_start++;
										
										my %d4hash = %{$d3hash{$d3}};
										foreach my $d4 (sort keys %d4hash){
											if(ref($d4hash{$d4}) eq "HASH"){
												print $log "\nMore than 8 digit precision detected $protein:$d1:$d2:$d3:$d4\n";
											}else{
												my $name = $d4hash{$d4};
												my $taxid = $leaves{$name};
												print NODES  "$taxid\t|\t$d3tax\t|\tspecies\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n";
											}	
										}
									}else{
										my $name = $d3hash{$d3};
										my $taxid = $leaves{$name};
										print NODES  "$taxid\t|\t$d2tax\t|\tgenus\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n";
									}	
								}
							}else{
								my $name = $d2hash{$d2};
								my $taxid = $leaves{$name};
								print NODES  "$taxid\t|\t$d1tax\t|\tfamily\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n";
							}
						}
					}else{
						my $name = $d1hash{$d1};
						my $taxid = $leaves{$name};
						print NODES  "$taxid\t|\t$ptax\t|\tphylum\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n";
					}
				}
			}else{
				my $name = $phash{$protein};
				my $taxid = $leaves{$name};
				print NODES  "$taxid\t|\t$nodes{$class}\t|\tkingdom\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n";
			}
		} 
	}
	close(NAMES);
	close(NODES);
	close ($log) if ($log_is_file);
}

__PACKAGE__->runCommandline() unless caller;
