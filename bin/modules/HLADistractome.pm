#!/usr/bin/env perl
package HLADistractome;
use strict;
use warnings;
use Getopt::Long;

(my $SCRIPT_NAME = $0) =~ s/.*\///;
my $version = "1.0";
my $creation_date = "1 May 2016";
my $last_updated = "14 Sep 2017";

my $usage="\n$SCRIPT_NAME v$version\n".
	  "\nDESCRIPTIONs:\n" .
	  "This modules creates a distractome by removing all transcripts in the HLA reference from a gene transcript reference fastq using a bed file specifying the genomic regions to be excluded from the distractome.This module can be called from within another perl script or run from the commandline\n\n". 
	  "\nUSAGE:\n" .
	  "perl $SCRIPT_NAME <options>\n" .
	  "\nRequired options:\n" .
	  "-transcript_fa|t\tFasta file containing the nucleotide sequence of the reference transcripts (required)\n" .
	  "-transcript_gtf|g\tGTF file corresponding to the transcript_fa (required)\n" .
	  "-exclude_bed|e\tBed file with the genomics regions to exclude (required)\n" .
	  "\nGeneral options:\n" .
	  "-output_fa|o\tName of distractome output file (default:distractome.fa)\n" .
	  "-log|l\t\tLog file(default:STDERR)\n" .
          "-help|h\t\tDisplays this message\n" .
          "\nAUTHORS:\n" .
	  "Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	  "Chad Brown\n" .
	  "\nCREATED:\n$creation_date\n" .
	  "\nLAST UPDATED:\n$last_updated\n" .
	  "\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
	  "\n";

my %opts=(output_fa=>"distractome.fa", transcript_fa=>"", transscript_gtf=>"");

my %exclude = ();

sub runCommandline{
	GetOptions(\%opts, qw(exclude_bed|e=s transcript_fa|t=s transcript_gtf|g=s output_fa|o=s log|l=s help|h));
	if($opts{help}){
		print STDERR "$usage";
		exit;
	}
	findExcludeGenes($opts{transcript_gtf},$opts{exclude_bed},$opts{log});
	createDistractome($opts{transcript_fa},$opts{output_fa},$opts{log});
}

sub findExcludeGenes{
	my $transcript_gtf = shift;
	$transcript_gtf = shift if ($transcript_gtf eq "HLADistractome");
	my $exclude_bed = shift;	
	my $logfile = shift;
	my %hla;
	my %nonhla;

        my $log;
        if($logfile){
                open($log, ">>", "$logfile") || die "Cannot open log for writing. '$logfile'\n";
        }else{
                $log = *STDERR;
        }
	
	my %genes;
	
	my $exclude_fh;
	if($exclude_bed ne "" && -e $exclude_bed){
		if($exclude_bed=~m/\.gz$/){
			open($exclude_fh, "gunzip -c $exclude_bed |") || (print $log "Cannot unzip and open exclude_bed. '$exclude_bed'\n" && exit 1);
		}else{
			open($exclude_fh, "$exclude_bed");	
		}	
	}else{
		print $log "Gene exclusion bed $exclude_bed does not exist. Exiting\n";
		exit 1;
	}
	
	##Read exclude bed
	my %exclude_coord;
	while(<$exclude_fh>){
		chomp;
		my @cols = split /\t/, $_;
		$exclude_coord{$cols[0]}{$cols[1]}=$cols[2];
	}
	close ($exclude_fh);

	##Read transcript gtf	
	my $transcript_fh;
	if($transcript_gtf ne "" && -e $transcript_gtf){
		if($transcript_gtf=~m/\.gz$/){
			open($transcript_fh, "gunzip -c $transcript_gtf |") || (print $log "Cannot unzip and open transcript gtf. '$transcript_gtf'\n" && exit 1);
		}else{
			open($transcript_fh, "$transcript_gtf");	
		}	
	}else{
		print $log "Transcripts GTF $transcript_gtf does not exist. Exiting\n";
		exit 1;
	}
	
	while(<$transcript_fh>){
		chomp;
		if(substr($_, 0 ,1) ne "#"){
			my @cols = split /\t/, $_;
			if(defined $exclude_coord{$cols[0]}){
					SEARCH:foreach my $start (sort keys %{$exclude_coord{$cols[0]}}){
						#The transcript overlaps or encompasses the exclude region
						if(($cols[3] >= $start && $cols[3] <= $exclude_coord{$cols[0]}{$start}) || ($cols[4] >= $start && $cols[4] <= $exclude_coord{$cols[0]}{$start}) || ($cols[3] <= $start && $cols[4] >= $exclude_coord{$cols[0]}{$start})){
							my @annot = split /;/, $cols[8];
							my %annot_hash = ();
							foreach my $field (@annot){
								$field =~ s/^ //;
								$field =~ s/"//g;
								my @parts = split / /, $field;
								$annot_hash{$parts[0]}=$parts[1];
							}
							my $name = $annot_hash{"gene_name"};
							if($name ne ""){
								$exclude{$name}=1;
								last SEARCH;	
							}
						}
					}
			}
		}
	}
	close($transcript_fh);
	return(\%exclude);
}

sub createDistractome{
	my $transcript_fa = shift;
	$transcript_fa = shift if ($transcript_fa eq "HLADistractome");
	my $output_fa = shift;
	my $logfile = shift;

        my $log;
	my $log_is_file = 0;
        if($logfile){
                open($log, ">>$logfile");
		$log_is_file = 1;
        }else{
                $log = *STDERR;
        }
	
	my $transcript_fh;
	if($transcript_fa ne "" && -e $transcript_fa){
		if($transcript_fa=~m/\.gz$/){
			open($transcript_fh, "gunzip -c $transcript_fa |") || (print $log "Cannot unzip and open transcript file. '$transcript_fa'\n" && exit 1);
		}else{
			open($transcript_fh, "$transcript_fa");	
		}	
	}else{
		print $log "Transcripts fa $transcript_fa does not exist. Exiting\n";
		exit 1;
	}
	my $keep = 0;
	my $sequence = "";
	open(OUT, ">$output_fa") || (print $log "Unable to open $output_fa for writing.Exiting.\n" and exit 1);
	while(<$transcript_fh>){
		chomp;
		if(substr($_,0,1) eq ">"){
			print OUT "$sequence\n" if ($keep == 1);
			$keep=1;
			$sequence="";
			my @parts = split /\|/, $_;
			foreach my $id (@parts){
				if(defined($exclude{$id})){
					$keep = 0;
					last;
				}
			}
			if($keep==1){
				$_=~s/^>/>kraken:taxid\|2\|/;
				print OUT "$_\n";
			}
		}else{
			$sequence .= $_;
		}
	}
	print OUT "$sequence\n" if ($keep == 1);
	close($transcript_fh);
	close(OUT);
	close ($log) if ($log_is_file == 1);
}
__PACKAGE__->runCommandline() unless caller;

