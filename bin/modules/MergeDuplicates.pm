#!/usr/bin/env perl
package MergeDuplicates;
use strict;
use warnings;
use Getopt::Long;

(my $SCRIPT_NAME = $0) =~ s/.*\///;
my $version = "1.0";
my $creation_date = "1 May 2016";
my $last_updated = "12 Jan 2017";

my $usage="\n$SCRIPT_NAME v$version\n".
	  "\nDESCRIPTIONs:\n" .
	  "This module creates merges duplicate sequences in a fasta file and creates an allele map file containing information about the merged alleles and whether each allele is common. This module can be run from commandline or directly through function calls from another perl script.\n\n" . 
	  "\nUSAGE:\n" .
	  "perl $SCRIPT_NAME <options>\n" .
	  "\nRequired options:\n" .
	  "-reference_fa|r\tFasta file containing the nucleotide sequences of alleles in the hla reference (required)\n" .
	  "-output_prefix|o\tPrefix of output file (required)\n" .
	  "-cwd_alleles|c\tLocation of a list of Common and Well Documented HLA alleles(required)\n" .
	  "\nGeneral options:\n" .
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
	my %opts = (reference_fasta=>"", output_prefix=>"", cwd_alleles=>"");
	GetOptions(\%opts, qw(help|h reference_fasta|r=s output_prefix|o=s cwd_alleles|c=s log|l=s));
	if($opts{help}){
		print $usage;
		exit;
	}elsif($opts{reference_fasta} eq ""){
		print STDERR "Please specify a reference fasta. $usage";
		exit 1;
	}elsif($opts{output_prefix} eq ""){
		print STDERR "Please specify an output_prefix. $usage";
		exit 1;
	}elsif($opts{cwd_alleles} eq ""){
		print STDERR "Please specify the location of a list of cwd alleles. $usage";
		exit 1;
	}else{
		mergeDuplicates($opts{reference_fasta}, $opts{output_prefix}, $opts{cwd_alleles}, $opts{log});
	}
	

}

sub mergeDuplicates{
	my $reference_fasta = shift || "";
	$reference_fasta = shift if($reference_fasta eq "MergeDuplicates");
	my $output_prefix = shift || "";
	my $cwd_alleles = shift || "";
	my $logfile = shift || "";
	
	my $log;
	my $log_is_file = 0; 
	if($logfile ne ""){
		open($log, ">>$logfile");
		$log_is_file = 1;
	}else{
		$log = *STDERR;
	}
	

	##The CWD allele file will be used to create a column specifying if each allele is common and well documented
	my %cwd = ();
	if(-e $cwd_alleles){
		open(CWD, $cwd_alleles);
		while(<CWD>){
			chomp;
			my @cols = split /\t/, $_;
			$cwd{$cols[0]}=1;
		}	
		close CWD;
	}else{
		print $log "CWD allele file $cwd_alleles does not exist.\n";
		exit 1; 
	}

	if(-e $reference_fasta){
		open(IN, $reference_fasta) || (print $log "Could not open reference fasta: $reference_fasta\n" && exit 1);
	}else{
		print $log "Reference fasta $reference_fasta does not exist.\n";
		exit 1; 
	}

	my $name="";
	my $seq;
	my %sequences;
	my $order=0;

	##Read in fasta file and store in hash based on seq
	while(<IN>){
		chomp;
		if(substr($_,0,1) eq ">"){	
			$_ =~ s/>//;
			if($name ne ""){
				$sequences{$seq}{order}=$order;
				$sequences{$seq}{cnt}++;
				my @cols = split /[\t ]/, $name;
				$cols[0] =~ s/HLA://;
				$cols[0] =~ s/IPD://;
				my $acwd = $cwd{$cols[0]}||0;
				$sequences{$seq}{cwd} = $acwd if(!defined($sequences{$seq}{cwd}) || $sequences{$seq}{cwd} != 1);
				$sequences{$seq}{name}.="$name,";
			}
			$order++;	
			$name=$_;
			$seq = "";
		}else{
			$seq .= $_;
		}	
	}
	$sequences{$seq}{order}=$order;
	$sequences{$seq}{cnt}++;
	$sequences{$seq}{name}.="$name,";
	my @cols = split /[\t ]/, $name;
	$cols[0] =~ s/HLA://;
	$cols[0] =~ s/IPD://;
	my $acwd = $cwd{$cols[0]}||0;
	$sequences{$seq}{cwd} = $acwd if(!defined($sequences{$seq}{cwd}) || $sequences{$seq}{cwd} != 1);
	
	close(IN);
	open(OUT, ">$output_prefix.merged.fa") || (print $log "Cannot open $output_prefix.merged.fa for writing" && exit 1);
	open(MAP, ">$output_prefix.merged.allele_map.txt") || (print $log "Cannot open $output_prefix.allele_map.txt for writing" && exit 1);
	
	#Iterates through sequence hash and prints out sequences and map file entry
	foreach my $seq (sort {$sequences{$a}{order}<=>$sequences{$b}{order}} keys %sequences){	
		$sequences{$seq}{name} =~ s/,$//;
		if($sequences{$seq}{cnt}>1){
			$sequences{$seq}{name}=~s/>//;
			my @names = split /,/, $sequences{$seq}{name};
			my @parts = split / /, $names[0];
			$parts[0]=~s/HLA://;
			$parts[0]=~s/IPD://;
			my $comment = "$parts[0]($parts[1]) sequence equivalent to ";	
		
			for (my $i=1; $i<=$#names;$i++){
				my @parts = split / /, $names[$i];
				$parts[0]=~s/HLA://;
				$parts[0]=~s/IPD://;
				$comment .= "$parts[0]($parts[1]);";
			}
			$comment =~ s/;$//;
			print OUT ">$sequences{$seq}{name}\n$seq\n";
			$parts[0]=~s/HLA://;
			$parts[0]=~s/IPD://;
			print MAP "$parts[0]\t$parts[1]\t$comment\t$sequences{$seq}{cwd}\n";
		}else{
			$sequences{$seq}{name} =~ s/>//;
			my @parts = split / /, $sequences{$seq}{name};
			my $comment = "";
			$parts[0]=~s/HLA://;
			$parts[0]=~s/IPD://;
			print OUT ">$sequences{$seq}{name}\n$seq\n";
			print MAP "$parts[0]\t$parts[1]\t$comment\t$sequences{$seq}{cwd}\n";
		}
	}
	
	close OUT;
	close MAP;
	close ($log) if ($log_is_file);
}
__PACKAGE__->runCommandline() unless caller;

