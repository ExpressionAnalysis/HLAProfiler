#!/usr/bin/env perl
package ReadCounter;
use strict;
use warnings;
use Getopt::Long;

(my $SCRIPT_NAME = $0) =~ s/.*\///;
my $version = "1.0";
my $creation_date = "6 Mar 2011";
my $last_updated = "14 Sep 2017";

my $usage="\n$SCRIPT_NAME v$version\n".
	  "\nDESCRIPTIONs:\n" .
	  "This module counts the number of occurences of each sequence in a fastq file and stores the counts to file. The module can be run from commandline or called directly from another perl script.\n" . 
	  "\nUSAGE:\n" .
	  "perl $SCRIPT_NAME <options>\n" .
	  "\nRequired options:\n" .
	  "-input_fq|i\tFastq file (required)\n" .
	  "-output_file|o\tName of output file (required)\n" .
	  "\nGeneral options:\n" .
	  "-log|l\tLog file(default:STDERR)\n" .
          "-help|h\tDisplays this message\n" .
          "\nAUTHORS:\n" .
	  "Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	  "Chad Brown\n" .
	  "Erik Aronesty\n" .
	  "\nCREATED:\n$creation_date\n" .
	  "\nLAST UPDATED:\n$last_updated\n" .
	  "\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
	  "\n";

sub runCommandline{
	my %opts = (input_fq=>"", output_file=>"");
	GetOptions(\%opts, qw(help|h input_fq|i=s output_file|o=s log|l=s));
	if($opts{help}){
		print $usage;
		exit;
	}elsif($opts{input_fq} eq ""){
		print STDERR "Please specify a fastq file. $usage";
		exit 1;
	}elsif($opts{output_file} eq ""){
		print STDERR "Please specify an output_file. $usage";
		exit 1;
	}else{
		count_reads($opts{input_fq}, $opts{output_file}, $opts{log});
	}
	

}


sub count_reads{
	my $infile = shift;
	if($infile eq "ReadCounter"){
		$infile = shift;
	}
	my $outfile = shift;
	my $logfile = shift;
	my $log;
	my $log_is_file = 0; 
	if($logfile){
		open($log, ">>$logfile");
		$log_is_file = 1;
	}else{
		$log = *STDERR;
	}

		
	open(INFILE, $infile =~ /gz$/ ? "gunzip -c $infile |" : $infile) or (print $log "Unable to open $infile for reading\n" and exit 1);
	open (OUTFILE, ">$outfile")  or (print $log "Unable to open $outfile for writing\n" and exit 1);
	
	my %unique_hash = ();
	my $num_total_reads = 0;
	my $uniq_count = 0;
	my $line_count = 0;
	my $uniq_fraction = 0;
	my $go = 0;
	
	while( my $line = <INFILE> ) {
		$line=<INFILE>;
		chomp $line;
		$unique_hash{$line}++;	##The value of each key is the number of times that sequence was hit
		scalar <INFILE>;
		scalar <INFILE>;
		$line_count++;
	}
	
	my @temp = keys(%unique_hash);
	
	@temp = sort {$unique_hash{$b} <=> $unique_hash{$a} || $a cmp $b} @temp;
	
	for (@temp) {
		print OUTFILE "$_\t$unique_hash{$_}\n";
	}
	
	$uniq_count = @temp;
	$uniq_fraction = $uniq_count/$line_count;
	
	close (INFILE);
	close (OUTFILE);
	close ($log) if $log_is_file;
}
__PACKAGE__->runCommandline() unless caller;
