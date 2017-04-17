#!/usr/bin/env perl
package RunKraken;
use strict;
use warnings;
use Getopt::Long;

(my $SCRIPT_NAME = $0) =~ s/.*\///;
my $version = "1.0";
my $creation_date = "1 Oct 2016";
my $last_updated = "12 Jan 2017";

my $usage="\n$SCRIPT_NAME v$version\n".
	  "\nDESCRIPTIONs:\n" .
	  "This module constructs and executes the system call to Kraken. The module can be run from commandline or called directly from another perl script.\n" . 
	  "\nUSAGE:\n" .
	  "perl $SCRIPT_NAME <options>\n" .
	  "\nRequired options:\n" .
	  "-mode|m\t\t\tMode of Kraken operation (required)\n" .
	  "\t\t\t   Modes available:\n" . 
	  "\t\t\t\tbuild\tBuilds a Kraken database\n" . 
	  "\t\t\t\tfilter\tFilters reads using an existing Kraken database\n" . 
	  "\t\t\t\tlibrary\tAdds the input fasta to the Kraken database\n" . 
	  "-database_name|db\tName of the kraken_database (required)\n" .
	  "\nBuild specific operations:\n" . 
	  "-build_options|b\tA list of build parameters for the Kraken database (default:'')\n" .
	  "\nFilter specific operations:\n" . 
	  "-fastq1|fq1\t\t\tFastq file read1 (required)\n" .
	  "-fastq2|fq2\t\t\tFastq file read2 (required)\n" .
	  "-output_prefix|op\t\tOutput prefix (required)\n" .
	  "-print_classifications|p\tIf a filename is specified will print classifications to that file. (default:\"\")\n" .
	  "\nLibrary specific operations:\n" . 
	  "-fasta|fa\tFasta file to add to library(required)\n" .
	  "\nGeneral options:\n" .
	  "-database_directory|dd\tThe directory containing the Kraken database(default:'.')\n" .
	  "-threads|c\t\tNumber of threads to use(default:1)\n" .
	  "-kraken_path|c\t\tThe location of the Kraken executables (default:'.')\n" .
	  "-log|l\t\t\tLog file(default:STDERR)\n" .
          "-help|h\t\t\tDisplays this message\n" .
          "\nAUTHORS:\n" .
	  "Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	  "Chad Brown:chad.brown\@q2labsolutions.com\n" .
	  "\nCREATED:\n$creation_date\n" .
	  "\nLAST UPDATED:\n$last_updated\n" .
	  "\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
	  "\n";

sub runCommandline{
	my %opts = (mode=>"", database_directory=>".", threads=>1, kraken_path=>".", build_options=>"", fasta=>"", fastq1=>"", fastq2=>"", database_name=>"", output_prefix=>"", print=>"");
	GetOptions(\%opts, qw(help|h mode|m=s database_name|db=s database_directory|dd=s print_classifications|p=s fastq1|fq1=s fastq2|fq2=s output_prefix|op=s threads|c=i kraken_path|kp=s build_options|b=s fasta|fa=s log|l=s));
	if($opts{help}){
		print $usage;
		exit;
	}elsif($opts{mode} eq ""){
		print STDERR "Please specify an operation mode. $usage";
		exit 1;
	}elsif($opts{database_name} eq ""){
		print STDERR "Please specify a database name. $usage";
		exit 1;
	}else{
		if($opts{mode} eq "filter"){
			if($opts{fastq1} eq ""){
				print STDERR "Please specify a read 1 fastq files. $usage"; 
				exit 1;
			}elsif($opts{fastq2} eq ""){
				print STDERR "Please specify a read 2 fastq files. $usage"; 
				exit 1;
			}elsif($opts{output_prefix} eq ""){
				print STDERR "Please specify and output prefix. $usage";
				exit 1;
			}else{
				filterReads($opts{database_directory},$opts{database_name},$opts{output_prefix}, $opts{threads},$opts{fastq1},$opts{fastq2},$opts{kraken_path},$opts{print_classifications}, $opts{log});
			}
		}	
		if($opts{mode} eq "build"){
			buildKraken($opts{database_directory},$opts{database_name},$opts{build_options}, $opts{fasta},$opts{threads},$opts{kraken_path},$opts{log});
		}	
		if($opts{mode} eq "library"){
			if($opts{fastq1} eq ""){
				print STDERR "Please specify a fasta files to add to library. $usage"; 
				exit 1;
			}else{
				addLibraryToKraken($opts{database_directory},$opts{database_name},$opts{fasta}, $opts{kraken_path},$opts{log});
			}
		}
	}	
}


sub buildKraken{
	my $database_dir = shift;
	$database_dir = shift if ($database_dir eq "RunKraken");
	my $database_name = shift;
	my $options = shift;
	my $threads = shift;
	my $path = shift; 
	my $logfile = shift;
	
        my $log;
	my $is_log_file = 0;
        if($logfile){
                open($log, ">>$logfile");
		$is_log_file = 1;
        }else{
                $log = *STDERR;
        }

	##Remove any " or ' surrounding the options
	$options=~s/^"//;	
	$options=~s/"$//;	
	$options=~s/^'//;	
	$options=~s/'$//;	
	print $log "$path/kraken-build --build $options --threads=$threads --db $database_dir/$database_name\n";
	open (IN, "$path/kraken-build --build $options --threads=$threads --db $database_dir/$database_name 2>&1|");
	while(<IN>){
		print $log $_;
	}
	close ($log) if ($is_log_file == 1);
	return($?);
}

sub addLibraryToKraken{
	my $database_dir = shift;
	$database_dir = shift if ($database_dir eq "RunKraken");
	my $database_name = shift;
	my $file = shift;
	my $path = shift; 
	my $logfile = shift;
	
        my $log;
	my $is_log_file = 0;
        if($logfile){
                open($log, ">>$logfile");
		$is_log_file = 1;
        }else{
                $log = *STDERR;
        }
	
	print $log "$path/kraken-build --add-to-library $file --db $database_dir/$database_name\n";
	open (IN, "$path/kraken-build --add-to-library $file --db $database_dir/$database_name 2>&1 |");
	while(<IN>){
		print $log $_;
	}
	close ($log) if ($is_log_file == 1);
	return($?);
}

sub filterReads{
	my $database_dir = shift;
	$database_dir = shift if ($database_dir eq "RunKraken");
	my $database_name = shift;
	my $filtered_output = shift;
	my $threads = shift;
	my $fq1 = shift;
	my $fq2 = shift;
	my $path = shift;
	my $base_dir  = "$database_dir/$database_name";
	my $print = shift;
	my $logfile = shift;
	
        my $log;
	my $is_log_file = 0;
        if($logfile){
                open($log, ">>$logfile");
		$is_log_file = 1;
        }else{
                $log = *STDERR;
        }
	
	my $suppress_flag = "-S";
	$suppress_flag = "" if ($print ne "");	
	print $log "$path/classify -t $threads -f -x ${filtered_output}. -p $base_dir/taxonomy/divisions.txt -d $base_dir/database.kdb -n $base_dir/taxonomy/nodes.dmp -i $base_dir/database.idx $suppress_flag -2 $fq1 $fq2 2>>$logfile\n";
	my $output = `$path/classify -t $threads -f -x ${filtered_output}. -p $base_dir/taxonomy/divisions.txt -d $base_dir/database.kdb -n $base_dir/taxonomy/nodes.dmp -i $base_dir/database.idx $suppress_flag -2 $fq1 $fq2 2>>$logfile`;
	
	my $ofh;
	my $output_is_file = 0;
	if ($print ne ""){
		open($ofh,">$print");	
		$output_is_file = 1;
	}else{
		$ofh = *STDOUT;	
	}
	print $ofh $output;
	close ($ofh) if ($output_is_file==1);
	close ($log) if ($is_log_file == 1);
	return($?)
}
__PACKAGE__->runCommandline() unless caller;
