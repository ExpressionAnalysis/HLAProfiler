#!/usr/bin/env perl
package SimulateReads;
use strict;
use warnings;
use Getopt::Long;
use List::Util;
use Parallel::ForkManager;
use Module::Load;

(my $SCRIPT_NAME = $0) =~ s/.*\///;
my $version = "1.0";
my $creation_date = "1 Oct 2016";
my $last_updated = "12 Jan 2017";

my $usage="\n$SCRIPT_NAME v$version\n".
	  "\nDESCRIPTIONs:\n" .
	  "This module simulates reads from each sequence in a fasta file or from a single sequence. The reads do not contain errors and assume a stranded protcol. A Pareto distribution is used to determine insert size of each read. The module can be run from commandline or called directly from another perl script.\n" . 
	  "\nUSAGE:\n" .
	  "perl $SCRIPT_NAME <options>\n" .
	  "\nRequired options:\n" .
	  "-reference|r\tReference fasta for simulation (required)\n" .
	  "\t-OR-\n" .
	  "-sequence|sq\tThe reference sequence for simulation (required)\n" .
	  "-name|fn\tThe name of the reference sequence (required)\n" .
	  "\nSimulation Options:\n" . 
	  "-numReads|n\tNumber of reads to simulate per sequence (default:500000)\n" .
	  "-readLength|l\tThe length of the reads to simulate (default:50)\n" .
	  "-outPrefix|o\tThe prefix of the output files (default:\"simulatedReads\")\n" .
	  "-max_insert|m\tThe maximum insert length (default:1000)\n" .
	  "-seed|se\tThe seed for random sampling (default:1234)\n" .
	  "-scale|sc\tThe scale factor for the Pareto distribution estimating insert size (default:80)\n" .
	  "-shape|sh\tThe shape factor for the Pareto distribution estimating insert size (default:0.7\n" .
	  "\nGeneral options:\n" .
	  "-scripts_dir|sd\t\tParent directory containing modules directory with HLAProfiler modules (default:'.')\n" .
	  "-threads|c\t\tNumber of threads to use(default:1)\n" .
	  "-log|l\t\t\tLog file(default:STDERR)\n" .
          "-help|h\t\t\tDisplays this message\n" .
          "\nAUTHORS:\n" .
	  "Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	  "Chad Brown:chad.brown\@q2labsolutions.com\n" .
	  "\nCREATED:\n$creation_date\n" .
	  "\nLAST UPDATED:\n$last_updated\n" .
	  "\nCopyright. Q2 Solutions|EA Genomics. 2016\n" ;

my %opts = (name=>"", sequence=>"", threads=>1, reference=>"", numReads=>500000, readLength=>50, outPrefix=>"simulatedReads", max_insert=>1000, seed=>1234, scale=>80, shape=>0.7, scripts_dir=>'.');

sub runCommandline{
	my $package = shift;
	GetOptions(\%opts, qw(threads|t=s sequence|sq=s  name|fn=s reference|r=s numReads|n=s readLength|l=s outPrefix|o=s max_insert|m=s seed|se=s scale|sc=s shape|sh=s log|l=s scripts_dir|sd=s help|h));
	
	if($opts{help}){
		print $usage;
		exit;
	}elsif($opts{reference} ne ""){
		load ("$opts{scripts_dir}/modules/SequenceFunctions.pm");
		$opts{shape} = 1.0/$opts{shape};
		simulateReads();
	}elsif($opts{sequence} ne "" && $opts{name} ne ""){
		load ("$opts{scripts_dir}/modules/SequenceFunctions.pm");
		$opts{shape} = 1.0/$opts{shape};
		makeSim($opts{name}, $opts{sequence},$opts{outPrefix},$opts{log});
	}else{
		print STDERR "Please specify a file of reference sequences OR a reference sequence and name.\n$usage\n";
		exit 1;
	}
}


##This allows another script to set global variables and then make repeated calls to makeSim
sub setSimulationOptions{
	$opts{reference}=shift;
	$opts{reference}=shift if ($opts{reference} eq "SimulateReads");
	$opts{numReads}=shift;
	$opts{readLength}=shift;
	$opts{max_insert}=shift;
	$opts{scale}=shift;
	$opts{shape}=shift;
	$opts{shape} = 1.0/$opts{shape};
	$opts{seed}=shift;
	$opts{threads}=shift;
	$opts{scripts_dir}=shift;
	$opts{outPrefix}=shift || "simulatedReads";
	load ("$opts{scripts_dir}/modules/SequenceFunctions.pm");
	return(\%opts);
}
	
## Simulate the reads
sub makeSim {
	##Set the seed for rand
	srand $opts{seed};

	## sequence name
	my $fn = shift;
	if($fn eq "SimulateReads"){
		$fn = shift;
	}

	## sequence
	my $fastaString = shift;
	my $out_prefix = shift;
	my $logfile = shift;
	
	my $log;
	my $is_log_file = 0; 
	if($logfile){
		open($log, ">>$logfile");
		$is_log_file=1;
	}else{
		$log = *STDERR;
	}

	my $r2s; my $r1s;
		
	## Remove standard prefixes from name 
	$fn =~ s/>HLA://;
	$fn =~ s/>IPD://;
	## Remove everything after the first space
	$fn =~ s/ .+//;

	open(OUT1, ">$out_prefix.${fn}_1.fastq") || (print $log "Could not open $out_prefix.${fn}_1.fastq for writing\n" && exit 1);
	open(OUT2, ">$out_prefix.${fn}_2.fastq") || (print $log "Could not open $out_prefix.${fn}_2.fastq for writing\n" && exit 1);;
	
	for(my $i=0; $i < $opts{numReads}; $i++) {
		my $length = length($fastaString);
		my $valid = 0;
		my $cnt = 0;
		my $curr_max_insert = $opts{max_insert};
		while(!$valid){ ##Iterates until a valid simulation is created
			my $inserts = $opts{scale} / (rand(1)**$opts{shape});
			$inserts = $inserts < $curr_max_insert ? $inserts : $curr_max_insert; ## Sets insert size to $curr_max_insert if it is larger than $curr_max_insert
			$inserts = int($inserts + 0.5);
			$r2s = int( (rand(1)*$length) + 0.5);
			$r1s = $r2s+$opts{readLength}-1+$inserts;
			if( ($r1s+$opts{readLength}) <= $length){
				$valid = 1;
			} elsif($cnt > 50){
				$curr_max_insert = $curr_max_insert*0.5; ##This splits the maximum insert size in half after 100 invalid tries to prevent and infinite loop for smaller reference
				$cnt = 0;
			} else {
				$cnt = $cnt+1;
			}
		}
		my $readName = "\@Sim_".$i."_R1_".$r1s."_R2_".$r2s;
		my $read1 = uc substr($fastaString, $r1s, $opts{readLength});
		$read1 = SequenceFunctions->revcomp($read1);
		my $read2 = uc substr($fastaString, $r2s, $opts{readLength});
		my $qual = "I" x $opts{readLength}; ## Gives all bases qual of 40
		print OUT1 "$readName\n$read1\n+\n$qual\n";
		print OUT2 "$readName\n$read2\n+\n$qual\n";			
	}
	close(OUT1);
	close(OUT2);
	close ($log) if ($is_log_file == 1);
	return($fn);
}

sub simulateReads{
	my $dh;
	open($dh, $opts{reference}) || (print STDERR "Cannot open reference file: $opts{reference}\n" && exit 1);
	my $line  = <$dh>;
	chomp $line;
	my $fn = $line;
	my $fastaString = "";
	my $fm = Parallel::ForkManager->new($opts{threads});
	FASTA:while($line  = <$dh>){
		chomp $line;
		# need lines 30 - 60 in a subroutine
		if(substr($line, 0, 1) eq ">") {
			my $prev_fn = $fn;
			my $prev_fastaString = $fastaString;
			$fn = $line;
			$fastaString = "";
			#my $pid = $fm->start and next FASTA;
				makeSim($prev_fn, $prev_fastaString, $opts{outPrefix}, $opts{log});
			#$fm->finish;
		} else {
			$fastaString=$fastaString.$line;
		}
	}
	#$fm->wait_all_children;
	makeSim($fn, $fastaString, $opts{outPrefix}, $opts{log});
}

__PACKAGE__->runCommandline() unless caller;
