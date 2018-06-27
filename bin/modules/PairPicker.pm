#!/usr/bin/env perl
package PairPicker;
use strict;
use warnings;
use Getopt::Long;
use Module::Load;

(my $SCRIPT_NAME = $0) =~ s/.*\///;
my $version = "1.0";
my $creation_date = "1 Nov 2016";
my $last_updated = "14 Sep 2017";

my $usage = "\n$SCRIPT_NAME v$version\n" .
	    "\nDESCRIPTION:\n" .
	    "This module does pairwise comparisons of allele pairs to count the paired-end reads that are perfect matches to only one of the pairs in the comparison. This module can run as a standalone program or invoked from within another perl script.\n" .
	    "\nUSAGE\n" .
	    "perl $SCRIPT_NAME <options>\n" .
	    "\nRequired options:\n" .
	    "-input|i\t\tThe prefix of the paired-end fastq files.(required)\n" .
	    "-gene|g\t\t\tThe gene from which the alleles originate. (requred)\n" .
	    "-reference_file|r\tThe reference fasta containing the alleles. (required)\n" .
	    "-candidates|c\t\tTwo candidate alleles from the same gene separated by a pipe, i.e. HLA00001|HLA00002. (required)\n" .
	    "\nGeneral options:\n" .
	    "-output|o\t\tOutput file.(default:STDOUT)\n" .
	    "-min_quality|q\t\tThe contribution of base to the read count at a given position is weighted by a quality score Q which depends on base quality q." .
	    "\n\t\t\tWhen q < min_qual, Q=0; min_qual<=q<35, Q=(q - min_qual)/(35 - min_qual); q>35, Q=1 (default:20)\n" .
	    "-scale|sc\t\tThis controls the scale of reads found in both pairs of the comparison. For a scale of 1 the read is counted as Q (see -min_quality) while for a scale of 0 read is counted as 0.(default:0)\n" .
	    "-ambiguous_reads|a\tThis includes reads that map ambigiously among genes for the same class. For example, for Gene A this would also include reads classified as ClassI.(default:OFF)\n" .
	    "-scripts_dir|sd\t\tThe parent directory of the modules folder containing the HLAProfiler modules (default:'..')\n" .
	    "-help|h\t\t\tDisplays this message\n" .
            "\nAUTHORS:\n" . 
	    "Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	    "Chad Brown\n" .
	    "\nCREATED:\n$creation_date\n" .
	    "\nLAST UPDATED:\n$last_updated\n" .
	    "\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
	    "\n";

my %opts = (input=>"",scale=>0, reference_file=>"", gene=>"", min_quality=>20);
my %candidate_reference;
my @candidate_pairs;
my @candidate_counts;

sub runCommandline{
	GetOptions(\%opts, qw(input|i=s gene|g=s min_quality|q=s output|o=s candidates|c=s reference_file|r=s scripts_dir|sd=s scale|sc=s ambiguous|a=s help|h));

	if($opts{help}){
		print $usage;
		exit;
	}elsif ($opts{input} eq ""){
		print STDERR "Please specify an input.\n$usage\n";
		exit 1;
	}elsif ($opts{gene} eq ""){
		print STDERR "Please specify a gene.\n$usage\n";
		exit 1;
	}elsif ($opts{reference_file} eq ""){
		print STDERR "Please specify a reference sequence file.\n$usage\n";
		exit 1;
	}
	load "$opts{scripts_dir}/modules/SequenceFunctions.pm";
	printResults(pickPairs($opts{input}, $opts{reference_file}, $opts{candidates}, $opts{ambiguous}), $opts{output});
}

sub runModule{
	my $input = shift;
	$input = shift if($input eq "PairPicker");
	my $reference = shift;
	my $candidates = shift;
	$opts{gene} = shift;
	$opts{min_quality} = shift;
	$opts{scale} = shift;
	my $scripts_dir = shift;
	my $ambiguous = shift;
	
	%candidate_reference = ();
	@candidate_pairs = ();
	@candidate_counts = ();
	load "${scripts_dir}/modules/SequenceFunctions.pm";
	return(pickPairs($input, $reference, $candidates, $ambiguous));
}

sub pickPairs{
	my $input = shift;
	$input = shift if($input eq "PairPicker");
	my $reference = shift;
	my $candidate_list = shift;
	my $include_ambiguous = shift || 0;
		
	my @candidates = split /[,\|]/, $candidate_list;
	@candidate_pairs = split /[,]/, $candidate_list;
	
	if($#candidates == -1){
		print STDERR "Candidate list is empty!\n";
		exit 1;
	}
	
	##Create empty count matrix
	for (my $i=0; $i<=$#candidate_pairs; $i++){
		push(@candidate_counts, [(0)x($#candidate_pairs +1)]);
	}	
	
	##Create a hash of candidates
	my %candidate_hash;
	foreach my $candidate (@candidates){
		$candidate_hash{$candidate}=0;
	}
	
	##Read in reference and only store the sequence for candidates
	open(IN, $reference) || (print STDERR  "Cannot open reference file: $reference\n" and exit 1);
	while(<IN>){
		my $line1 = $_;
		my $line2 = <IN>;
		chomp $line1;
		chomp $line2;
	
		my @parts = split /[\t ]/, $line1;
		my $allele = substr($parts[0],1);
		$allele=~s/HLA\://;
		$allele=~s/IPD\://;
		#Allows for either the accession or allele name to be used in candidate list
		if (defined($candidate_hash{$allele})){
			$candidate_reference{$allele}=$line2;
		}elsif (defined($candidate_hash{$parts[1]})){
			$candidate_reference{$parts[1]}=$line2;	
		}
	}	
	
	
	my $unique_read_count = 0;
	$unique_read_count += readFile("$input.$opts{gene}");

	## This allows for the reads shared between genes to be included to help in cases where homology between genes is reducing the power of PairPicker to resolve between alleles.
	if($include_ambiguous){
		my %class = (A=>"ClassI", B=>"ClassI", C=>"ClassI", DMA=>"ClassII", DMB=>"ClassII", DOA=>"ClassII", DOB=>"ClassII", DPA1=>"ClassII",DPA2=>"ClassII",DPA2=>"ClassII", DPB1=>"ClassII", DPB2=>"ClassII", DQA1=>"ClassII", DQB1=>"ClassII", DQB2=>"ClassII", DQB3=>"ClassII", DRA=>"ClassII", DRB1=>"ClassII", DRB2=>"ClassII", DRB3=>"ClassII", DRB4=>"ClassII", DRB5=>"ClassII", DRB6=>"ClassII", DRB7=>"ClassII", DRB8=>"ClassII", DRB9=>"ClassII", E=>"ClassI", F=>"ClassI", G=>"ClassI", H=>"ClassI", HFE=>"Non-HLA", J=>"ClassI", K=>"ClassI", L=>"ClassI", MICA=>"Non-HLA", MICB=>"Non-HLA", TAP1=>"Non-HLA", TAP2=>"Non-HLA", V=>"ClassI", N=>"ClassI", P=>"ClassI",S=>"ClassI",T=>"ClassI",U=>"ClassI",W=>"ClassI",X=>"ClassI",Y=>"ClassI",Z=>"ClassI", PSMB9=>"Non-HLA", PSMB8=>"Non-HLA",MICC=>"Non-HLA",MICD=>"Non-HLA",MICE=>"Non-HLA");
		$unique_read_count += readFile("$input.$class{$opts{gene}}");
	}
	return(\@candidate_counts,\@candidate_pairs, $unique_read_count);
}
	
sub readFile{
	my $input = shift;
	
	my $fh1;
	my $fh2;
	if (-e "${input}_1.fastq.gz"){

		open($fh1, "<", "gunzip -c ${input}_1.fastq.gz |") || (print STDERR "Cannot open file ${input}_1.fastq.gz\n" and exit 1);
		open($fh2, "<", "gunzip -c ${input}_2.fastq.gz |") || (print STDERR "Cannot open file ${input}_2.fastq.gz\n" and exit 1);

	}elsif (-e "${input}_1.fastq"){
		open($fh1, "<", "${input}_1.fastq") || (print STDERR "Cannot open file ${input}_1.fastq\n" and exit 1);
		open($fh2, "<", "${input}_2.fastq") || (print STDERR "Cannot open file ${input}_2.fastq\n" and exit 1);
	}elsif (-e "${input}_1.fq.gz"){
		open($fh1, "<", "gunzip -c ${input}_1.fq.gz |") || (print STDERR "Cannot open file ${input}_1.fq.gz\n" and exit 1);
		open($fh2, "<", "gunzip -c ${input}_2.fq.gz |") || (print STDERR "Cannot open file ${input}_2.fq.gz\n" and exit 1);
	}elsif (-e "${input}_1.fq"){
		open($fh1, "<", "${input}_1.fq") || (print STDERR "Cannot open file ${input}_1.fq\n" and exit 1);
		open($fh2, "<", "${input}_2.fq") || (print STDERR "Cannot open file ${input}_2.fq\n" and exit 1);
	}else{
		print STDERR "Cannot find fastq file with prefix ${input}\n";
		exit 1;
	}
	my $unique_read_count = 0;	
	while(<$fh1>){
		my $header1=$_;
		my $sequence1=<$fh1>;
		my $filler1=<$fh1>;
		my $quality1=<$fh1>;
		
		my $header2=<$fh2>;
		my $sequence2=<$fh2>;
		my $filler2=<$fh2>;
		my $quality2=<$fh2>;
		chomp($sequence1);
		chomp($sequence2);
		chomp($quality1);
		chomp($quality2);	
	
		my @alleles = ();
		my $allele_cnt = 0;
		my %containsRead =();
		my $pair_count=0;
		my $quality_read = 1;
		
		my @qual1  = split //, $quality1;
		my @qual2  = split //, $quality2;
		
		my $min_quality = 100;
		##The read is only considered if all bases are above the minimum quality
		if($#qual1 == $#qual2){
			QUAL:for(my $i=0; $i<=$#qual1; $i++){
				my $q1 = ord($qual1[$i]) -33;
				my $q2 = ord($qual2[$i]) -33;
				if ($q1 < $opts{min_quality} || $q2<$opts{min_quality}){
					$quality_read = 0;
					last QUAL;
				}else{
					$min_quality = $q1 if ($q1 <$min_quality);
					$min_quality = $q2 if ($q2 <$min_quality);
				}
			}
		}else{
			my $max_index = $#qual1;
			$max_index = $#qual2 if ($#qual1 < $#qual2);
			QUAL:for(my $i=0; $i<=$max_index; $i++){
				my $q1 = "";
				my $q2 = ""; 
				$q1 =  ord($qual1[$i]) -33 if ($i <= $#qual1);
				$q2 =  ord($qual2[$i]) -33 if ($i <= $#qual2);
				
				if (($i<=$#qual1 && $q1 < $opts{min_quality}) ||($i<=$#qual2 && $q2< $opts{min_quality})){
					$quality_read = 0;
					last QUAL;
				}else{
					$min_quality = $q1 if ($q1 <$min_quality);
					$min_quality = $q2 if ($q2 <$min_quality);	
				}
			}
		}


		if($quality_read == 1){
			my $qual_score = 1;
			#Weights the contribution of the reads based on the lowest base quality. 
			if($min_quality < 35){
				my $diff = 35 - $opts{min_quality};
				$qual_score = ($min_quality-$opts{min_quality})/$diff;
			}
			#Find pairs containing reads
			foreach my $pair (@candidate_pairs){
				my @pair_parts = split /\|/, $pair;
				my $candidate1 = $pair_parts[0];
				my $candidate2 = $pair_parts[1];
				my $in_pair = 0;
				if($candidate_reference{$candidate1}=~m/$sequence1/){
					my $sequence2r = SequenceFunctions->revcomp($sequence2);
					if($candidate_reference{$candidate1}=~m/$sequence2r/){
						$in_pair = 1;
					}
				}
				if($candidate_reference{$candidate1}=~m/$sequence2/){
					my $sequence1r = SequenceFunctions->revcomp($sequence1);
					if($candidate_reference{$candidate1}=~m/$sequence1r/){
						$in_pair = 1;
					}
				}
				if($candidate_reference{$candidate2}=~m/$sequence1/){
					my $sequence2r = SequenceFunctions->revcomp($sequence2);
					if($candidate_reference{$candidate2}=~m/$sequence2r/){
						$in_pair = 1;
					}
				}
				if($candidate_reference{$candidate2}=~m/$sequence2/){
					my $sequence1r = SequenceFunctions->revcomp($sequence1);
					if($candidate_reference{$candidate2}=~m/$sequence1r/){
						$in_pair = 1;
					}
				}	
				$containsRead{$pair}=$in_pair;
				$pair_count += $in_pair;
			}
			if($pair_count != 0){
				my $unique_pair = 0;
				#Do pairwise comparisons of all candidate pairs
				for (my $i=0; $i<=$#candidate_pairs;$i++){
					for (my $j=0; $j<=$#candidate_pairs;$j++){
						if($i<$j){
							my $p1 = $containsRead{$candidate_pairs[$i]};
							my $p2 = $containsRead{$candidate_pairs[$j]};
							if($p1 == 1 && $p2 == 1){
								$candidate_counts[$i][$j]+=(1*$opts{scale}*$qual_score);		
								$candidate_counts[$j][$i]+=(1*$opts{scale}*$qual_score);		
							}elsif($p1 == 1 && $p2 == 0){
								$candidate_counts[$i][$j]+=$qual_score;		
								$unique_pair = 1;	
							}elsif($p1 == 0 && $p2 == 1){
								$candidate_counts[$j][$i]+=$qual_score;		
								$unique_pair = 1;	
							}
						}
					
					}
				}
				$unique_read_count += ($unique_pair*$qual_score);
			}
		}
	}
	return($unique_read_count);
}

sub printResults{
	my $cc_ref = shift;
	my @candidate_counts = @$cc_ref;
	my $cp_ref = shift;
	my @candidate_pairs = @$cp_ref;
	my $unique_read_count = shift;
	my $output = shift;
	my $oh = *STDOUT;
	if($output ne ""){
		open($oh, ">", "$output") || (print STDERR "Cannot open output file '$output'\n" && exit 1);	
	}
	for (my $i=0; $i<=$#candidate_pairs;$i++){
		print $oh "\t$candidate_pairs[$i]";
	}
	print $oh "\n";
	for (my $i=0; $i<=$#candidate_pairs;$i++){
		print $oh "$candidate_pairs[$i]";
		for (my $j=0; $j<=$#candidate_pairs;$j++){
			print $oh "\t$candidate_counts[$i][$j]";
		}
		print $oh "\n";
	}
}

__PACKAGE__->runCommandline() unless caller;
