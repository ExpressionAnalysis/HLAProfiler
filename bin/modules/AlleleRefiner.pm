#!/usr/bin/env perl
package AlleleRefiner;
use strict;
use warnings;
use Getopt::Long;
use Module::Load;

(my $SCRIPT_NAME = $0) =~ s/.*\///;
my $version = "1.0";
my $creation_date = "1 Nov 2016";
my $last_updated = "18 Jan 2018";

my $usage = "\n$SCRIPT_NAME v$version\n" .
	    "\nDESCRIPTION:\n" .
	    "A perl module created to interface with HLAPredict.pm. This module aligns reads from the specified paired-end fastq file to the reference sequences of a predicted alleles pair for a gene from the HLA database and predicts any variations between the reference sequence and the observed reads. This module can run as a standalone program or invoked from within another perl script.\n" .
	    "\nUSAGE\n" .
	    "perl $SCRIPT_NAME <options>\n" .
	    "\nRequired options:\n" .
	    "-input|i\t\tThe prefix of the paired-end fastq files.(required)\n" .
	    "-output|o\t\tOutput prefix.(required)\n" .
	    "-gene|g\t\t\tThe gene from which the alleles originate. (requred)\n" .
	    "-reference_file|r\tThe reference fasta containing the alleles. (required)\n" .
	    "-candidates|c\t\tTwo candidate alleles from the same gene separated by a pipe, i.e. HLA00001|HLA00002. (required)\n" .
	    "\nGeneral options:\n" .
	    "-min_qual|q\t\tThe contribution of base to the read count at a given position is weighted by a quality score Q which depends on base quality q." .
	    "\n\t\t\tWhen q < min_qual, Q=0; min_qual<=q<35, Q=(q - min_qual)/(35 - min_qual); q>35, Q=1 (default:20)\n" .
	    "-counts_threshold|t\tThe minimum number of non-reference counts required for a position to be considered a mismatch (default:10)\n" .
	    "-kmer|l\t\t\tThe size of the seed kmer for alignment (default:20)\n" .
	    "-mismatch|m\t\tThe maximum number of mismaatches allowed in a single read during alignment (default:1)\n" .
	    "-minimum_frequency|f\t\tThe minimum frequency, calculated as [non-reference count/(non-reference count+ reference count)] required for a position to be considered a reference mismatch (default:.75)\n" .
	    "-scripts_dir|sd\t\tThe parent directory of the modules folder containing the HLAProfiler modules (default:'..')\n" .
	    "-help|h\t\t\tDisplays this message\n" .
	    "\nOutput options\n" .
	    "-print_coverage|pc\tPrints the coverage for each position in each allele to file (default:off)\n" .
	    "-print_mismatches|pm\tPrints the bases with possible mismatches to file\n" .
	    "-print_kmers|pk\t\tPrints the sequence of read pairs that do not match perfectly to the reference\n" .
            "\nAUTHORS:\n" . 
	    "Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	    "Chad Brown\n" .
	    "\nCREATED:\n$creation_date\n" .
	    "\nLAST UPDATED:\n$last_updated\n" .
	    "\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
	    "\n";

my %opts = (gene=>"",ouput=>"", input=>"", min_qual=>20, reference_file=>"", minimum_frequency=>0.75,counts_threshold=>10, debug=>0,print_coverage=>0,print_mismatches=>0,print_kmers=>0, kmer=>20, mismatch=>1, scripts_dir=>"..", help=>0);

my %allele_sizes = ();
my %reference_hash =();
my %candidate_reference = ();
my %gene_references = ();
my %extra_kmers;
my %mismatches = ();
my %reference_counts = ();
my @candidates = ();

sub runCommandline{
	GetOptions(\%opts, qw(gene|g=s input|i=s output|o=s min_qual|q=s candidates|c=s reference_file|r=s counts_threshold|t=s kmer|k=s print_coverage|pc print_mismatches|pm print_kmers|pk mismatch|m=s minimum_frequency|f=s scripts_dir|sd=s help));
	
	
	if($opts{help}){
		print "$usage";
		exit;
	}elsif ($opts{input} eq ""){
		print STDERR "Please specify an input.\n$usage\n";
		exit 1;
	}elsif ($opts{output} eq ""){
		print STDERR "Please specify an output prefix.\n$usage\n";
		exit 1;
	}elsif ($opts{gene} eq ""){
		print STDERR "Please specify a gene.\n$usage\n";
		exit 1;
	}elsif ($opts{reference_file} eq ""){
		print STDERR "Please specify a reference sequence file.\n$usage\n";
		exit 1;
	}elsif ($opts{candidates} eq ""){
		print STDERR "Please specify two candidates alleles from the same gene.\n$usage\n";
		exit 1;
	}
	load "$opts{scripts_dir}/modules/SequenceFunctions.pm";
	my ($allele_ref, $cov_ref) = detectMutations();
	my %alleles = %{$allele_ref};
	my %coverage = %{$cov_ref};
	foreach my $allele (sort {$a cmp $b} keys %coverage){
		print "$allele\t$coverage{$allele}{mcon}\t$coverage{$allele}{pcon}\t$coverage{$allele}{mcov}\t$coverage{$allele}{pcov}\n"
	}

	foreach my $allele (sort {$a cmp $b} keys  %alleles){
		if($alleles{$allele}{replace}==1){
			print "Based on the observed data $allele is really predicted to be $alleles{$allele}{allele} $alleles{$allele}{name}.\n>$alleles{$allele}{allele} $alleles{$allele}{name} " . length($alleles{$allele}{seq}) . " bp\n$alleles{$allele}{seq}\n";
		}
	}
}

sub runModule{
	$opts{gene}= shift;
	$opts{gene}=shift if ($opts{gene} eq "AlleleRefiner");
	$opts{input}=shift;
	$opts{candidates}=shift;
	$opts{reference_file}=shift;
	$opts{min_qual}=shift;
	$opts{counts_threshold}=shift;
	$opts{kmer}=shift;
	$opts{mismatch}=shift;
	$opts{minimum_frequency}=shift;	
	$opts{scripts_dir}=shift;
	$opts{output} = shift;
	$opts{print_mismatches} = shift;
	$opts{print_kmers} = shift;
	$opts{print_coverage} = shift;
	
	return detectMutations();
}

sub detectMutations{
	my @inputs = split/,/, $opts{input};
	@candidates = split /[,\|]/, $opts{candidates};
	my @candidate_pairs = split /[,]/, $opts{candidates};
	print STDERR "Candidate list is empty!\n" and exit 1 if ($#candidates<0);
	
	my %candidate_counts;
	foreach my $candidate (@candidate_pairs){
		$candidate_counts{$candidate}=0;
	}
	
	my %candidate_hash;
	foreach my $candidate (@candidates){
		$candidate_hash{$candidate}=0;
	}
	open(IN, $opts{reference_file}) || (print STDERR "Cannot open reference file: $opts{reference_file}\n" && exit 1);
	
	#Read in each reference sequence and creates a hash with the position of all k-mers from each allele
	while(<IN>){
		my $line1 = $_;
		my $line2 = <IN>;
		chomp $line1;
		chomp $line2;
	
		my @parts = split /[\t ]/, $line1;
		my $allele = substr($parts[0],1);
		$allele=~s/HLA\://;
		$allele=~s/IPD\://;
		$parts[1]=~m/([A-Z0-9a-z]+)\*/;
		my $gene = $1;	
		if($gene eq $opts{gene}){
			$gene_references{$allele}{name}=$parts[1];
			$gene_references{$allele}{seq}=$line2;
		}	
		#This allows for either accession or allele name to be used as a candidate. HLAProfiler uses accession ids only.
		if (defined($candidate_hash{$allele})){
			$allele_sizes{$allele}=$parts[2];
			addAlleleToReference($allele,$line2,$opts{kmer},\%reference_hash,\%candidate_reference)
		}elsif (defined($candidate_hash{$parts[1]})){
			$allele_sizes{$parts[1]}=$parts[2];
			addAlleleToReference($parts[1],$line2,$opts{kmer},\%reference_hash,\%candidate_reference)
		}
	}	
	
	foreach my $input (@inputs){
		my $fh1;
		my $fh2;
		if (-e "${input}_1.fastq.gz"){
			open($fh1, "<", "gunzip -c ${input}_1.fastq.gz |") || (print STDERR "Cannot open file ${input}_1.fastq.gz\n" && exit 1);
			open($fh2, "<", "gunzip -c ${input}_2.fastq.gz |") || (print $STDERR "Cannot open file ${input}_2.fastq.gz\n" & exit 1);
		}elsif (-e "${input}_1.fastq"){
			open($fh1, "<", "${input}_1.fastq") || (print STDERR "Cannot open file ${input}_1.fastq\n"  & exit 1);
			open($fh2, "<", "${input}_2.fastq") || (print STDERR "Cannot open file ${input}_2.fastq\n"  & exit 1);
		}elsif (-e "${input}_1.fq.gz"){
			open($fh1, "<", "gunzip -c ${input}_1.fq.gz |") || (print STDERR "Cannot open file ${input}_1.fq.gz\n" & exit 1);
			open($fh2, "<", "gunzip -c ${input}_2.fq.gz |") || (print STDERR "Cannot open file ${input}_2.fq.gz\n" & exit 1);
		}elsif (-e "${input}_1.fq"){
			open($fh1, "<", "${input}_1.fq") || (print STDERR "Cannot open file ${input}_1.fq\n"  & exit 1);
			open($fh2, "<", "${input}_2.fq") || (print STDERR "Cannot open file ${input}_2.fq\n"  & exit 1);
		}else{
			print STDERR "Cannot find fastq file with prefix ${input}\n";
			exit 1;
		}
	
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
			chomp($quality1);
			chomp($quality2);
			chomp($sequence2);
			
			my @alleles = ();
			my $allele_cnt = 0;
			foreach my $pair (@candidate_pairs){
				my @pair_parts = split /\|/, $pair;
				my $candidate1 = $pair_parts[0];
				my $candidate2 = $pair_parts[1];
				my $in_pair = 0;
				
				##These reads map perfectly to one of the alleles in the pair 
				if($candidate_reference{$candidate1}=~m/$sequence1/){
					my $sequence2r = SequenceFunctions->revcomp($sequence2);
					if($candidate_reference{$candidate1}=~m/$sequence2r/){
						alignRead($sequence1,$sequence2,1,$quality1, $quality2, $opts{kmer}, $opts{mismatch}, $opts{min_qual},\%reference_hash,\%candidate_reference,\%mismatches,\%reference_counts,\@candidates);
						$in_pair=1;
					}
				}
				if($candidate_reference{$candidate1}=~m/$sequence2/){
					my $sequence1r = SequenceFunctions->revcomp($sequence1);
					if($candidate_reference{$candidate1}=~m/$sequence1r/){
						alignRead($sequence1,$sequence2,1,$quality1, $quality2, $opts{kmer}, $opts{mismatch}, $opts{min_qual},\%reference_hash,\%candidate_reference,\%mismatches,\%reference_counts,\@candidates);
						$in_pair=1;
					}
				}
				if($candidate_reference{$candidate2}=~m/$sequence1/){
					my $sequence2r = SequenceFunctions->revcomp($sequence2);
					if($candidate_reference{$candidate2}=~m/$sequence2r/){
						alignRead($sequence1,$sequence2,1,$quality1, $quality2, $opts{kmer}, $opts{mismatch}, $opts{min_qual},\%reference_hash,\%candidate_reference,\%mismatches, \%reference_counts,\@candidates);
						$in_pair=1;
					}
				}
				if($candidate_reference{$candidate2}=~m/$sequence2/){
					my $sequence1r = SequenceFunctions->revcomp($sequence1);
					if($candidate_reference{$candidate2}=~m/$sequence1r/){
						alignRead($sequence1,$sequence2,1,$quality1, $quality2, $opts{kmer}, $opts{mismatch}, $opts{min_qual},\%reference_hash,\%candidate_reference,\%mismatches, \%reference_counts,\@candidates);
						$in_pair=1;
					}
				}
				
				if($in_pair ==0){
					##The read does no map perfectly
					$extra_kmers{$sequence1}=$sequence2;
					alignRead($sequence1,$sequence2,1,$quality1, $quality2, $opts{kmer}, $opts{mismatch}, $opts{min_qual},\%reference_hash,\%candidate_reference,\%mismatches,\%reference_counts,\@candidates);
				}
			}
			foreach my $allele (@alleles){
				$candidate_counts{$allele}++ if ($allele_cnt <= $opts{allele_threshold} && $allele_cnt > 0);
			}
		}
	}
	if($opts{print_kmers}){	
		open(KOUT, ">$opts{output}.extra_kmers.txt");
		foreach my $kmer(sort {$extra_kmers{$b} cmp $extra_kmers{$a}} keys %extra_kmers){
			print KOUT "$kmer\t$extra_kmers{$kmer}\n";
		}
		close(KOUT);
	}

	for my $allele (keys %mismatches){
		my $length = length($candidate_reference{$allele});
		#This attempts to assemble and report sequences that may fall at the beginning or end of the current allele
		if(exists $mismatches{$allele}{$length}){
			$mismatches{$allele}{$length}=compileComplexMutations($mismatches{$allele}{$length});
		}
		if(exists $mismatches{$allele}{1}){
			$mismatches{$allele}{1}=compileComplexMutations($mismatches{$allele}{1});	
		}
	}

	if($opts{print_coverage}){	
		open(OUT, ">$opts{output}.coverage.txt");
	}
	my %coverage_stats = ();
	for my $allele (sort {$a cmp $b} keys %reference_counts){
		my $max_continuous = 0;
		my %pos_hash = %{$reference_counts{$allele}};
		my $continuous = 0;
		my $covered = 0;
		for (my $pos=1; $pos<=$allele_sizes{$allele}; $pos++){
			my $ref_count = $reference_counts{$allele}{$pos} || 0;
			if($ref_count == 0){
				$max_continuous = $continuous if ($continuous > $max_continuous);
				$continuous = 0;
			}else{
				$continuous++;
				$covered++;
			}
			print OUT "$allele\t$pos\t$ref_count\n" if($opts{print_coverage});
		}
		$max_continuous = $continuous if ($continuous > $max_continuous);
		my $perc = 0;
		$perc = sprintf("%.2f", $max_continuous/$allele_sizes{$allele}*100) if($allele_sizes{$allele} != 0);
		my $cov_perc = 0;
		$cov_perc = sprintf("%.2f", $covered/$allele_sizes{$allele}*100) if($allele_sizes{$allele} != 0);
		$coverage_stats{$allele}{mcon}=$max_continuous;
		$coverage_stats{$allele}{pcon}=$perc;
		$coverage_stats{$allele}{mcov}=$covered;
		$coverage_stats{$allele}{pcov}=$cov_perc;
	}
	if($opts{print_coverage}){ 
		close(OUT);
	}
	
	if($opts{print_mismatches}){
		open(MOUT, ">$opts{output}.mismatches.txt");
		print MOUT "Allele\tPositions\tMismatchBase\tMismatchCoverage\tReferenceCoverage\n";
	}

	my %allele_hash=();
	foreach my $allele (@candidates){
		$allele_hash{$allele}{replace}=0;
	}
	for my $allele (sort {$a cmp $b} keys %mismatches){
		my %real_mutations = (); 
		my %positions = %{$mismatches{$allele}};
		foreach my $pos (sort {$a<=>$b} keys %positions){
			my %bases = %{$positions{$pos}};
			foreach my $base (sort {$a cmp $b} keys %bases){
				my $ref_counts = $reference_counts{$allele}{$pos} || 0;
				if ($base ne "N" && $bases{$base}>$opts{counts_threshold}){
					print MOUT "$allele\t$pos\t$base\t$bases{$base}\t$ref_counts\n"  if($opts{print_mismatches});
					my $freq = $bases{$base}/($ref_counts+$bases{$base});
					if(length($base)==1 && $freq>=$opts{minimum_frequency}){
						if(defined $real_mutations{$pos}{freq}){
							if($freq > $real_mutations{$pos}{freq}){
								$real_mutations{$pos}{base} = $base;
								$real_mutations{$pos}{freq} = $freq;
							}	
						}else{
							$real_mutations{$pos}{base} = $base;
							$real_mutations{$pos}{freq} = $freq;
						}
					}
				}
			}	
		}
		if(scalar keys %real_mutations > 0){
			##Uses the mismatches to guess a novel allele sequence and detremine if the sequence is novel or an updated version of an existing partial allele
			my ($new_allele,$new_name,$new_seq, $novel) = identifyNovelAllele($allele,\%real_mutations, \%candidate_reference,\%gene_references);
			$allele_hash{$allele}{allele}=$new_allele;
			$allele_hash{$allele}{name}=$new_name;
			$allele_hash{$allele}{seq}=$new_seq;
			$allele_hash{$allele}{replace}=1;
			$allele_hash{$allele}{novel}=$novel;
		}
	}

	if($opts{print_mismatches}){
		close MOUT;
	}
	
	return(\%allele_hash,\%coverage_stats);
}

sub alignRead{
	my $read1 = shift;
	$read1 = shift if ($read1 eq "AlleleRefiner");
	my $read2 = shift;
	my $count = shift;
	my $read1r = SequenceFunctions->revcomp($read1);
	my $read2r = SequenceFunctions->revcomp($read2);
	my $length1 = length($read1r);
	my $length2 = length($read2r);
	my $quality1 = shift;
	my $quality2 = shift;
	my $quality1r = SequenceFunctions->rev($quality1);
	my $quality2r = SequenceFunctions->rev($quality2);
	my %positions;
	my $kmer_size = shift;
	my $mthresh = shift;
	my $qthresh = shift;
	my $reference_hash = shift;
	my $candidate_ref = shift;
	my $mismatches = shift;
	my $reference_counts = shift;
	my $candidate_array_ref = shift;
	my @candidates = @$candidate_array_ref;
	
	#Finds possible mappings for each k-mer in read 1
	for(my $i=0; $i<=($length1-$kmer_size); $i++){
		my $k1 = substr($read1,$i,$kmer_size);
		my $k1r = substr($read1r,$i,$kmer_size);

		if(exists $reference_hash->{$k1}){
			foreach my $chrpos (keys %{$reference_hash->{$k1}}){
				my @parts = split /\:\:/, $chrpos;	
				my $pos = $parts[1]-$i;
				$positions{k1}{$parts[0]}{$pos}++;
			}	
		}
		if(exists $reference_hash->{$k1r}){
			foreach my $chrpos (keys %{$reference_hash->{$k1r}}){
				my @parts = split /\:\:/, $chrpos;	
				my $pos = $parts[1]-$i;
				$positions{k1r}{$parts[0]}{$pos}++;
			}	
		}
	}
	#Finds possible mappings for each k-mer in read 2
	for(my $i=0; $i<=($length2-$kmer_size); $i++){
		my $k2 = substr($read2,$i,$kmer_size);
		my $k2r = substr($read2r,$i,$kmer_size);
		if(exists $reference_hash->{$k2}){
			foreach my $chrpos (keys %{$reference_hash->{$k2}}){
				my @parts = split /\:\:/, $chrpos;	
				my $pos = $parts[1]-$i;
				$positions{k2}{$parts[0]}{$pos}++;
			}	
		}
		if(exists $reference_hash->{$k2r}){
			foreach my $chrpos (keys %{$reference_hash->{$k2r}}){
				my @parts = split /\:\:/, $chrpos;	
				my $pos = $parts[1]-$i;
				$positions{k2r}{$parts[0]}{$pos}++;
			}	
		}
	}
	
	my %alignments = ();
	$alignments{k1}={};
	$alignments{k2}={};
	$alignments{k1r}={};
	$alignments{k2r}={};
	my %mmc = ();
	$mmc{k1} = $mthresh;
	$mmc{k2} = $mthresh;
	$mmc{k1r} = $mthresh;
	$mmc{k2r} = $mthresh;
	my %k1 = (); 
	%k1 = %{$positions{k1}} if (defined $positions{k1}); 
	my %k2 = (); 
	%k2 = %{$positions{k2}} if (defined $positions{k2}); 
	my %k1r = (); 
	%k1r = %{$positions{k1r}} if (defined $positions{k1r}); 
	my %k2r = (); 
	%k2r = %{$positions{k2r}} if (defined $positions{k2r}); 

	##Extends the reads to tally the number of mismatches with the reference
	foreach my $allele (@candidates){
		if(defined($k1{$allele}) && defined($k2r{$allele})){
			my %k1_pos=%{$k1{$allele}};
			my %k2r_pos=%{$k2r{$allele}};
			if(scalar keys %k1_pos >0 && scalar keys %k2r_pos >0){
				$mmc{k1} = identifyMismatches($allele,$length1,$read1,\%k1_pos,$alignments{k1}, $mmc{k1},$quality1,$candidate_ref, $mthresh, $qthresh);		
				$mmc{k2r} = identifyMismatches($allele,$length2,$read2r,\%k2r_pos,$alignments{k2r}, $mmc{k2r}, $quality2r,$candidate_ref,  $mthresh, $qthresh);		
			}
		}
		if(defined($k2{$allele}) && defined($k1r{$allele})){
			my %k2_pos=%{$k2{$allele}};
			my %k1r_pos=%{$k1r{$allele}};
			if(scalar keys %k2_pos >0 && scalar keys %k1r_pos >0){
				$mmc{k2} = identifyMismatches($allele,$length2,$read2,\%k2_pos,$alignments{k2}, $mmc{k2}, $quality1r,$candidate_ref,  $mthresh, $qthresh);		
				$mmc{k1r} = identifyMismatches($allele,$length1,$read1r,\%k1r_pos,$alignments{k1r}, $mmc{k1r}, $quality2,$candidate_ref,  $mthresh, $qthresh);		
			}
		}
			
	}

	my @orientations = ("k1_k2r","k2_k1r");
	my $r1_counts = 0;
	my $r2_counts = 0;
	my %mapped_alleles = ();
	
	##This examines both reads to see which alleles are mapped. Only pairs with both reads mapping will be considered as mapped
	foreach my $allele (@candidates){
		if(exists ($alignments{'k1'}{$mmc{'k1'}}) && exists ($alignments{'k2r'}{$mmc{'k2r'}})){
			my $r1cnt = 0;
			my $r2cnt = 0;
			$r1cnt = $alignments{'k1'}{$mmc{'k1'}}{counts}{$allele} > 0 if (exists($alignments{'k1'}{$mmc{'k1'}}{counts}{$allele}));
			$r2cnt = $alignments{'k2r'}{$mmc{'k2r'}}{counts}{$allele} > 0 if (exists($alignments{'k2r'}{$mmc{'k2r'}}{counts}{$allele}));
			if($r1cnt > 0 && $r2cnt > 0){
				$r1_counts += $r1cnt;
				$r2_counts += $r2cnt;
				$mapped_alleles{$allele}=1;	
			}
		}
		if(exists ($alignments{'k1r'}{$mmc{'k1r'}}) && exists ($alignments{'k2'}{$mmc{'k2'}})){
			my $r1cnt = 0;
			my $r2cnt = 0;
			$r1cnt = $alignments{'k1r'}{$mmc{'k1r'}}{counts}{$allele} > 0 if (exists($alignments{'k1r'}{$mmc{'k1r'}}{counts}{$allele}));
			$r2cnt = $alignments{'k2'}{$mmc{'k2'}}{counts}{$allele} if (exists($alignments{'k2'}{$mmc{'k2'}}{counts}{$allele}));
			if($r1cnt > 0 && $r2cnt > 0){
				$r1_counts += $r1cnt;
				$r2_counts += $r2cnt;
				$mapped_alleles{$allele}=1;	
			}
		}
	}
	my %counts = ();
	$counts{'k1'} = $r1_counts;	
	$counts{'k1r'} = $r1_counts;	
	$counts{'k2'} = $r2_counts;	
	$counts{'k2r'} = $r2_counts;	
	
	##Calculates reference and mismatch counts for positions based on the best match. Counts are weighted based on the numer best matches and the quality of the base.
	foreach my $opairs(@orientations){
		my @pairs = split "_", $opairs;
		if(exists ($alignments{$pairs[0]}{$mmc{$pairs[0]}}) && exists ($alignments{$pairs[1]}{$mmc{$pairs[1]}})){
			foreach my $orient(@pairs){
				my $alignment_count = $counts{$orient};
				my $ratio = 0;
				$ratio = 1/$alignment_count if ($alignment_count != 0);
				
				#Count should always be 1 in the current framework
				my $score = $ratio * $count;
				
				foreach my $allele (keys %mapped_alleles){
					my %best_alignments = %{$alignments{$orient}{$mmc{$orient}}{alignments}{$allele}};
					foreach my $start_pos (keys %best_alignments){
						foreach my $mut_pos (sort {$a<=>$b} keys %{$best_alignments{$start_pos}}){
							my $base = $best_alignments{$start_pos}{$mut_pos};
							my $qual_score = $alignments{$orient}{$mmc{$orient}}{scores}{$allele}{$start_pos}{$mut_pos};
							if( !defined($alignments{$orient}{$mmc{$orient}}{scores}{$allele}{$start_pos}{$mut_pos})|| $alignments{$orient}{$mmc{$orient}}{scores}{$allele}{$start_pos}{$mut_pos} eq ""){
								$qual_score = 1;
							}
							$mismatches->{$allele}{$mut_pos}{$base}+=($score*$qual_score);		
						}
					}
				}
				foreach my $allele (keys %mapped_alleles){
					my %refpos = %{$alignments{$orient}{$mmc{$orient}}{ref_counts}{$allele}};
					foreach my $start_pos (keys %refpos){
						foreach my $ref_pos (sort {$a<=>$b} keys %{$refpos{$start_pos}}){
							my $qual_score = $alignments{$orient}{$mmc{$orient}}{scores}{$allele}{$start_pos}{$ref_pos};
							my $refcount = $refpos{$start_pos}{$ref_pos};
							$reference_counts->{$allele}{$ref_pos}+=($score*$refcount*$qual_score);		
						}
					}
				}
			}
		}
	}
}

sub identifyMismatches{
	my $allele = shift;
	$allele = shift if ($allele eq "AlleleRefiner");
	my $length = shift;
	my $read = shift;
	my $ref =shift;
	my $alignment_ref = shift;
	my $minimum_mismatch_count = shift;
	my $quality = shift;
	my $candidate_ref = shift;
	my $mthresh = shift;
	my $qthresh = shift;
	my %positions = %{$ref};

	foreach my $pos (keys %positions){
		my @seq = split //, $read;
		my @qual = split //, $quality;
		
		my @reference = split //, substr($candidate_ref->{$allele},$pos-1,$length);
		if($pos <=0){
			@reference = split //, substr($candidate_ref->{$allele},0,$length+$pos-1);
		}
		my $mismatch_count=0;
		my %mismatch_hash = ();
		my %possible_alignment;
		my %ref_hash =();
		##Truncates the sequences with overhang on front or end of reference sequence
		if(($#reference+1)<$length){
			if(($pos + $length) > length($candidate_ref->{$allele})){
				my @tmp_seq = ();
				my @tmp_qual = ();
				my $min_qual = 100;
				for (my $i=0;$i<=$#seq; $i++){
					if($i < $#reference){
						$tmp_seq[$i]=$seq[$i];
						$tmp_qual[$i]=$qual[$i];
					}else{
						#Add base of overhang to last position
						$tmp_seq[$#reference].=$seq[$i];
						my $q = ord($qual[$i]) - 33;
						#Stores the lowest quality of the entire overhang as the quality of last position
						if($q<$min_qual){
							$min_qual = $q;
							$tmp_qual[$#reference]=$qual[$i];
						}
					}
				}
				@seq=@tmp_seq;
				@qual=@tmp_qual;
			}elsif($pos <= 0){
				my @tmp_seq = ();
				my @tmp_qual = ();
				my $index = 1;
				my $min_qual = 100;
				for (my $i=0;$i<=$#seq; $i++){
					if($i <= ($length-$#reference-1)){
						#Add base of overhang to position 0
						$tmp_seq[0] .= $seq[$i];
						$tmp_qual[0] .= $qual[$i];
						my $q = ord($qual[$i]) - 33;
						#Stores the lowest quality of the entire overhang as the quality of 0 position
						if($q<$min_qual){
							$min_qual = $q;
							$tmp_qual[0]=$qual[$i];
						}
					}else{
						$tmp_seq[$index] = $seq[$i];
						$tmp_qual[$index]= $qual[$i];
						$index++;
					}
				}
				@seq=@tmp_seq;	
				@qual=@tmp_qual;	
				$pos = 1;
			}
		}
		my %score_hash = ();

		##Counts mismatches and references
		for(my $i=0; $i<=$#seq; $i++){
			my $mpos = $pos + $i;
			my $mqual = ord($qual[$i])-33;
			my $score = 0;
			#This section determines the weight of the Base based on quality. Base with quality below the specified quality are weighted 0 and 35 and above are weighted 1. 
			if($mqual >= 35){
				$score = 1;
			}elsif ($mqual >= $qthresh){
				my $diff = 35 - $qthresh;
				$score = ($mqual - $qthresh)/$diff;
			}
			$score_hash{$mpos} = $score;
			if($seq[$i] ne $reference[$i]){
				$mismatch_count++;
				$mismatch_hash{$mpos}=$seq[$i];
			}else{
				$ref_hash{$mpos}++;
			}
		}
		if($mismatch_count <= $mthresh){
			$minimum_mismatch_count = $mismatch_count if ($mismatch_count < $minimum_mismatch_count);
			$alignment_ref->{$mismatch_count}{alignments}{$allele}{$pos}=\%mismatch_hash;
			$alignment_ref->{$mismatch_count}{ref_counts}{$allele}{$pos}=\%ref_hash;
			$alignment_ref->{$mismatch_count}{counts}{$allele}++;
			$alignment_ref->{$mismatch_count}{scores}{$allele}{$pos}=\%score_hash;
		} 
	}
	return($minimum_mismatch_count);
}

sub compileComplexMutations{
	my $mutref = shift;
	$mutref = shift if ($mutref eq "AlleleRefiner");
	my %mut_hash = %{$mutref};
	my %new_mut_hash = ();
	my @mutations = sort {length($b)<=>length($a)} keys %mut_hash;
	my @ref=();
	push @ref, shift @mutations;
	$new_mut_hash{$ref[0]}=$mut_hash{$ref[0]};
	foreach my $mut (@mutations){
		my $in_ref = 0;
		my %in_ref_hash = ();
		foreach my $r (@ref){
			if($r =~ m/$mut/){
				$in_ref_hash{$r}=$mut_hash{$mut};
				$in_ref++;
			}
		}
		if($in_ref == 0){
			$new_mut_hash{$mut}=$mut_hash{$mut};	
			push @ref, $mut;	
		}else{
			my $score = $mut_hash{$mut} * (1/$in_ref);
			foreach my $r (keys %in_ref_hash){
				$new_mut_hash{$r}+=$score;
			}	
		}
	}
	return \%new_mut_hash;
}

sub addAlleleToReference{
	my $name = shift;
	$name = shift if ($name eq "AlleleRefiner");
	my $sequence = shift;
	my $kmer_size = shift;
	my $ref = shift;
	my $candidate_ref = shift;
	$candidate_ref->{$name}=$sequence;	
	my $length = length($sequence);
	for(my $i=0; $i<=($length-$kmer_size); $i++){
		my $pos = $i + 1;
		my $k = substr($sequence,$i,$kmer_size);
		$ref->{$k}{"$name\:\:$pos"}++;
	}
}

sub identifyNovelAllele{
	my $allele = shift;
	$allele = shift if ($allele eq "AlleleRefiner");
	my $mut_ref = shift;
	my %mutations = %{$mut_ref};
	my $candidate_ref = shift;
	my $ref_seq = $candidate_ref->{$allele};
	my $gene_ref = shift;	
	my %gene_references = %{$gene_ref};
	#Change the sequence to match mutation
	foreach my $pos (keys %mutations){
		my $ind = $pos -1;
		substr($ref_seq, $ind, 1) = $mutations{$pos}{base};
	}
	
	#Search to see if part of the novel sequence exactly matches an existing (and likely partial) allele
	my %matches;
	foreach my $allele (keys %gene_references){
		my $seq = $gene_references{$allele}{seq};
		if($ref_seq =~ m/$seq/){
			$matches{$allele}=length($seq);
		}	
	}
	my @matching_alleles = sort {$matches{$b}<=>$matches{$a}} keys %matches;
	my $allele_string;
	if($#matching_alleles >= 0){
		return("$matching_alleles[0]U","$gene_references{$matching_alleles[0]}{name}_updated", $ref_seq, 0);
	}else{
		return("${allele}N","$gene_references{$allele}{name}_novel", $ref_seq, 1)
	}
}

__PACKAGE__->runCommandline() unless caller;
