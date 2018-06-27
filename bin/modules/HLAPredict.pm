#!/usr/bin/env perl
package HLAPredict;
use strict;
use warnings;
use Getopt::Long;
use List::Util;
use Statistics::Basic;
use Module::Load;
use Storable;
use Parallel::ForkManager;

(my $SCRIPT_NAME = $0) =~ s/.*\///;
my $version = "1.0";
my $creation_date = "1 Oct 2016";
my $last_updated = "11 Jan 2017";

my $usage = "\n$SCRIPT_NAME v$version\n".
	"\nDESCRIPTION:\n" .
	"This module predicts HLATypes using a k-mer profiler. The module can be run directly from the commandline or directly from another perl script.\n" .
	"\nUSAGE:\n" .
	"perl $SCRIPT_NAME <options>\n" .
	"\nRequired Options:\n" .
	"-unique1|u\tFile containing counts for the filtered read1 reads (required)\n" .
	"-fastq_dir|fd\tFastq file containing filtered read1 reads (required)\n" .
	"-reference|r\tReference HLA fasta(required)\n" .
	"-profile|ph\tProfile perl hash stored in a Storable file.(required)\n" .
	"-allele_map|am\tFile mapping HLA allele accessions ids to HLA allele names(required)\n" .
	"\nHLA Typing Parameters:\n" .
	"-cntThresh|ct\t\tThe minimum counts required for an observed k-mer to be considered (default:2)\n" . 
	"-numSemifinalists|ns\tThe number of semifinalists considered for evaluating possible pairs of HLA Alleles (default:100)\n" .
	"-numFinalists|ns\t\tThe number of finalist allele pairs considered for HLA types (default:20)\n" .
	"-minReads|mr\t\tThe minimum number of reads required to predict profile (default:100)\n" . 
	"-PairPickerPower|pow\tThe weight of the effect of PairPicker.pm results on the final prediction score (default:0.25)\n" .
	"-minQuality|q\t\tOnly consider reads with all bases greater than this quality when running PairPicker.pm (default:20)\n" . 
	"\nAllele Refinement Parameters:\n" .
	"-allele_refinement|ar\tSpecifies the level to which the predicted alleles are to be refined based on the observed reads (default:all)\n" .
	"\t\t\tPossible values:\n" .
	"\t\t\t   refine_only\t\tRefines the allelle call by looking predicting the true allele sequence using observed reads and looking for a better match in the reference\n" .
	"\t\t\t   predict_only\t\tReports if the observe reads support a novel allele sequence not found in the reference\n" .
	"\t\t\t   refineAndPredict\tRefines the allele call (-refine_only) and report novel alleles (-novel_only)\n" . 
	"\t\t\t   all\t\t\tRefines the allele call (-refine_only) and report novel alleles (-novel_only), creates a profile for the refined/novel allele sequence and calculates prediction metrics.\n" . 
	"\t\t\t   none\t\t\tTurns off refinement and prediction.\n" . 
	"-sim_num_reads|snr\tnumber of reads to simulated per reference allele for k-mer profile creations.(default:500000)\n" .
	"-sim_read_length|srl\tlength of reads simulated for k-mer profile. Same as the length of the k-mers in the profile.(default:50)\n" .
	"-sim_max_insert|sm\tmaximum size of insert (default:1000)\n" .
	"-sim_scale|ssc\t\tscale of pareto distribution to determine insert size (default:80)\n" .
	"-sim_shape|ssh\t\tshape of pareto distribution to determine insert size (default:0.7)\n" .
	"-sim_seed|ssd\t\tseed of random number generator for simulation (default:1234)\n" .
	"-kraken_db|kdb\t\tBase directory containing the kraken database files\n" .
	"-kraken_path|kp\t\tDirectory containing kraken scripts\n" .
	"\nGeneral Options:\n" .
	"-out_dir|od\t\tlocation of output_directory (default:'.')\n" .
	"-output|o\t\toutput file (default:STDOUT)\n" .
	"-log|l\t\t\tName of log file (default:STDERR)\n" .
	"-threads|c\t\tNumber of threads (default:1)\n" .
	"-scripts_dir|sd\t\tThe base directory of the HLAProfiler scripts (default:'.')\n" .
	"-help|h\t\t\tDisplays this message\n" .
        "\nAUTHORS:\n" .
	"Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	"Chad Brown\n" .
	"\nCREATED:\n$creation_date\n" .
	"\nLAST UPDATED:\n$last_updated\n" .
	"\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
	"\n";	

my %opts = (sim_num_reads=>500000,sim_read_length=>50, sim_max_insert=>1000, sim_seed=>1234,sim_scale=>80, sim_shape=>0.7,pp_scale=>0,threads=>1, unique1=>"",fastq_dir=>"", log=>"", cntThresh=>2, scripts_dir=>".",minReads=>100, reference=>"", minQuality=>20, profile=>"", output=>"", numFinalists=>20,numSemifinalists=>100, PairPickerPower=>0.25, allele_refinement=>"all", out_dir => ".", kraken_db=>".", kraken_path=>".", allele_map=>"");
my %sig = ();
my %kmer_keys = ();
my %allele_cnts = ();
my %sig_kmers = ();
my %combined_kmers = ();
my $log;
sub run_Commandline{
	GetOptions(\%opts, qw(help|h sim_num_reads|snr=s sim_read_length|srl=s sim_max_insert|sm=s sim_scale|ssc=s sim_shape|ssh=s sim_seed|ssd=s allele_refinement|ar=s pp_scale|ps=s threads|c=s fastq_dir|fd=s unique1|u=s cntThresh|ct=s scripts_dir|sd=s reference|r=s minReads|mr=s minQuality|q=s output|o=s profile|ph=s numFinalists|nf=s numSemifinalists|ns=s PairPickerPower|pow=s allele_map|am=s log|l=s out_dir|od=s kraken_path|kp=s kraken_db|kdb=s));
	if ($opts{help}){
		print "$usage";
		exit
	}elsif($opts{reference} eq ""){
		print "Please specify a reference.\n$usage\n";
		exit 1;
	}elsif ($opts{unique1} eq ""){
		print "Please specify a unique counts file.\n$usage\n";
		exit 1;
	}elsif ($opts{fastq_dir} eq ""){
		print "Please specify a fastq directory.\n$usage\n";
		exit 1;
	}elsif ($opts{profile} eq ""){
		print "Please specify a profile file.\n$usage\n";
		exit 1;
	}elsif ($opts{allele_map} eq ""){
		print "Please specify an allele_map file.\n$usage\n";
		exit 1;
	}elsif(! -s $opts{allele_map}){
		print "$opts{allele_map} not valid file. Please specify a valid -allele_map\n";
		exit 1;
	}
	my $answer = predictHLAType($opts{unique1},$opts{fastq_dir},$opts{profile},$opts{reference},$opts{allele_map},$opts{log});
	my $out= *STDOUT;
	if($opts{output}){
		open($out, ">$opts{output}");
	}
	print $out "Allele1_Accession\tAllele2_Accession\tAllele1\tAllele2\tProportion_reads\tProportion_signal\tCorrelation\tError\tPair_score\tFinal_score\tAllele1 Comments\tAllele2 Comments\n";
	print $out "$answer";
	close($out);
}

sub setOptions{
	$opts{cntThresh}=shift;
	$opts{cntThresh}=shift if ($opts{cntThresh} eq "HLAPredict");
	$opts{numFinalists}=shift;
	$opts{numSemifinalists}=shift;
	$opts{PairPickerPower}=shift;
	$opts{minQuality} = shift;
	$opts{threads} = shift;
	$opts{pp_scale}=shift;
	$opts{out_dir}=shift;
	$opts{scripts_dir}=shift;
	$opts{allele_refinement}=shift;
	my $opts_ref = shift;
	my %simOpts = %{$opts_ref};
	foreach my $o (keys %simOpts){
		$opts{"sim_$o"}=$simOpts{$o};
	}
	$opts{kraken_db}=shift;
	$opts{kraken_path}=shift;
	$opts{minReads}=shift;
	
	return \%opts;
}

sub predictHLAType{
	$opts{unique1} = shift;
	$opts{unique1} = shift if ($opts{unique1} eq "HLAPredict");
	$opts{fastq_dir} = shift;
	$opts{profile} = shift;
	$opts{reference} = shift;
	$opts{allele_map} = shift;
	$opts{log} = shift;
	$opts{output_dir} = shift;       
 
        if($opts{log}){
                open($log, ">>$opts{log}");
        }else{
                $log = *STDERR;
	}
	
	load "$opts{scripts_dir}/modules/SequenceFunctions.pm";
	$opts{gene} = $opts{unique1};
	$opts{gene} =~ s/\_1\.uniq\.cnts$//;
	$opts{gene} =~ s/^.+\.//;
	
	my $prefix = $opts{unique1};
	$prefix =~s/\_1\.uniq\.cnts$//;
	
	## The allele map file is tab-delimeted. Accesssion<TAB>Allele<TAB>Comments containing merged sequences<TAB>1 or 0 specifying if the allele is common
	my %allele_map;
	open (IN, $opts{allele_map});
	while(<IN>){
		chomp;
		$_=~s/HLA://;
		$_=~s/IPD://;
		my @cols=split /\t/, $_;
		$allele_map{$cols[0]}{names}=$cols[1];
		$allele_map{$cols[0]}{comments}=$cols[2];
		$allele_map{$cols[0]}{cwd}=$cols[3];
	}
	close(IN);
	
	## Uses Storable::retrieve to load in the signaure file, which contains a perl hash
	my $sig_ref = "";
	$sig_ref = retrieve($opts{profile});
	%sig =();
	%kmer_keys = ();
	%allele_cnts = ();
	%sig_kmers = ();

	## %sig is a hash of allele accesion ids with each value being another hash with a k-mer as the key and corresponding read count as the value
	## To save on storage space, the k-mers in %sig have been replaced with a unique identifier 
	%sig = %{$sig_ref->{sig}};
	
	## Stores the total number of counts for each allele
	%allele_cnts = %{$sig_ref->{totals}};

	
	## key: unique identified value:kmer
	%kmer_keys = %{$sig_ref->{names}};

	## key: kmer value:unique identifier
	%sig_kmers = %{$sig_ref->{kmers}};
	
	my @allele_names = keys %sig;
	my $kmer1 = (keys %{$sig{$allele_names[0]}})[0];
	my $klen = length($kmer_keys{$kmer1});

	## Read in count files and determine the observed k-mer profile
	my ($read_cnt, $total_cnts, $mx) = readUnique12($prefix, $klen, \%sig_kmers, \%sig);

	if($read_cnt < $opts{minReads}){
		print $log "$opts{gene} has less than $opts{minReads} reads.\n";
		return "HLA\tHLA\t$opts{gene}\t$opts{gene}\t-\t-\t-\t-\t-\t-\tNot enough reads to make call\t\n";
	}
		
	if(! defined $sig{cnt}){
		print $log "$opts{gene} has no onbserved counts for the profile k-mers.\n";
		return "HLA\tHLA\t$opts{gene}\t$opts{gene}\t-\t-\t-\t-\t-\t-\tNot enough reads to make call\t\n";
	}
	## Identify the top 100 alleles whose profile is most accounted for by the reads
	my %frac = ();

	##As sequencing depth increases the effect of sequencing error on allele prioritization increases we will scale the kmer counts of each signature accordingly to prioritize variants
	my $lower_thresh = int($read_cnt/50000*$mx);
	my $upper_thresh = int($read_cnt/5000*$mx);
	$lower_thresh = 1 if ($lower_thresh < 1);
	$lower_thresh = 1; ##REMOVE AFTER DEBUGGING		
	my $range = $upper_thresh - $lower_thresh;
	
	## Iterate over alleles
	ALLELE:while (my ($allele, $asig) = each(%sig) ) {
		my $total = 0;
	
			if($allele eq "cnt"){ ## This is the observed profile and so cnt should not be considered as a candidate allele
				$total = -1;
			}else{
				##Iterate over k-mers
				while (my ($kmer, $value) = each(%{$sig{$allele}}) ) {
					if( defined $sig{cnt}{$kmer_keys{$kmer}} ) {
						my $scale = 1; 
						if($value >= $upper_thresh){
							#$scale = 1
						}elsif($value < $upper_thresh && $value > $lower_thresh){
							#$scale = ($value-$lower_thresh)/($range)
						}elsif($value <= $lower_thresh){
							#$scale = 0;
						}
						$total+=($value*$scale);
					}else{
						my $tmp = 0;
					}
				}
			}
			$frac{$allele}=$total;
	}
	
	## %frac contains the counts of kmers for the allele observed in the data 
	foreach my $allele (keys %frac){
		$frac{$allele}=$frac{$allele}/$allele_cnts{$allele} if ($allele ne "cnt"); #Normalize by count of observed kmers by total counts in the profile 
	}
	my $numSemifinalists = $opts{numSemifinalists};
	my @topAlleles = sort {$frac{$b} <=> $frac{$a}} keys %frac;

	if($#topAlleles<$numSemifinalists){
		$numSemifinalists=$#topAlleles;
	}else{
		my $index = $numSemifinalists;
		my $min_frac = $frac{$topAlleles[$numSemifinalists-1]};
		##This extends the number of semifinalists to include those below the hard cut-off but having the same fraction as the last semi finalist
		while($frac{$topAlleles[$index]} == $min_frac && $index<scalar keys %frac){
			$numSemifinalists++;
			$index++;
		}
	}
	
	@topAlleles = @topAlleles[0..($numSemifinalists-1)];

	#Create union set of primers to be used in correlations
	%combined_kmers=();
	findCombinedKmers(\%sig,\%kmer_keys,\@topAlleles, $lower_thresh);
	my %allele_pairs = ();
	my $offset = 1 /($read_cnt*2);
	
	##For each pair of possible alleles. Calculate:
	##Correlation of observed profile and profile profile
	##Proportion of signal accounted for by the observed reads
	##Proportion of reads that account for the signal
	##Error (1-proportion of reads)+(1-proportion of signal)+(1-correlation)
	for(my $i =0; $i<=$#topAlleles; $i++){
		for(my $j =0; $j<=$#topAlleles; $j++){
			if($i<=$j){
				my ($cor, $propReads, $propSignal) = getSSE($topAlleles[$i], $topAlleles[$j], $total_cnts, ($allele_cnts{$topAlleles[$i]} + $allele_cnts{$topAlleles[$j]}), $offset, \%sig);
				#my @result = ("$i$j", $topAlleles[$i], $topAlleles[$j],$cor, $propReads, $propSignal);
				#$fm->finish(0,\@result);
		
				$allele_pairs{"$i$j"}{allele1}=$topAlleles[$i];	
				$allele_pairs{"$i$j"}{allele2}=$topAlleles[$j];
				$allele_pairs{"$i$j"}{comment1}=$allele_map{$topAlleles[$j]}{comment};
				$allele_pairs{"$i$j"}{comment2}=$allele_map{$topAlleles[$j]}{comment};
				$allele_pairs{"$i$j"}{cor}=$cor;	
				$allele_pairs{"$i$j"}{propReads}=$propReads;	
				$allele_pairs{"$i$j"}{propSignal}=$propSignal;	
				$allele_pairs{"$i$j"}{error}=(1-$propSignal)+(1-$propReads)+(1-$cor);	
			}
		}
	}
	
	##Identifies the finalists and creates the candidate list of the PairPicker module
	my $pair_cnt=1;
	my $candidates = "";
	my @finalists=();
	my $prev_error = "";
	foreach my $pair (sort {$allele_pairs{$a}{error}<=>$allele_pairs{$b}{error}} keys %allele_pairs){
		if($pair_cnt>$opts{numFinalists} && $allele_pairs{$pair}{error} != $prev_error){
			last;
		}
		elsif($pair ne ""){
			$candidates .= "$allele_pairs{$pair}{allele1}|$allele_pairs{$pair}{allele2},";
			#push(@finalists, "$allele_pairs{$pair}{allele1}|$allele_pairs{$pair}{allele2}";
			$prev_error = $allele_pairs{$pair}{error};
			push(@finalists, "$pair");
			$pair_cnt++;
		}
	}
	$candidates=~s/,$//;
	
	
	my $fq_prefix = $prefix;
	$fq_prefix =~s/.*\///;
	$fq_prefix=~s/\.$opts{gene}$//;

	##The PairPicker module does pairwise comparison between candidate allele pairs to identify reads the map perfectly and uniquely to one allele pair but not the other. This returns a matrix. For each comparison of row(r) vs. col(c) the matrix contains the count of reads unique to r in [r,c] and unique to c in [c,r]
	load "$opts{scripts_dir}/modules/PairPicker.pm";
	my ($pp_results_ref,$pp_candidates_ref, $pp_unique_counts) = PairPicker->runModule("$opts{fastq_dir}/$fq_prefix", $opts{reference},$candidates,$opts{gene}, $opts{minQuality}, $opts{pp_scale}, $opts{scripts_dir}, 0);
	
	##Sums the columns and rows from the PairPicker matrix. The row counts are counts unique to the allele pair of the row for that allele pair and the column counts contain the counts   
	my $sums_ref = sumMatrix($pp_results_ref,\@finalists,($read_cnt*2)/10000);
	my %matrix_sums = %{$sums_ref};

	my $min_pp_reads = 10;
	my $max_pp_reads = 100;
	
	my $power = $opts{PairPickerPower};
	##PairPicker matrices with low counts are more sensitive to sequencing error can falsely increase the score of incorrect alleles 
	if($pp_unique_counts < $min_pp_reads){
		$power = 0;	
	}elsif($pp_unique_counts < $max_pp_reads){
		$power = ($pp_unique_counts/$max_pp_reads)*$power;
	}	

	##Calculates the pp_score and then the final score
	my %final_scores = ();
	foreach my $pair (keys %matrix_sums){
		$allele_pairs{$pair}{error} = .000001 if ($allele_pairs{$pair}{error}==0);
		$allele_pairs{$pair}{pp_score} = ($matrix_sums{$pair}{row}/$matrix_sums{$pair}{col})**$power;
		if($power == 0){
			$allele_pairs{$pair}{pp_score} = 1;
		}
		$matrix_sums{$pair}{col}=1 if ($matrix_sums{$pair}{col}==0);
		$final_scores{$pair} = $allele_pairs{$pair}{pp_score}/$allele_pairs{$pair}{error};
	}
	
	my $max_score = 0;
	
	my %answers = ();
	my $answer = "";
	my $answer_cnt = 0;
	my $refine = 0;
	my $predict = 0;
	my $recalculate = 0;

	if($opts{allele_refinement} eq "all"){
		$refine = 1;
		$predict = 1;
		$recalculate = 1;
	}elsif($opts{allele_refinement} eq "refine_only"){
		$refine = 1;
	}elsif($opts{allele_refinement} eq "predict_only"){
		$predict = 1;
	}elsif($opts{allele_refinement} eq "refineAndPredict"){
		$refine = 1; 
		$predict = 1;
	}elsif($opts{allele_refinement} eq "none"){
	}else{
		print $log "-allele_refinement entry $opts{allele_refinement} is not a valid option. Using default of \"all\"\n"; 
		$refine = 1;
		$predict = 1;
		$recalculate = 1;
	}
 
	my @topPairs = sort {$final_scores{$b}<=>$final_scores{$a}} keys %final_scores;
	
	## This will run AlleleRefiner which predicts updated or novel allele sequences based on the observed data
	if($refine == 1 || $predict ==1){
		my $pair = $topPairs[0];
		load "$opts{scripts_dir}/modules/AlleleRefiner.pm";
		my ($allele_ref, $cov_ref) = AlleleRefiner->runModule($opts{gene},"$opts{fastq_dir}/$fq_prefix.$opts{gene}","$allele_pairs{$pair}{allele1}|$allele_pairs{$pair}{allele2}",$opts{reference},20,2,20,1,.75,$opts{scripts_dir},"$opts{out_dir}/$opts{gene}.refinement" ,1,0,0);
		my $original_allele1 = $allele_pairs{$pair}{allele1};
		my $original_allele2 = $allele_pairs{$pair}{allele2};
		my $new_name1 = $allele_map{$original_allele1}{names};
		my $new_name2 = $allele_map{$original_allele2}{names};
		my $new_allele1 = $original_allele1;
		my $new_allele2 = $original_allele2;
		my $rcomment1 = "";
		my $rcomment2 = "";

		if($allele_ref->{$original_allele1}{replace}==1){
			$new_allele1 = $allele_ref->{$original_allele1}{allele};
			if($allele_ref->{$original_allele1}{novel} ==1 && $predict==1){
				$new_name1 = $allele_ref->{$original_allele1}{name};	
				$rcomment1 =" Predicted novel allele based on differences between top allele $original_allele1 and observed data.";
			}elsif($allele_ref->{$original_allele1}{novel} ==0){
				$new_name1 = $allele_ref->{$original_allele1}{name};	
				$rcomment1 =" $new_allele1 is a perfect subset of predicted allele1 sequence based on $original_allele1 and observed reads."; 	
			}	
		}
		if($allele_ref->{$original_allele2}{replace}==1){
			$new_allele2 = $allele_ref->{$original_allele2}{allele};
			if($allele_ref->{$original_allele2}{novel} ==1 && $predict==1){	
				$new_name2 = $allele_ref->{$original_allele2}{name};	
				$rcomment2 =" Predicted novel allele based on differences between top allele $original_allele2 and observed data.";
			}elsif($allele_ref->{$original_allele2}{novel} ==0){
				$new_name2 = $allele_ref->{$original_allele2}{name};	
				$rcomment2 =" $new_allele2 is a perfect subset of predicted allele2 sequence based on $original_allele2 and observed reads.";
			}	
		}
		if($recalculate == 1 && ($allele_ref->{$original_allele1}{replace}==1 ||$allele_ref->{$original_allele2}{replace}==1)){
			##Add new sequence to profile
			my ($cor, $propReads, $propSignal, $new_reference) = createProfileAndCalculateSignal($opts{gene},$new_allele1, $new_allele2, $new_name1, $new_name2, $allele_ref->{$original_allele1}{replace}, $allele_ref->{$original_allele2}{replace}, $allele_ref->{$original_allele1}{seq}, $allele_ref->{$original_allele2}{seq}, $total_cnts, $offset,\%sig, \%kmer_keys, \%sig_kmers,\%allele_cnts, $opts{reference}); 
			
			my $pair = "novel";	
			$allele_map{$new_allele1}{names}=$new_name1;
			$allele_map{$new_allele2}{names}=$new_name2;
			$allele_map{$new_allele1}{comment}=$rcomment1;
			$allele_map{$new_allele2}{comment}=$rcomment2;
			$allele_pairs{$pair}{allele1}=$new_allele1;	
			$allele_pairs{$pair}{allele2}=$new_allele2;
			$allele_pairs{$pair}{comment1}=$rcomment1;
			$allele_pairs{$pair}{comment2}=$rcomment2;
			$allele_pairs{$pair}{cor}=$cor;	
			$allele_pairs{$pair}{propReads}=$propReads;	
			$allele_pairs{$pair}{propSignal}=$propSignal;	
			$allele_pairs{$pair}{error}=(1-$propSignal)+(1-$propReads)+(1-$cor);
			$candidates .= ",$new_allele1|$new_allele2";
			push @finalists, "novel";	
			
			##Re-run pair picker with new allele pair
			($pp_results_ref,$pp_candidates_ref, $pp_unique_counts) = PairPicker->runModule("$opts{fastq_dir}/$fq_prefix", $new_reference,$candidates, $opts{gene}, $opts{minQuality}, $opts{pp_scale}, $opts{scripts_dir}, 0);
			$sums_ref = sumMatrix($pp_results_ref,\@finalists,($read_cnt*2)/10000);
			%matrix_sums = %{$sums_ref};
			
			$power = $opts{PairPickerPower};
			if($pp_unique_counts < 10){
				$power = 0;	
			}elsif($pp_unique_counts < 100){
				$power = (($pp_unique_counts-10)/90)*$power;
			}	

			%final_scores = ();
			foreach my $pair (keys %matrix_sums){
				$allele_pairs{$pair}{error} = .000001 if ($allele_pairs{$pair}{error}==0);
				$allele_pairs{$pair}{pp_score} = ($matrix_sums{$pair}{row}/$matrix_sums{$pair}{col})**$power;
				if($power == 0){
					$allele_pairs{$pair}{pp_score} = 1;
				}
				$matrix_sums{$pair}{col}=1 if ($matrix_sums{$pair}{col}==0);
				$final_scores{$pair} = $allele_pairs{$pair}{pp_score}/$allele_pairs{$pair}{error};		
			}		
		}elsif($allele_ref->{$original_allele1}{replace}==1 ||$allele_ref->{$original_allele2}{replace}==1){
			$allele_map{$original_allele1}{names}=$new_name1;
			$allele_map{$original_allele2}{names}=$new_name2;
			$allele_pairs{$pair}{comment1}.="$rcomment1";
			$allele_pairs{$pair}{comment2}.="$rcomment2";
			$allele_pairs{$pair}{allele1}=$new_allele1;
			$allele_pairs{$pair}{allele2}=$new_allele2;
		}
	}

	@topPairs = sort {$final_scores{$b}<=>$final_scores{$a}} keys %final_scores;

	##When the proportion of signal is low, there is an enrichment of cases where the rare allele just edges out the correct common allele. In these cases more weight is given to the common allele
	if($allele_pairs{$topPairs[0]}{propSignal}<.98 && ($final_scores{$topPairs[0]} - $final_scores{$topPairs[1]})<0.5){
		foreach my $pair (@topPairs){
			my $common = 0;
			my $cwd1 = $allele_map{$allele_pairs{$pair}{allele1}}{cwd} || 0;
			my $cwd2 = $allele_map{$allele_pairs{$pair}{allele2}}{cwd} || 0;
			$common++ if ($cwd1 ==1);
			$common++ if ($cwd2 ==1);
			my $boost = $final_scores{$pair} * (0.03*$common);
			$final_scores{$pair}+=$boost;
		}
		@topPairs = sort {$final_scores{$b}<=>$final_scores{$a}} keys %final_scores;
	}

	foreach my $pair (@topPairs){
		if($final_scores{$pair} >= $max_score){
			$max_score = $final_scores{$pair};
			$answers{$pair}{allele1}=$allele_map{$allele_pairs{$pair}{allele1}}||$allele_pairs{$pair}{allele1};
			$answers{$pair}{allele2}=$allele_map{$allele_pairs{$pair}{allele2}}||$allele_pairs{$pair}{allele2};
			$answers{$pair}{propSig}=$allele_pairs{$pair}{propSignal};
			$answers{$pair}{propReads}=$allele_pairs{$pair}{propReads};
			$answers{$pair}{cor}=$allele_pairs{$pair}{cor};
			$answers{$pair}{error}=$allele_pairs{$pair}{error};
			$answers{$pair}{pp_score}=$allele_pairs{$pair}{pp_score};
			$answers{$pair}{final_score}=$final_scores{$pair};
			my $comment1 = $allele_pairs{$pair}{comment1} || "";
			my $comment2 = $allele_pairs{$pair}{comment2} || "";
			$answer .= "$allele_pairs{$pair}{allele1}\t$allele_pairs{$pair}{allele2}\t$answers{$pair}{allele1}{names}\t$answers{$pair}{allele2}{names}\t$answers{$pair}{propReads}\t$answers{$pair}{propSig}\t$answers{$pair}{cor}\t$answers{$pair}{error}\t$answers{$pair}{pp_score}\t$answers{$pair}{final_score}\t$comment1\t$comment2\n";
			$answer_cnt++;	
		}else{
			last;	
		}	
	}
	if($answer eq "" || $max_score == 0){
		$answer = "HLA\tHLA\t$opts{gene}\t$opts{gene}\t-\t-\t-\t-\t-\t-\tNo call made\t\n";
	}
	return($answer);
}

sub createProfileAndCalculateSignal{
	my $gene = shift;
	$gene = shift if ($gene eq "HLAPredict");
	my $allele1 = shift;
	my $allele2 = shift;
	my $name1 = shift;
	my $name2 = shift;
	my $replace1 = shift;
	my $replace2 = shift;
	my $seq1 = shift;
	my $seq2 = shift;
	my $total_cnts = shift;
	my $offset = shift;
	my $sig = shift;
	my $kmer_keys = shift;
	my $sig_kmers = shift;
	my $allele_cnts = shift; 
	my $reference = shift;

	load "$opts{scripts_dir}/modules/SimulateReads.pm";
	load "$opts{scripts_dir}/modules/RunKraken.pm";
	load "$opts{scripts_dir}/modules/ReadCounter.pm";
	my $rlogfile = "$opts{out_dir}/$gene.refinement.log";
	my $rlog;
	open($rlog, ">>$rlogfile");
	if(! (-e "$opts{out_dir}/refinement_files" && -d "$opts{out_dir}/refinement_files")){
		print $rlog "Making refinement directory (mkdir $opts{out_dir}/refinement_files)...";
		if(! mkdir "$opts{out_dir}/refinement_files"){
			print $rlog "Fatal Error: Error creating directory $opts{out_dir}/refinement_files\n";
			exit 1;
		}else{
			print $rlog "DONE\n";
		}
	}
	if(! (-e "$opts{out_dir}/refinement_files/simulated" && -d "$opts{out_dir}/refinement_files/simulated")){
		print $rlog "Making refinement simulation directory (mkdir $opts{out_dir}/refinement_files/simulated)...";
		if(! mkdir "$opts{out_dir}/refinement_files/simulated"){
			print $rlog "Fatal Error: Error creating directory $opts{out_dir}/refinement_files/simulated\n";
			exit 1;
		}else{
			print $rlog "DONE\n";
		}
	}
						
	open(OUT,">$opts{out_dir}/refinement_files/refined_reference.$gene.fa"); 	
	$opts{kraken_db} =~ s/\/$//;
	my $database_name = $opts{kraken_db};
	$database_name =~ s/.*\///;
	$opts{kraken_db} =~ m/(.*)\/.+$/;
	my $database_dir = $1;

	##Set simulation options 
	SimulateReads->setSimulationOptions("$opts{out_dir}/refinement_files/refined_alleles.fa", $opts{sim_num_reads}, $opts{sim_read_length},$opts{sim_max_insert}, $opts{sim_scale}, $opts{sim_shape}, $opts{sim_seed}, $opts{threads},$opts{scripts_dir},"simulated");
	
	if($replace1 == 1){
		print OUT ">$allele1 $name1\n$seq1\n";
		##Simulate reads from the new allele
		SimulateReads->makeSim($allele1, $seq1, "$opts{out_dir}/refinement_files/simulated/simulated");
		
		##Filter reads with kraken
		RunKraken->filterReads($database_dir, $database_name, "$opts{out_dir}/refinement_files/simulated/simulated.$allele1",$opts{threads}, "$opts{out_dir}/refinement_files/simulated/simulated.${allele1}_1.fastq","$opts{out_dir}/refinement_files/simulated/simulated.${allele1}_2.fastq",$opts{kraken_path},"", $rlogfile);
		
		##Count reads
		ReadCounter->count_reads("$opts{out_dir}/refinement_files/simulated/simulated.$allele1.${gene}_1.fq","$opts{out_dir}/refinement_files/simulated/simulated.$allele1.${gene}_1.uniq.cnts",$rlogfile); 
		ReadCounter->count_reads("$opts{out_dir}/refinement_files/simulated/simulated.$allele1.${gene}_2.fq","$opts{out_dir}/refinement_files/simulated/simulated.$allele1.${gene}_2.uniq.cnts",$rlogfile); 
		
		##Add reads to the profile
		##Putting read 2 first to make order of read in determineProfile
		addReadsToProfile("$opts{out_dir}/refinement_files/simulated/simulated.$allele1.${gene}_1.uniq.cnts",$allele1,1,$sig_kmers,$kmer_keys,\%combined_kmers,$sig,$allele_cnts);
		addReadsToProfile("$opts{out_dir}/refinement_files/simulated/simulated.$allele1.${gene}_2.uniq.cnts",$allele1,0,$sig_kmers,$kmer_keys,\%combined_kmers,$sig,$allele_cnts);
	}
	
	if($replace2 == 1){
		print OUT ">$allele2 $name2\n$seq2\n";			
		##Simulate reads from the new allele
		SimulateReads->makeSim($allele2, $seq2, "$opts{out_dir}/refinement_files/simulated/simulated");
		
		##Filter reads with kraken
		RunKraken->filterReads($database_dir, $database_name, "$opts{out_dir}/refinement_files/simulated/simulated.$allele2",$opts{threads}, "$opts{out_dir}/refinement_files/simulated/simulated.${allele2}_1.fastq","$opts{out_dir}/refinement_files/simulated/simulated.${allele2}_2.fastq",$opts{kraken_path},"", $rlogfile);
		
		##Count reads
		ReadCounter->count_reads("$opts{out_dir}/refinement_files/simulated/simulated.$allele2.${gene}_1.fq","$opts{out_dir}/refinement_files/simulated/simulated.$allele2.${gene}_1.uniq.cnts",$rlogfile); 
		ReadCounter->count_reads("$opts{out_dir}/refinement_files/simulated/simulated.$allele2.${gene}_2.fq","$opts{out_dir}/refinement_files/simulated/simulated.$allele2.${gene}_2.uniq.cnts",$rlogfile); 
		
		##Add reads to the profile
		##Putting read 2 first to make order of read in determineProfile
		addReadsToProfile("$opts{out_dir}/refinement_files/simulated/simulated.$allele2.${gene}_1.uniq.cnts",$allele2,1,$sig_kmers,$kmer_keys,\%combined_kmers,$sig,$allele_cnts);
		addReadsToProfile("$opts{out_dir}/refinement_files/simulated/simulated.$allele2.${gene}_2.uniq.cnts",$allele2,0,$sig_kmers,$kmer_keys,\%combined_kmers,$sig,$allele_cnts);
	
	}

	if(open(IN, $reference)){
		while(<IN>){
			print OUT $_;
		}
		close(IN);
	}
		
	close(OUT);
	##Calculate the correlation, proportion of reads, proportion of signal for the new allele pair
	my ($cor, $propReads, $propSignal) = getSSE($allele1, $allele2, $total_cnts, ($allele_cnts->{$allele1} + $allele_cnts->{$allele2}), $offset, $sig);
	return ($cor, $propReads, $propSignal,"$opts{out_dir}/refinement_files/refined_reference.${gene}.fa");
}

sub addReadsToProfile{
	my $file = shift;
	$file = shift if ($file eq "HLAPredict");
	my $allele = shift;
	my $revcomp = shift;
	my $sig_kmers = shift;
	my $kmer_keys = shift;
	my $combined_kmers = shift;
	my $sig = shift;
	my $allele_cnts = shift;
	
	my $kmer_cnt = scalar keys %{$sig_kmers};
	open(IN, $file);
	while(<IN>){
		chomp;
		my @parts = split /\t/, $_;
		my $key = "";
		$parts[0] = SequenceFunctions->revcomp($parts[0]) if ($revcomp == 1);
		if(defined $sig_kmers->{$parts[0]}){
		$key = $sig_kmers->{$parts[0]};
	}else{
		$key = $kmer_cnt;
		$kmer_keys->{$key}=$parts[0];
		$sig_kmers->{$parts[0]}=$key;
		$combined_kmers->{$parts[0]}=$key;
		$kmer_cnt++;
	}
	$sig->{$allele}{$key}+=$parts[1];
	$allele_cnts->{$allele}+=$parts[1];
}
close(IN);	
}

sub sumMatrix{
	my $pp_results_ref = shift;
	$pp_results_ref = shift if($pp_results_ref eq "HLAPredict");
	my $finalist_ref = shift;
	my @finalists = @$finalist_ref;
	my $noise_weight = shift;
	my %matrix_sums;
	my @pp_results = @$pp_results_ref;
	for(my $i = 0; $i<=$#pp_results; $i++){
		for(my $j = 0; $j<=$#pp_results; $j++){
			if($i < $j){
				$matrix_sums{$finalists[$i]}{row}+=($noise_weight+$pp_results[$i][$j]);
				$matrix_sums{$finalists[$i]}{col}+=($noise_weight+$pp_results[$j][$i]);
				$matrix_sums{$finalists[$j]}{row}+=($noise_weight+$pp_results[$j][$i]);
				$matrix_sums{$finalists[$j]}{col}+=($noise_weight+$pp_results[$i][$j]);
			}
		}
	}
	return(\%matrix_sums);
}

sub readUnique12{
	my $prefix = shift;
	$prefix = shift if ($prefix eq "HLAPredict");
	my $klen = shift;
	my $sig_kmers = shift;
	my $sig = shift;
	my $read_total = 0;
	my $total_cnts = 0;
	my $mx;
	
	for(my $read=1; $read<=2; $read++) {
		my $fn = "${prefix}_$read.uniq.cnts";
		open(IN, $fn) || die "Cannot open unique file: $fn\n";
		READ:while(<IN>){
			chomp;
			my @parts = split /\t/, $_;
			$parts[0] = SequenceFunctions->revcomp($parts[0]) if ($read == 1);
			my $key = $parts[0];
			$mx = length($key) - $klen + 1;
			my $value = $parts[1];
			my $key_rc = SequenceFunctions->revcomp($key);				
			for(my $i=0; $i<$mx; $i++) {
				my $substr = substr ($key, $i, $klen);
				my $substr_rc = substr ($key_rc, $i, $klen);
				if( defined $sig_kmers->{$substr} ) {
					if( defined $sig->{cnt}{$substr} ) {
						$sig->{cnt}{$substr} += $value;
					} else {
						$sig->{cnt}{$substr} = $value;
					}
					$total_cnts += $value;
				}elsif(defined $sig_kmers->{$substr_rc}){
					if( defined $sig->{cnt}{$substr_rc} ) {
						$sig->{cnt}{$substr_rc} += $value;
					} else {
						$sig->{cnt}{$substr_rc} = $value;
					}
					$total_cnts += $value;
				}elsif(defined $sig->{cnt}{$substr}){
					$sig->{cnt}{$substr}+=$value;
					$total_cnts += $value;
				}elsif(defined $sig->{cnt}{$substr_rc}){
					$sig->{cnt}{$substr_rc}+=$value;
					$total_cnts += $value;
				}else{
					$sig->{cnt}{$substr}=$value;
					$total_cnts += $value;
				}	
			}
			$read_total+=$value;
		}
		close(IN);	
	}
	$read_total = $read_total/2;
	return($read_total, $total_cnts,$mx);
}

sub findCombinedKmers{
	my $sig = shift;
	$sig = shift if ($sig eq "HLAPredict");
	my $kmer_keys = shift;
	my $allele_ref = shift;
	my $thresh = shift;
	my @topAlleles = @{$allele_ref};
	%combined_kmers=();	
	foreach my $key (keys %{$sig->{cnt}}){
			$combined_kmers{$key}= "" if($sig->{cnt}{$key} >= $thresh);
	}
	foreach my $allele (@topAlleles){
		my $key = "";
		foreach $key (keys %{$sig->{$allele}}){
			if($key ne ""){
				$combined_kmers{$kmer_keys->{$key}}=$key;
			}
		}
	}
	return \%combined_kmers;
}

sub getSSE{
	my $cand1 = shift;
	$cand1 = shift if ($cand1 eq "HLAPredict");
	my $cand2 = shift;
	my $obs_total = shift;
	my $pred_total = shift;
	my $offset = shift;
	my $profile = shift;		
	if(!defined $cand1 || !defined $cand2){
		print STDERR "$opts{gene} Cand1 $cand1 or Cand2 $cand2 not defined\n";
		return (0,0,0);
	}
	
	
	my $inpred_cnt = 0;
	my $inobs_sum = 0;
	my $inboth_cnt = 0;
	my $inboth_sum = 0;

	my @obs_cnt=();
	my @pred_cnt=();

	my $cnt_thresh = $opts{cntThresh};	
	foreach my $k (keys %combined_kmers){
		##Determines proportions
		my $ocnt = $profile->{cnt}{$k} || 0;
		my $pred1 = $profile->{$cand1}{$combined_kmers{$k}} || 0;
		my $pred2 = $profile->{$cand2}{$combined_kmers{$k}} || 0;
		my $pcnt = $pred1 + $pred2;

		if($pcnt > 0 && $ocnt >= $cnt_thresh){
			$inboth_cnt++;
			$inboth_sum+=$ocnt;
			$inpred_cnt++;
			$inobs_sum+=$ocnt;
		}elsif($pcnt > 0){
			$inpred_cnt++;
		}elsif($ocnt >= $cnt_thresh){
			$inobs_sum+=$ocnt;
		}
		##Create arrays for calculating correlation
		push (@obs_cnt, log($ocnt/$obs_total + $offset));
		push (@pred_cnt, log($pcnt/$pred_total + $offset));
	}
	my $cor = 0+ Statistics::Basic::correlation(Statistics::Basic::vector(@pred_cnt),Statistics::Basic::vector(@obs_cnt));
	
	my $propReads = 0;
	my $propSig = 0;
	if($inobs_sum != 0){
		$propReads = $inboth_sum/$inobs_sum;	
	}
	if($inpred_cnt != 0){
		$propSig = $inboth_cnt/$inpred_cnt;	
	}
	return($cor, $propReads, $propSig);
}

__PACKAGE__->run_Commandline() unless caller;
