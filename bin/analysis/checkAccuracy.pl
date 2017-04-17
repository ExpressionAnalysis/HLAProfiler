#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;

(my $SCRIPT_NAME = $0) =~ s/.*\///;
my $version = "1.0";
my $creation_date = "2 Feb 2017";
my $last_updated = "2 Feb 2017";

my $usage="\n$SCRIPT_NAME v$version\n".
	  "\nDESCRIPTIONs:\n" .
	  "This script accepts output from various HLA callers and calculates concordance/accuracy with the supplied truth set.\n\n" . 
	  "\nUSAGE:\n" .
	  "perl $SCRIPT_NAME <options>\n" .
	  "\nRequired options:\n" .
	  "-truth|t\t\tTab delimited file containing the true calls for samples. This file contains the sample name in the first column and alleles starting in third column. A header is required and should contain only the gene name for allele columns.(required)\n" .
	  "-predictions|p\t\tTab delimited, 2 column file containing the location of predictions for each sample. The first column should be the sample name.(required)\n" .
	  "\nGeneral options:\n" .
	  "-caller|c\t\tThe caller used to generate the prediction files (default:'HLAProfiler')\n" .
	  "\tPossible options\n" .
	  "\tHLAForest\n" .
	  "\tHLAMiner`\n" .
	  "\tHLAProfiler\n" .
	  "\tOptiType`\n" .
	  "\tPHLAT`\n" .
	  "\tseq2HLA\n" .
	  "-ambiguous_alleles|a\tIf there is ambiguity in the truth caused by sequencing of a subset of exons, supply the path to the IMGT ambiguous allele file here. This will check predictions possible truth alleles. (default:'')\n" .
	  "-log|l\t\tLog file(default:STDERR)\n" .
	  "-allele_subset|s\tPath to a tab-delimeted file containing a subset of alleles to count, with sample in column1 and the truth allele in column2 (default:'')\n" .
          "-help|h\t\tDisplays this message\n" .
          "\nAUTHORS:\n" .
	  "Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	  "\nCREATED:\n$creation_date\n" .
	  "\nLAST UPDATED:\n$last_updated\n" .
	  "\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
	  "\n";

my %opts = (truth=>"",allele_subset=>"", predictions=>"", ambiguous_alleles=>"", caller=>"HLAProfiler");
GetOptions(\%opts, qw(help|h truth|t=s allele_subset|s=s predictions|p=s ambiguous_alleles|a=s caller|c=s));

if($opts{help}){
	print "$usage";
	exit;	
}elsif($opts{truth} eq ""){
	print STDERR "Please specify a truth file.\n$usage";
	exit 1;
}elsif($opts{predictions} eq ""){
	print STDERR "Please specify a predictions file.\n$usage";
	exit 1;
}

my $log;
if($opts{log}){
	open($log, ">$opts{log}");
}else{
	$log=*STDERR;
}

my %truth;
my %ambiguous_alleles;
my %stats;
my %subset=();

if($opts{allele_subset} ne ""){
	my $sfh;
	open($sfh, "$opts{allele_subset}");
	while(<$sfh>){
		chomp;
		my @cols = split /\t/, $_;
		$subset{$cols[0]}{$cols[1]}=1;		
	}
	close($sfh);
}

populateAmbiguousAllele($opts{ambiguous_alleles}, \%ambiguous_alleles) if ($opts{ambiguous_alleles} ne '');
populateTruth($opts{truth}, \%truth, \%ambiguous_alleles);
my $predictions = "";	
if($opts{caller} eq "HLAProfiler"){
	$predictions = loadHLAProfilerPredictions($opts{predictions}, \%ambiguous_alleles);
}elsif($opts{caller} eq "HLAMiner"){
	$predictions = loadHLAMinerPredictions($opts{predictions}, \%ambiguous_alleles);
}elsif($opts{caller} eq "HLAForest"){
	$predictions = loadHLAForestPredictions($opts{predictions}, \%ambiguous_alleles);
}elsif($opts{caller} eq "OptiType"){
	$predictions = loadOptiTypePredictions($opts{predictions}, \%ambiguous_alleles);
}elsif($opts{caller} eq "PHLAT"){
	$predictions = loadPhlatPredictions($opts{predictions}, \%ambiguous_alleles);
}elsif($opts{caller} eq "seq2HLA"){
	$predictions = loadSeq2HLAPredictions($opts{predictions}, \%ambiguous_alleles);
}else{
	print STDERR "-c $opts{caller} is not a valid HLATyping algorithm option.$usage\n";
	exit 1;
}

my $wrong = checkConcordance($predictions, \%truth, \%stats);
summarizeStats(\%stats);
summarizeSubsetStats(\%stats) if ($opts{allele_subset} ne "");
printWrongAnswers($wrong);

sub populateAmbiguousAllele{
	##Some alleles are identical at exon 2 and 3 (commonly used for orthogonal typing). IMGT groups these into groups. This will allow for comparisons to determine if the detected allele is from the same group as the truth even though the call at 4-digit precision is incorrect.
	my $ambiguous_allele_file = shift;
	my $ambiguous = shift;
	
	my $afh;
	open($afh,  $ambiguous_allele_file);
	while(<$afh>){
		chomp;
		my @cols = split /[ \t]/, $_;
		$cols[0] =~ s/G$//;
		for(my $i = 1; $i <=$#cols; $i++){
			$ambiguous->{$cols[$i]}=$cols[0];
		}
	}
	
}

sub loadHLAProfilerPredictions{
	my $file = shift;
	my $ambiguous = shift;
	my $fh;
	my %predictions = ();
	open($fh, "$file") || (print $log "Cannot open predictions file $file\n" && exit 1);
	while(<$fh>){
		chomp;
		my @cols = split /\t/, $_;
		my $pfh;
		open($pfh, $cols[1]) || (print $log "Cannot open predictions file for $cols[0] $cols[1]\n" && exit 1);
		my $line = <$pfh>;
		my %scores = ();
		LINE:while(my $line = <$pfh>){
			chomp $line;
			my @parts = split /\t/, $line;
			#Remove added parts for novel or updated alleles
			my $allele_string1 = $parts[2];
			my $allele_string2 = $parts[3];
			$parts[2] =~ s/_.*//;
			$parts[3] =~ s/_.*//;
			next LINE if($parts[2] eq "" || $parts[3] eq "");
				
			#Check for ambiguous alleles
			if(defined $ambiguous->{$parts[2]}){
				$parts[2] = $ambiguous->{$parts[2]};
			}
			if(defined $ambiguous->{$parts[3]}){
				$parts[3] = $ambiguous->{$parts[3]};
			}
	
			my @a1_parts = split /[\*:]/, $parts[2];				
			my @a2_parts = split /[\*:]/, $parts[3];
			
			my $gene1 = shift @a1_parts;
			my $gene2 = shift @a2_parts;

			#In case the allele predictions become out of order (not sorted with most highest score first).THIS HAS HAPPENED. These next two lines will wipe the data for the gene clean and allow repeat the process for samples of hte same gene with a higher score
			next LINE if(defined $scores{$gene1} && $parts[9] < $scores{$gene1});
			$scores{$gene1}=$parts[9];
			$predictions{$cols[0]}{$gene1}=(); 
			$predictions{$cols[0]}{$gene1}{string}{1}=$allele_string1;		
			$predictions{$cols[0]}{$gene2}{string}{2}=$allele_string2;		
			if($#a1_parts == 0){ ##2 digit
				$predictions{$cols[0]}{$gene1}{2}{1}{$a1_parts[0]}=1;
				$predictions{$cols[0]}{$gene1}{digits}{2}{1} = 1;
				$predictions{$cols[0]}{$gene1}{exact}{1}{$a1_parts[0]}=1;
			}elsif($#a1_parts == 1){ ##4 digit
				$predictions{$cols[0]}{$gene1}{digits}{2}{1} = 1;
				$predictions{$cols[0]}{$gene1}{digits}{4}{1} = 1;
				$predictions{$cols[0]}{$gene1}{2}{1}{"$a1_parts[0]"}=1;
				$predictions{$cols[0]}{$gene1}{4}{1}{"$a1_parts[0]:$a1_parts[1]"}=1;
				$predictions{$cols[0]}{$gene1}{exact}{1}{"$a1_parts[0]:$a1_parts[1]"}=1;
			}elsif($#a1_parts == 2){ ##6 digit
				$predictions{$cols[0]}{$gene1}{digits}{2}{1} = 1;
				$predictions{$cols[0]}{$gene1}{digits}{4}{1} = 1;
				$predictions{$cols[0]}{$gene1}{digits}{6}{1} = 1;
				$predictions{$cols[0]}{$gene1}{2}{1}{"$a1_parts[0]"}=1;
				$predictions{$cols[0]}{$gene1}{4}{1}{"$a1_parts[0]:$a1_parts[1]"}=1;
				$predictions{$cols[0]}{$gene1}{6}{1}{"$a1_parts[0]:$a1_parts[1]:$a1_parts[2]"}=1;
				$predictions{$cols[0]}{$gene1}{exact}{1}{"$a1_parts[0]:$a1_parts[1]:$a1_parts[2]"}=1;
			}elsif($#a1_parts == 3){ ##8 digit
				$predictions{$cols[0]}{$gene1}{digits}{2}{1} = 1;
				$predictions{$cols[0]}{$gene1}{digits}{4}{1} = 1;
				$predictions{$cols[0]}{$gene1}{digits}{6}{1} = 1;
				$predictions{$cols[0]}{$gene1}{digits}{8}{1} = 1;
				$predictions{$cols[0]}{$gene1}{2}{1}{"$a1_parts[0]"}=1;
				$predictions{$cols[0]}{$gene1}{4}{1}{"$a1_parts[0]:$a1_parts[1]"}=1;
				$predictions{$cols[0]}{$gene1}{6}{1}{"$a1_parts[0]:$a1_parts[1]:$a1_parts[2]"}=1;
				$predictions{$cols[0]}{$gene1}{8}{1}{"$a1_parts[0]:$a1_parts[1]:$a1_parts[2]:$a1_parts[3]"}=1;
				$predictions{$cols[0]}{$gene1}{exact}{1}{"$a1_parts[0]:$a1_parts[1]:$a1_parts[2]:$a1_parts[3]"}=1;
			}
			if($#a2_parts == 0){ ##2 digit
				$predictions{$cols[0]}{$gene2}{2}{2}{$a2_parts[0]}=1;
				$predictions{$cols[0]}{$gene2}{digits}{2}{2} = 1;
				$predictions{$cols[0]}{$gene2}{exact}{2}{$a2_parts[0]}=1;
			}elsif($#a2_parts == 1){ ##4 digit
				$predictions{$cols[0]}{$gene2}{digits}{4}{2} = 1;
				$predictions{$cols[0]}{$gene2}{2}{2}{"$a2_parts[0]"}=1;
				$predictions{$cols[0]}{$gene2}{4}{2}{"$a2_parts[0]:$a2_parts[1]"}=1;
				$predictions{$cols[0]}{$gene2}{exact}{2}{"$a2_parts[0]:$a2_parts[1]"}=1;
			}elsif($#a2_parts == 2){ ##6 digit
				$predictions{$cols[0]}{$gene2}{digits}{6}{2} = 1;
				$predictions{$cols[0]}{$gene2}{2}{2}{"$a2_parts[0]"}=1;
				$predictions{$cols[0]}{$gene2}{4}{2}{"$a2_parts[0]:$a2_parts[1]"}=1;
				$predictions{$cols[0]}{$gene2}{6}{2}{"$a2_parts[0]:$a2_parts[1]:$a2_parts[2]"}=1;
				$predictions{$cols[0]}{$gene2}{exact}{2}{"$a2_parts[0]:$a2_parts[1]:$a2_parts[2]"}=1;
			}elsif($#a2_parts == 3){ ##8 digit
				$predictions{$cols[0]}{$gene2}{digits}{8}{2} = 1;
				$predictions{$cols[0]}{$gene2}{2}{2}{"$a2_parts[0]"}=1;
				$predictions{$cols[0]}{$gene2}{4}{2}{"$a2_parts[0]:$a2_parts[1]"}=1;
				$predictions{$cols[0]}{$gene2}{6}{2}{"$a2_parts[0]:$a2_parts[1]:$a2_parts[2]"}=1;
				$predictions{$cols[0]}{$gene2}{8}{2}{"$a2_parts[0]:$a2_parts[1]:$a2_parts[2]:$a2_parts[3]"}=1;
				$predictions{$cols[0]}{$gene2}{exact}{2}{"$a2_parts[0]:$a2_parts[1]:$a2_parts[2]:$a2_parts[3]"}=1;
			}
							
		}
		close($pfh);
	}
	close($fh);
	return (\%predictions);
}
sub loadOptiTypePredictions{
	my $file = shift;
	my $ambiguous = shift;
	my $fh;
	my %predictions = ();
	open($fh, "$file") || (print $log "Cannot open predictions file $file\n" && exit 1);
	while(<$fh>){
		chomp;
		my @cols = split /\t/, $_;
		my $pfh;
		open($pfh, $cols[1]) || (print $log "Cannot open predictions file for $cols[0] $cols[1]\n" && exit 1);
		my $line = <$pfh>;
		LINE:while(my $line = <$pfh>){
			chomp $line;
			my @parts = split /\t/, $line;
			for(my $i =1; $i<=6; $i++){
				#Check for ambiguous alleles
				if(defined $ambiguous->{$parts[$i]}){
					$parts[2] = $ambiguous->{$parts[$i]};
				}
	
				my @a_parts = split /[\*:]/, $parts[$i];				
				
				my $gene = shift @a_parts;
				my $num = 1;
				$num = 2 if (exists $predictions{$cols[0]}{$gene}{exact}{1});
				$predictions{$cols[0]}{$gene}{string}{$num}=$parts[$i];	
				if($#a_parts == 0){ ##2 digit
					$predictions{$cols[0]}{$gene}{2}{$num}{$a_parts[0]}=1;
					$predictions{$cols[0]}{$gene}{digits}{2}{$num} = 1;
					$predictions{$cols[0]}{$gene}{exact}{$num}{$a_parts[0]}=1;
				}elsif($#a_parts == 1){ ##4 digit
					$predictions{$cols[0]}{$gene}{digits}{2}{$num} = 1;
					$predictions{$cols[0]}{$gene}{digits}{4}{$num} = 1;
					$predictions{$cols[0]}{$gene}{2}{$num}{"$a_parts[0]"}=1;
					$predictions{$cols[0]}{$gene}{4}{$num}{"$a_parts[0]:$a_parts[1]"}=1;
					$predictions{$cols[0]}{$gene}{exact}{$num}{"$a_parts[0]:$a_parts[1]"}=1;
				}elsif($#a_parts == 2){ ##6 digit
					$predictions{$cols[0]}{$gene}{digits}{2}{$num} = 1;
					$predictions{$cols[0]}{$gene}{digits}{4}{$num} = 1;
					$predictions{$cols[0]}{$gene}{digits}{6}{$num} = 1;
					$predictions{$cols[0]}{$gene}{2}{$num}{"$a_parts[0]"}=1;
					$predictions{$cols[0]}{$gene}{4}{$num}{"$a_parts[0]:$a_parts[1]"}=1;
					$predictions{$cols[0]}{$gene}{6}{$num}{"$a_parts[0]:$a_parts[1]:$a_parts[2]"}=1;
					$predictions{$cols[0]}{$gene}{exact}{$num}{"$a_parts[0]:$a_parts[1]:$a_parts[2]"}=1;
				}elsif($#a_parts == 3){ ##8 digit
					$predictions{$cols[0]}{$gene}{digits}{2}{$num} = 1;
					$predictions{$cols[0]}{$gene}{digits}{4}{$num} = 1;
					$predictions{$cols[0]}{$gene}{digits}{6}{$num} = 1;
					$predictions{$cols[0]}{$gene}{digits}{8}{$num} = 1;
					$predictions{$cols[0]}{$gene}{2}{$num}{"$a_parts[0]"}=1;
					$predictions{$cols[0]}{$gene}{4}{$num}{"$a_parts[0]:$a_parts[1]"}=1;
					$predictions{$cols[0]}{$gene}{6}{$num}{"$a_parts[0]:$a_parts[1]:$a_parts[2]"}=1;
					$predictions{$cols[0]}{$gene}{8}{$num}{"$a_parts[0]:$a_parts[1]:$a_parts[2]:$a_parts[3]"}=1;
					$predictions{$cols[0]}{$gene}{exact}{$num}{"$a_parts[0]:$a_parts[1]:$a_parts[2]:$a_parts[3]"}=1;
				}
			}
		}
		close($pfh);
	}
	close($fh);
	return(\%predictions);
}
sub loadHLAForestPredictions {
	my $file = shift;
	my $ambiguous = shift;
	my $fh;
	my %predictions = ();
	open($fh, "$file") || (print $log "Cannot open predictions file $file\n" && exit 1);
	while(<$fh>){
		chomp;
		my @cols = split /\t/, $_;
		my $pfh;
		open($pfh, $cols[1]) || (print $log "Cannot open predictions file for $cols[0] $cols[1]\n" && exit 1);
		LINE:while(my $line = <$pfh>){
			chomp ($line);
			my @parts = split /\t/, $line;
			$parts[$#parts] =~ s/ROOT://;	
			$parts[$#parts] =~ s/:/\*/;	
			if(defined $ambiguous->{$parts[$#parts]}){
				$parts[$#parts] = $ambiguous->{$parts[$#parts]};
			}

			my @a_parts = split /[\*:]/, $parts[$#parts];				
			
			my $gene = shift @a_parts;
			my $num = 1;
			$num = 2 if (exists $predictions{$cols[0]}{$gene}{exact}{1});
			$predictions{$cols[0]}{$gene}{string}{$num}=$parts[$#parts];		
			if($#a_parts == 0){ ##2 digit
				$predictions{$cols[0]}{$gene}{2}{$num}{$a_parts[0]}=1;
				$predictions{$cols[0]}{$gene}{digits}{2}{$num} = 1;
				$predictions{$cols[0]}{$gene}{exact}{$num}{$a_parts[0]}=1;
			}elsif($#a_parts == 1){ ##4 digit
				$predictions{$cols[0]}{$gene}{digits}{2}{$num} = 1;
				$predictions{$cols[0]}{$gene}{digits}{4}{$num} = 1;
				$predictions{$cols[0]}{$gene}{2}{$num}{"$a_parts[0]"}=1;
				$predictions{$cols[0]}{$gene}{4}{$num}{"$a_parts[0]:$a_parts[1]"}=1;
				$predictions{$cols[0]}{$gene}{exact}{$num}{"$a_parts[0]:$a_parts[1]"}=1;
			}elsif($#a_parts == 2){ ##6 digit
				$predictions{$cols[0]}{$gene}{digits}{2}{$num} = 1;
				$predictions{$cols[0]}{$gene}{digits}{4}{$num} = 1;
				$predictions{$cols[0]}{$gene}{digits}{6}{$num} = 1;
				$predictions{$cols[0]}{$gene}{2}{$num}{"$a_parts[0]"}=1;
				$predictions{$cols[0]}{$gene}{4}{$num}{"$a_parts[0]:$a_parts[1]"}=1;
				$predictions{$cols[0]}{$gene}{6}{$num}{"$a_parts[0]:$a_parts[1]:$a_parts[2]"}=1;
				$predictions{$cols[0]}{$gene}{exact}{$num}{"$a_parts[0]:$a_parts[1]:$a_parts[2]"}=1;
			}elsif($#a_parts == 3){ ##8 digit
				$predictions{$cols[0]}{$gene}{digits}{2}{$num} = 1;
				$predictions{$cols[0]}{$gene}{digits}{4}{$num} = 1;
				$predictions{$cols[0]}{$gene}{digits}{6}{$num} = 1;
				$predictions{$cols[0]}{$gene}{digits}{8}{$num} = 1;
				$predictions{$cols[0]}{$gene}{2}{$num}{"$a_parts[0]"}=1;
				$predictions{$cols[0]}{$gene}{4}{$num}{"$a_parts[0]:$a_parts[1]"}=1;
				$predictions{$cols[0]}{$gene}{6}{$num}{"$a_parts[0]:$a_parts[1]:$a_parts[2]"}=1;
				$predictions{$cols[0]}{$gene}{8}{$num}{"$a_parts[0]:$a_parts[1]:$a_parts[2]:$a_parts[3]"}=1;
				$predictions{$cols[0]}{$gene}{exact}{$num}{"$a_parts[0]:$a_parts[1]:$a_parts[2]:$a_parts[3]"}=1;
			}
		}
		close($pfh);
	}
	close($fh);
	return(\%predictions);
}
sub loadHLAMinerPredictions {
	my $file = shift;
	my $ambiguous = shift;
	my $fh;
	my %predictions = ();
	open($fh, "$file") || (print $log "Cannot open predictions file $file\n" && exit 1);
	while(<$fh>){
		chomp;
		my @cols = split /\t/, $_;
		my $pfh;
		open($pfh, $cols[1]) || (print $log "Cannot open predictions file for $cols[0] $cols[1]\n" && exit 1);
		my $gene = "";
		my $num = 0;
		my $min_score = 0;
		LINE:while(my $line = <$pfh>){
			chomp ($line);
			if($line =~ m/HLA-(.*)$/){
				$gene = $1;
				$min_score = 0;
				$predictions{$cols[0]}{$gene}{string}{1}="NC";
				$predictions{$cols[0]}{$gene}{string}{2}="NC";
			}elsif($line =~ m/Prediction #(\d)/){
				$num = $1;
				$min_score = 0;
			}elsif($line =~ m/.*\*\d+:\d+/){		
				$line =~ s/^[ \t]*//;
				my @parts = split /,/, $line;
				next LINE if ($parts[1] < $min_score);
				$parts[0]=~s/P$//;
				$min_score = $parts[1];
				if(defined $ambiguous->{$parts[0]}){
					$parts[0] = $ambiguous->{$parts[0]};
				}
				my @a_parts = split /[\*:]/, $parts[0];				
				shift @a_parts;
			
				if($predictions{$cols[0]}{$gene}{string}{$num} ne "NC"){
					$predictions{$cols[0]}{$gene}{string}{$num}.="/$parts[0]";	
				}else{
					$predictions{$cols[0]}{$gene}{string}{$num}=$parts[0];	
				}		
				if($#a_parts == 0){ ##2 digit
					$predictions{$cols[0]}{$gene}{2}{$num}{$a_parts[0]}=1;
					$predictions{$cols[0]}{$gene}{digits}{2}{$num} = 1;
					$predictions{$cols[0]}{$gene}{exact}{$num}{$a_parts[0]}=1;
				}elsif($#a_parts == 1){ ##4 digit
					$predictions{$cols[0]}{$gene}{digits}{2}{$num} = 1;
					$predictions{$cols[0]}{$gene}{digits}{4}{$num} = 1;
					$predictions{$cols[0]}{$gene}{2}{$num}{"$a_parts[0]"}=1;
					$predictions{$cols[0]}{$gene}{4}{$num}{"$a_parts[0]:$a_parts[1]"}=1;
					$predictions{$cols[0]}{$gene}{exact}{$num}{"$a_parts[0]:$a_parts[1]"}=1;
				}elsif($#a_parts == 2){ ##6 digit
					$predictions{$cols[0]}{$gene}{digits}{2}{$num} = 1;
					$predictions{$cols[0]}{$gene}{digits}{4}{$num} = 1;
					$predictions{$cols[0]}{$gene}{digits}{6}{$num} = 1;
					$predictions{$cols[0]}{$gene}{2}{$num}{"$a_parts[0]"}=1;
					$predictions{$cols[0]}{$gene}{4}{$num}{"$a_parts[0]:$a_parts[1]"}=1;
					$predictions{$cols[0]}{$gene}{6}{$num}{"$a_parts[0]:$a_parts[1]:$a_parts[2]"}=1;
					$predictions{$cols[0]}{$gene}{exact}{$num}{"$a_parts[0]:$a_parts[1]:$a_parts[2]"}=1;
				}elsif($#a_parts == 3){ ##8 digit
					$predictions{$cols[0]}{$gene}{digits}{2}{$num} = 1;
					$predictions{$cols[0]}{$gene}{digits}{4}{$num} = 1;
					$predictions{$cols[0]}{$gene}{digits}{6}{$num} = 1;
					$predictions{$cols[0]}{$gene}{digits}{8}{$num} = 1;
					$predictions{$cols[0]}{$gene}{2}{$num}{"$a_parts[0]"}=1;
					$predictions{$cols[0]}{$gene}{4}{$num}{"$a_parts[0]:$a_parts[1]"}=1;
					$predictions{$cols[0]}{$gene}{6}{$num}{"$a_parts[0]:$a_parts[1]:$a_parts[2]"}=1;
					$predictions{$cols[0]}{$gene}{8}{$num}{"$a_parts[0]:$a_parts[1]:$a_parts[2]:$a_parts[3]"}=1;
					$predictions{$cols[0]}{$gene}{exact}{$num}{"$a_parts[0]:$a_parts[1]:$a_parts[2]:$a_parts[3]"}=1;
				}
			}
		}
		close($pfh);
	}
	close($fh);
	return(\%predictions);
}

sub loadSeq2HLAPredictions{
	my $file = shift;
	my $ambiguous = shift;
	my $fh;
	my %predictions = ();
	open($fh, "$file") || (print $log "Cannot open predictions file $file\n" && exit 1);
	while(<$fh>){
		chomp;
		my @cols = split /\t/, $_;
		my $file2 = $cols[1];
		$file2 =~ s/ClassI/ClassII/;
		my $pfh;
		open($pfh, "cat $cols[1] $file2 |") || (print $log "Cannot open predictions file for $cols[0] $cols[1] and $file2\n" && exit 1);
		LINE:while(my $line = <$pfh>){
			chomp $line;
			my @parts = split /\t/, $line;
			next LINE if($parts[1] eq "" || $parts[3] eq "");
			next LINE if(substr($line,0,1) eq "#");
			
			$parts[1] =~ s/'//;	
			$parts[3] =~ s/'//;	
			
			#Check for ambiguous alleles
			if(defined $ambiguous->{$parts[3]}){
				$parts[3] = $ambiguous->{$parts[3]};
			}
			if(defined $ambiguous->{$parts[1]}){
				$parts[1] = $ambiguous->{$parts[1]};
			}
			
			
			my @a1_parts = split /[\*:]/, $parts[1];				
			my @a2_parts = split /[\*:]/, $parts[3];
			
			my $gene1 = shift @a1_parts;
			my $gene2 = shift @a2_parts;

			$predictions{$cols[0]}{$gene1}{string}{1}=$parts[1];		
			$predictions{$cols[0]}{$gene2}{string}{2}=$parts[3];		
			if($#a1_parts == 0){ ##2 digit
				$predictions{$cols[0]}{$gene1}{2}{1}{$a1_parts[0]}=1;
				$predictions{$cols[0]}{$gene1}{digits}{2}{1} = 1;
				$predictions{$cols[0]}{$gene1}{exact}{1}{$a1_parts[0]}=1;
			}elsif($#a1_parts == 1){ ##4 digit
				$predictions{$cols[0]}{$gene1}{digits}{2}{1} = 1;
				$predictions{$cols[0]}{$gene1}{digits}{4}{1} = 1;
				$predictions{$cols[0]}{$gene1}{2}{1}{"$a1_parts[0]"}=1;
				$predictions{$cols[0]}{$gene1}{4}{1}{"$a1_parts[0]:$a1_parts[1]"}=1;
				$predictions{$cols[0]}{$gene1}{exact}{1}{"$a1_parts[0]:$a1_parts[1]"}=1;
			}elsif($#a1_parts == 2){ ##6 digit
				$predictions{$cols[0]}{$gene1}{digits}{2}{1} = 1;
				$predictions{$cols[0]}{$gene1}{digits}{4}{1} = 1;
				$predictions{$cols[0]}{$gene1}{digits}{6}{1} = 1;
				$predictions{$cols[0]}{$gene1}{2}{1}{"$a1_parts[0]"}=1;
				$predictions{$cols[0]}{$gene1}{4}{1}{"$a1_parts[0]:$a1_parts[1]"}=1;
				$predictions{$cols[0]}{$gene1}{6}{1}{"$a1_parts[0]:$a1_parts[1]:$a1_parts[2]"}=1;
				$predictions{$cols[0]}{$gene1}{exact}{1}{"$a1_parts[0]:$a1_parts[1]:$a1_parts[2]"}=1;
			}elsif($#a1_parts == 3){ ##8 digit
				$predictions{$cols[0]}{$gene1}{digits}{2}{1} = 1;
				$predictions{$cols[0]}{$gene1}{digits}{4}{1} = 1;
				$predictions{$cols[0]}{$gene1}{digits}{6}{1} = 1;
				$predictions{$cols[0]}{$gene1}{digits}{8}{1} = 1;
				$predictions{$cols[0]}{$gene1}{2}{1}{"$a1_parts[0]"}=1;
				$predictions{$cols[0]}{$gene1}{4}{1}{"$a1_parts[0]:$a1_parts[1]"}=1;
				$predictions{$cols[0]}{$gene1}{6}{1}{"$a1_parts[0]:$a1_parts[1]:$a1_parts[2]"}=1;
				$predictions{$cols[0]}{$gene1}{8}{1}{"$a1_parts[0]:$a1_parts[1]:$a1_parts[2]:$a1_parts[3]"}=1;
				$predictions{$cols[0]}{$gene1}{exact}{1}{"$a1_parts[0]:$a1_parts[1]:$a1_parts[2]:$a1_parts[3]"}=1;
			}
			if($#a2_parts == 0){ ##2 digit
				$predictions{$cols[0]}{$gene2}{2}{2}{$a2_parts[0]}=1;
				$predictions{$cols[0]}{$gene2}{digits}{2}{2} = 1;
				$predictions{$cols[0]}{$gene2}{exact}{2}{$a2_parts[0]}=1;
			}elsif($#a2_parts == 1){ ##4 digit
				$predictions{$cols[0]}{$gene2}{digits}{4}{2} = 1;
				$predictions{$cols[0]}{$gene2}{2}{2}{"$a2_parts[0]"}=1;
				$predictions{$cols[0]}{$gene2}{4}{2}{"$a2_parts[0]:$a2_parts[1]"}=1;
				$predictions{$cols[0]}{$gene2}{exact}{2}{"$a2_parts[0]:$a2_parts[1]"}=1;
			}elsif($#a2_parts == 2){ ##6 digit
				$predictions{$cols[0]}{$gene2}{digits}{6}{2} = 1;
				$predictions{$cols[0]}{$gene2}{2}{2}{"$a2_parts[0]"}=1;
				$predictions{$cols[0]}{$gene2}{4}{2}{"$a2_parts[0]:$a2_parts[1]"}=1;
				$predictions{$cols[0]}{$gene2}{6}{2}{"$a2_parts[0]:$a2_parts[1]:$a2_parts[2]"}=1;
				$predictions{$cols[0]}{$gene2}{exact}{2}{"$a2_parts[0]:$a2_parts[1]:$a2_parts[2]"}=1;
			}elsif($#a2_parts == 3){ ##8 digit
				$predictions{$cols[0]}{$gene2}{digits}{8}{2} = 1;
				$predictions{$cols[0]}{$gene2}{2}{2}{"$a2_parts[0]"}=1;
				$predictions{$cols[0]}{$gene2}{4}{2}{"$a2_parts[0]:$a2_parts[1]"}=1;
				$predictions{$cols[0]}{$gene2}{6}{2}{"$a2_parts[0]:$a2_parts[1]:$a2_parts[2]"}=1;
				$predictions{$cols[0]}{$gene2}{8}{2}{"$a2_parts[0]:$a2_parts[1]:$a2_parts[2]:$a2_parts[3]"}=1;
				$predictions{$cols[0]}{$gene2}{exact}{2}{"$a2_parts[0]:$a2_parts[1]:$a2_parts[2]:$a2_parts[3]"}=1;
			}
							
		}
		close($pfh);
	}
	close($fh);
	return (\%predictions);
}
sub loadPhlatPredictions{
	my $file = shift;
	my $ambiguous = shift;
	my $fh;
	my %predictions = ();
	open($fh, "$file") || (print $log "Cannot open predictions file $file\n" && exit 1);
	while(<$fh>){
		chomp;
		my @cols = split /\t/, $_;
		my $pfh;
		open($pfh, $cols[1]) || (print $log "Cannot open predictions file for $cols[0] $cols[1]\n" && exit 1);
		my $line = <$pfh>;
		LINE:while(my $line = <$pfh>){
			chomp $line;
			my @parts = split /\t/, $line;
			next LINE if($parts[1] eq "" || $parts[2] eq "");
			
			#Check for ambiguous alleles
			if(defined $ambiguous->{$parts[2]}){
				$parts[2] = $ambiguous->{$parts[2]};
			}
			if(defined $ambiguous->{$parts[1]}){
				$parts[1] = $ambiguous->{$parts[1]};
			}
	
			my @a1_parts = split /[\*:]/, $parts[1];				
			my @a2_parts = split /[\*:]/, $parts[2];
			
			my $gene1 = shift @a1_parts;
			my $gene2 = shift @a2_parts;

			$predictions{$cols[0]}{$gene1}{string}{1}=$parts[1];		
			$predictions{$cols[0]}{$gene2}{string}{2}=$parts[2];		
			if($#a1_parts == 0){ ##2 digit
				$predictions{$cols[0]}{$gene1}{2}{1}{$a1_parts[0]}=1;
				$predictions{$cols[0]}{$gene1}{digits}{2}{1} = 1;
				$predictions{$cols[0]}{$gene1}{exact}{1}{$a1_parts[0]}=1;
			}elsif($#a1_parts == 1){ ##4 digit
				$predictions{$cols[0]}{$gene1}{digits}{2}{1} = 1;
				$predictions{$cols[0]}{$gene1}{digits}{4}{1} = 1;
				$predictions{$cols[0]}{$gene1}{2}{1}{"$a1_parts[0]"}=1;
				$predictions{$cols[0]}{$gene1}{4}{1}{"$a1_parts[0]:$a1_parts[1]"}=1;
				$predictions{$cols[0]}{$gene1}{exact}{1}{"$a1_parts[0]:$a1_parts[1]"}=1;
			}elsif($#a1_parts == 2){ ##6 digit
				$predictions{$cols[0]}{$gene1}{digits}{2}{1} = 1;
				$predictions{$cols[0]}{$gene1}{digits}{4}{1} = 1;
				$predictions{$cols[0]}{$gene1}{digits}{6}{1} = 1;
				$predictions{$cols[0]}{$gene1}{2}{1}{"$a1_parts[0]"}=1;
				$predictions{$cols[0]}{$gene1}{4}{1}{"$a1_parts[0]:$a1_parts[1]"}=1;
				$predictions{$cols[0]}{$gene1}{6}{1}{"$a1_parts[0]:$a1_parts[1]:$a1_parts[2]"}=1;
				$predictions{$cols[0]}{$gene1}{exact}{1}{"$a1_parts[0]:$a1_parts[1]:$a1_parts[2]"}=1;
			}elsif($#a1_parts == 3){ ##8 digit
				$predictions{$cols[0]}{$gene1}{digits}{2}{1} = 1;
				$predictions{$cols[0]}{$gene1}{digits}{4}{1} = 1;
				$predictions{$cols[0]}{$gene1}{digits}{6}{1} = 1;
				$predictions{$cols[0]}{$gene1}{digits}{8}{1} = 1;
				$predictions{$cols[0]}{$gene1}{2}{1}{"$a1_parts[0]"}=1;
				$predictions{$cols[0]}{$gene1}{4}{1}{"$a1_parts[0]:$a1_parts[1]"}=1;
				$predictions{$cols[0]}{$gene1}{6}{1}{"$a1_parts[0]:$a1_parts[1]:$a1_parts[2]"}=1;
				$predictions{$cols[0]}{$gene1}{8}{1}{"$a1_parts[0]:$a1_parts[1]:$a1_parts[2]:$a1_parts[3]"}=1;
				$predictions{$cols[0]}{$gene1}{exact}{1}{"$a1_parts[0]:$a1_parts[1]:$a1_parts[2]:$a1_parts[3]"}=1;
			}
			if($#a2_parts == 0){ ##2 digit
				$predictions{$cols[0]}{$gene2}{2}{2}{$a2_parts[0]}=1;
				$predictions{$cols[0]}{$gene2}{digits}{2}{2} = 1;
				$predictions{$cols[0]}{$gene2}{exact}{2}{$a2_parts[0]}=1;
			}elsif($#a2_parts == 1){ ##4 digit
				$predictions{$cols[0]}{$gene2}{digits}{4}{2} = 1;
				$predictions{$cols[0]}{$gene2}{2}{2}{"$a2_parts[0]"}=1;
				$predictions{$cols[0]}{$gene2}{4}{2}{"$a2_parts[0]:$a2_parts[1]"}=1;
				$predictions{$cols[0]}{$gene2}{exact}{2}{"$a2_parts[0]:$a2_parts[1]"}=1;
			}elsif($#a2_parts == 2){ ##6 digit
				$predictions{$cols[0]}{$gene2}{digits}{6}{2} = 1;
				$predictions{$cols[0]}{$gene2}{2}{2}{"$a2_parts[0]"}=1;
				$predictions{$cols[0]}{$gene2}{4}{2}{"$a2_parts[0]:$a2_parts[1]"}=1;
				$predictions{$cols[0]}{$gene2}{6}{2}{"$a2_parts[0]:$a2_parts[1]:$a2_parts[2]"}=1;
				$predictions{$cols[0]}{$gene2}{exact}{2}{"$a2_parts[0]:$a2_parts[1]:$a2_parts[2]"}=1;
			}elsif($#a2_parts == 3){ ##8 digit
				$predictions{$cols[0]}{$gene2}{digits}{8}{2} = 1;
				$predictions{$cols[0]}{$gene2}{2}{2}{"$a2_parts[0]"}=1;
				$predictions{$cols[0]}{$gene2}{4}{2}{"$a2_parts[0]:$a2_parts[1]"}=1;
				$predictions{$cols[0]}{$gene2}{6}{2}{"$a2_parts[0]:$a2_parts[1]:$a2_parts[2]"}=1;
				$predictions{$cols[0]}{$gene2}{8}{2}{"$a2_parts[0]:$a2_parts[1]:$a2_parts[2]:$a2_parts[3]"}=1;
				$predictions{$cols[0]}{$gene2}{exact}{2}{"$a2_parts[0]:$a2_parts[1]:$a2_parts[2]:$a2_parts[3]"}=1;
			}
							
		}
		close($pfh);
	}
	close($fh);
	return (\%predictions);
}

sub populateTruth{
	my $truth_file = shift;
	my $truth = shift;
	my $ambiguous = shift;
	
	##Open truth file
	my $tfh;
	open($tfh, $truth_file);

	##Header will be used to deduce gene name if necessary
	my $truth_header = <$tfh>;
	chomp($truth_header);

	my @truth_cols = split /[\t ]/, $truth_header;
	my %genes = ();
	while(<$tfh>){
		chomp;
		$_=~s/"//g;
		$_=~s/\r//g;
		my @cols = split /[ \t]/, $_;

		##Iterate through these truth calls.
		for(my $i = 2; $i<=$#cols; $i++){
			my @alleles = split /\//, $cols[$i];
			my $gene = $truth_cols[$i];
			my $aname = "";
			if($cols[$i] =~ m/(.*)\*/){
				$gene = $1;
				$aname=$cols[$i];
			}else{
				$aname="$gene\*$cols[$i]";
			}
			$genes{$gene} = 1;
			my $num = 1;
			$num = 2 if(exists $truth->{$cols[0]}{$gene}{exact}{1});
			$truth->{$cols[0]}{$gene}{string}{$num} = $cols[$i];
			$truth->{$cols[0]}{$gene}{subset}{$num} = 0;
			if(defined $subset{$cols[0]}{$aname}){
				$truth->{$cols[0]}{$gene}{subset}{$num} = 1;
			}
			my $subset = 0;
			##Some truths may have multiple alleles possible for truth. Iterate through these.
			foreach my $allele (@alleles){
				if(exists $ambiguous->{$allele}){
					$allele = $ambiguous->{$allele};
				}elsif(exists $ambiguous->{"$gene*$allele"}){
					$allele = $ambiguous->{"$gene*$allele"};
				
				}
				my @allele_parts = split /[:\*]/, $allele;
				if($allele =~ m/^$gene\*/){
					shift(@allele_parts);
				}
				if($#allele_parts == 0){ ##2 digit
					$truth->{$cols[0]}{$gene}{2}{$num}{$allele_parts[0]}=1;
					$truth->{$cols[0]}{$gene}{digits}{2}{$num} = 1;
					$truth->{$cols[0]}{$gene}{exact}{$num}{$allele_parts[0]}=1;
				}elsif($#allele_parts == 1){ ##4 digit
					$truth->{$cols[0]}{$gene}{digits}{2}{$num} = 1;
					$truth->{$cols[0]}{$gene}{digits}{4}{$num} = 1;
					$truth->{$cols[0]}{$gene}{2}{$num}{"$allele_parts[0]"}=1;
					$truth->{$cols[0]}{$gene}{4}{$num}{"$allele_parts[0]:$allele_parts[1]"}=1;
					$truth->{$cols[0]}{$gene}{exact}{$num}{"$allele_parts[0]:$allele_parts[1]"}=1;
				}elsif($#allele_parts == 2){ ##6 digit
					$truth->{$cols[0]}{$gene}{digits}{2}{$num} = 1;
					$truth->{$cols[0]}{$gene}{digits}{4}{$num} = 1;
					$truth->{$cols[0]}{$gene}{digits}{6}{$num} = 1;
					$truth->{$cols[0]}{$gene}{2}{$num}{"$allele_parts[0]"}=1;
					$truth->{$cols[0]}{$gene}{4}{$num}{"$allele_parts[0]:$allele_parts[1]"}=1;
					$truth->{$cols[0]}{$gene}{6}{$num}{"$allele_parts[0]:$allele_parts[1]:$allele_parts[2]"}=1;
					$truth->{$cols[0]}{$gene}{exact}{$num}{"$allele_parts[0]:$allele_parts[1]:$allele_parts[2]"}=1;
				}elsif($#allele_parts == 3){ ##8 digit
					$truth->{$cols[0]}{$gene}{digits}{2}{$num} = 1;
					$truth->{$cols[0]}{$gene}{digits}{4}{$num} = 1;
					$truth->{$cols[0]}{$gene}{digits}{6}{$num} = 1;
					$truth->{$cols[0]}{$gene}{digits}{8}{$num} = 1;
					$truth->{$cols[0]}{$gene}{2}{$num}{"$allele_parts[0]"}=1;
					$truth->{$cols[0]}{$gene}{4}{$num}{"$allele_parts[0]:$allele_parts[1]"}=1;
					$truth->{$cols[0]}{$gene}{6}{$num}{"$allele_parts[0]:$allele_parts[1]:$allele_parts[2]"}=1;
					$truth->{$cols[0]}{$gene}{8}{$num}{"$allele_parts[0]:$allele_parts[1]:$allele_parts[2]:$allele_parts[3]"}=1;
					$truth->{$cols[0]}{$gene}{exact}{$num}{"$allele_parts[0]:$allele_parts[1]:$allele_parts[2]:$allele_parts[3]"}=1;
				}
			}
		}
	}
	close($tfh);
}

sub checkConcordance{
	my $answers = shift;
	my $truth = shift;
	my $stats = shift;
	my %wrong = ();	
	foreach my $sample (keys %{$answers}){
		foreach my $gene (keys %{$answers->{$sample}}){
			if(exists ($truth->{$sample}) && exists ($truth->{$sample}{$gene})){
				#Calculate how many truth alleles (0,1 or both) the sample has containing each precision
				my $tnum_in2_1 = $truth->{$sample}{$gene}{digits}{2}{1} || 0; 
				my $tnum_in4_1 = $truth->{$sample}{$gene}{digits}{4}{1} || 0; 
				my $tnum_in6_1 = $truth->{$sample}{$gene}{digits}{6}{1} || 0; 
				my $tnum_in8_1 = $truth->{$sample}{$gene}{digits}{8}{1} || 0; 
				my $tnum_in2_2 = $truth->{$sample}{$gene}{digits}{2}{2} || 0; 
				my $tnum_in4_2 = $truth->{$sample}{$gene}{digits}{4}{2} || 0; 
				my $tnum_in6_2 = $truth->{$sample}{$gene}{digits}{6}{2} || 0; 
				my $tnum_in8_2 = $truth->{$sample}{$gene}{digits}{8}{2} || 0; 
				my $tnum_in2 = $tnum_in2_1 + $tnum_in2_2;
				my $tnum_in4 = $tnum_in4_1 + $tnum_in4_2;
				my $tnum_in6 = $tnum_in6_1 + $tnum_in6_2;
				my $tnum_in8 = $tnum_in8_1 + $tnum_in8_2;

				#Calculate how many alleles (0,1, or both) of each precision are predicted by the sample
				my $num_in2_1 = $predictions->{$sample}{$gene}{digits}{2}{1} || 0; 
				my $num_in4_1 = $predictions->{$sample}{$gene}{digits}{4}{1} || 0; 
				my $num_in6_1 = $predictions->{$sample}{$gene}{digits}{6}{1} || 0; 
				my $num_in8_1 = $predictions->{$sample}{$gene}{digits}{8}{1} || 0; 
				my $num_in2_2 = $predictions->{$sample}{$gene}{digits}{2}{2} || 0; 
				my $num_in4_2 = $predictions->{$sample}{$gene}{digits}{4}{2} || 0; 
				my $num_in6_2 = $predictions->{$sample}{$gene}{digits}{6}{2} || 0; 
				my $num_in8_2 = $predictions->{$sample}{$gene}{digits}{8}{2} || 0;
				if($answers->{$sample}{$gene}{string}{1} eq "NC"){
					$num_in2_1 =1;	
					$num_in4_1 =1;	
					$num_in6_1 =1;	
					$num_in8_1 =1;	
				}	 
				if($answers->{$sample}{$gene}{string}{2} eq "NC"){
					$num_in2_2 =1;	
					$num_in4_2 =1;	
					$num_in6_2 =1;	
					$num_in8_2 =1;	
				}	 
				my $num_in2 = $num_in2_1 + $num_in2_2;
				my $num_in4 = $num_in4_1 + $num_in4_2;
				my $num_in6 = $num_in6_1 + $num_in6_2;
				my $num_in8 = $num_in8_1 + $num_in8_2;
				
				if($num_in2 > 0 && $tnum_in2 > 0){
					my ($matches,$match1,$match2) = checkTruth($gene, $truth->{$sample}{$gene}{2},$answers->{$sample}{$gene}{2}{1},$answers->{$sample}{$gene}{2}{2});
					my $possible = $tnum_in2;
					$possible = $tnum_in2 if($num_in2 <$tnum_in2);
					if($matches >= 2 && $possible == 2){
						$stats->{$gene}{2}{correct}+=2;
						$stats->{$gene}{2}{total}+=2;
					}elsif($matches ==1 && $possible >= 1){
						$stats->{$gene}{2}{correct}+=1;
						$stats->{$gene}{2}{total}+=$possible;
						$wrong{2}{$gene}{$sample}="2 digit: Truth $truth->{$sample}{$gene}{string}{1}|$truth->{$sample}{$gene}{string}{2} does not match $answers->{$sample}{$gene}{string}{1}|$answers->{$sample}{$gene}{string}{2}\t1" if ($possible > 1);
					}elsif ($possible > 0){
						$stats->{$gene}{2}{correct}+=0;
						$stats->{$gene}{2}{total}+=$possible;
						$wrong{2}{$gene}{$sample}="2 digit: Truth $truth->{$sample}{$gene}{string}{1}|$truth->{$sample}{$gene}{string}{2} does not match $answers->{$sample}{$gene}{string}{1}|$answers->{$sample}{$gene}{string}{2}\t0";
					}else{
						$stats->{$gene}{2}{correct}+=0;
						$stats->{$gene}{2}{total}+=0;
					}
					if($opts{allele_subset} ne ""){
						if($truth->{$sample}{$gene}{subset}{1} == 1){
							$stats->{$gene}{subset}{2}{correct}+=$match1;
							$stats->{$gene}{subset}{2}{total}+=1;
						}
						if($truth->{$sample}{$gene}{subset}{2} == 1){
							$stats->{$gene}{subset}{2}{correct}+=$match2;
							$stats->{$gene}{subset}{2}{total}+=1;
						}
					}
				}	
				if($num_in4 > 0 && $tnum_in4 > 0){
					my ($matches,$match1,$match2) = checkTruth($gene, $truth->{$sample}{$gene}{4},$answers->{$sample}{$gene}{4}{1},$answers->{$sample}{$gene}{4}{2});
					my $possible = $tnum_in4;
					$possible = $tnum_in4 if($num_in4 <$tnum_in4);
					if($matches >= 2 && $possible == 2){
						$stats->{$gene}{4}{correct}+=2;
						$stats->{$gene}{4}{total}+=2;
					}elsif($matches ==2 && $possible == 1){
						$stats->{$gene}{4}{correct}+=1;
						$stats->{$gene}{4}{total}+=$possible;
					}elsif($matches ==1 && $possible >= 1){
						$stats->{$gene}{4}{correct}+=1;
						$stats->{$gene}{4}{total}+=$possible;
						$wrong{4}{$gene}{$sample}="4 digit: Truth $truth->{$sample}{$gene}{string}{1}|$truth->{$sample}{$gene}{string}{2} does not match $answers->{$sample}{$gene}{string}{1}|$answers->{$sample}{$gene}{string}{2}\t1" if ($possible > 1);
					}elsif ($possible > 0){
						$stats->{$gene}{4}{correct}+=0;
						$stats->{$gene}{4}{total}+=$possible;
						$wrong{4}{$gene}{$sample}="4 digit: Truth $truth->{$sample}{$gene}{string}{1}|$truth->{$sample}{$gene}{string}{2} does not match $answers->{$sample}{$gene}{string}{1}|$answers->{$sample}{$gene}{string}{2}\t0";
					}else{
						$stats->{$gene}{4}{correct}+=0;
						$stats->{$gene}{4}{total}+=0;
					}
					if($opts{allele_subset} ne ""){
						if($truth->{$sample}{$gene}{subset}{1} == 1){
							$stats->{$gene}{subset}{4}{correct}+=$match1;
							$stats->{$gene}{subset}{4}{total}+=1;
						}
						if($truth->{$sample}{$gene}{subset}{2} == 1){
							$stats->{$gene}{subset}{4}{correct}+=$match2;
							$stats->{$gene}{subset}{4}{total}+=1;
						}
						my $subset = "";
						my $match = "does not match";
						my $sm1 = "0";
						my $sm2 = "0";
						if($truth->{$sample}{$gene}{subset}{1} == 1){
							$stats->{$gene}{subset}{exact}{correct}+=$match1;
							$stats->{$gene}{subset}{exact}{total}+=1;
							$subset = "$truth->{$sample}{$gene}{string}{1}/";
							$match = "matches" if ($match1 == 1);
							$sm1 = $match1;
						}
						if($truth->{$sample}{$gene}{subset}{2} == 1){
							$stats->{$gene}{subset}{exact}{correct}+=$match2;
							$stats->{$gene}{subset}{exact}{total}+=1;
							$subset .= $truth->{$sample}{$gene}{string}{2};
							$match = "matches" if ($match2 == 1);
							$sm1 = $match2;
						}

						if($subset ne ""){
							$subset =~ s/\///;
							$wrong{subset}{$gene}{$sample}="Subset: Subset truth $subset $match $answers->{$sample}{$gene}{string}{1}|$answers->{$sample}{$gene}{string}{2}\t$sm1$sm2"; 
						}
					}	
				}	
				if($num_in6 > 0 && $tnum_in6 > 0){
					my ($matches,$match1,$match2) = checkTruth($gene, $truth->{$sample}{$gene}{6},$answers->{$sample}{$gene}{6}{1},$answers->{$sample}{$gene}{6}{2});
					my $possible = $tnum_in6;
					$possible = $num_in6 if($num_in6 <$tnum_in6);
					if($matches >= 2 && $possible == 2){
						$stats->{$gene}{6}{correct}+=2;
						$stats->{$gene}{6}{total}+=2;
					}elsif($matches ==1 && $possible >= 1){
						$stats->{$gene}{6}{correct}+=1;
						$stats->{$gene}{6}{total}+=$possible;
						$wrong{6}{$gene}{$sample}="6 digit: Truth $truth->{$sample}{$gene}{string}{1}|$truth->{$sample}{$gene}{string}{2} does not match $answers->{$sample}{$gene}{string}{1}|$answers->{$sample}{$gene}{string}{2}\t1" if ($possible > 1);
					}elsif ($possible > 0){
						$stats->{$gene}{6}{correct}+=0;
						$stats->{$gene}{6}{total}+=$possible;
						$wrong{6}{$gene}{$sample}="6 digit: Truth $truth->{$sample}{$gene}{string}{1}|$truth->{$sample}{$gene}{string}{2} does not match $answers->{$sample}{$gene}{string}{1}|$answers->{$sample}{$gene}{string}{2}\t0";
					}else{
						$stats->{$gene}{6}{correct}+=0;
						$stats->{$gene}{6}{total}+=0;
					}	
					if($opts{allele_subset} ne ""){
						if($truth->{$sample}{$gene}{subset}{1} == 1){
							$stats->{$gene}{subset}{6}{correct}+=$match1;
							$stats->{$gene}{subset}{6}{total}+=1;
						}
						if($truth->{$sample}{$gene}{subset}{2} == 1){
							$stats->{$gene}{subset}{6}{correct}+=$match2;
							$stats->{$gene}{subset}{6}{total}+=1;
						}
					}
				}	
				if($num_in8 > 0 && $tnum_in8 > 0){
					my ($matches,$match1,$match2) = checkTruth($gene, $truth->{$sample}{$gene}{8},$answers->{$sample}{$gene}{8}{1},$answers->{$sample}{$gene}{8}{2});
					my $possible = $tnum_in8;
					$possible = $num_in8 if($num_in8 <$tnum_in8);
					if($matches >= 2 && $possible == 2){
						$stats->{$gene}{8}{correct}+=2;
						$stats->{$gene}{8}{total}+=2;
					}elsif($matches ==1 && $possible >= 1){
						$stats->{$gene}{8}{correct}+=1;
						$stats->{$gene}{8}{total}+=1;
						$wrong{8}{$gene}{$sample}="8 digit: Truth $truth->{$sample}{$gene}{string}{1}|$truth->{$sample}{$gene}{string}{2} does not match $answers->{$sample}{$gene}{string}{1}|$answers->{$sample}{$gene}{string}{2}\t1" if ($possible > 1);
					}elsif ($possible > 0){
						$stats->{$gene}{8}{correct}+=0;
						$stats->{$gene}{8}{total}+=$possible;
						$wrong{8}{$gene}{$sample}="8 digit: Truth $truth->{$sample}{$gene}{string}{1}|$truth->{$sample}{$gene}{string}{2} does not match $answers->{$sample}{$gene}{string}{1}|$answers->{$sample}{$gene}{string}{2}\t0";
					}else{
						$stats->{$gene}{8}{correct}+=0;
						$stats->{$gene}{8}{total}+=0;
					}
					if($opts{allele_subset} ne ""){
						if($truth->{$sample}{$gene}{subset}{1} == 1){
							$stats->{$gene}{subset}{8}{correct}+=$match1;
							$stats->{$gene}{subset}{8}{total}+=1;
						}
						if($truth->{$sample}{$gene}{subset}{2} == 1){
							$stats->{$gene}{subset}{8}{correct}+=$match2;
							$stats->{$gene}{subset}{8}{total}+=1;
						}
					}
				}
				my ($matches,$match1,$match2) = checkTruth($gene, $truth->{$sample}{$gene}{exact},$answers->{$sample}{$gene}{exact}{1},$answers->{$sample}{$gene}{exact}{2});
				$stats->{$gene}{exact}{correct}+=$matches;
				$stats->{$gene}{exact}{total}+=2;
				
				if($matches <2){
					$wrong{exact}{$gene}{$sample}="Exact: Truth $truth->{$sample}{$gene}{string}{1}|$truth->{$sample}{$gene}{string}{2} does not match $answers->{$sample}{$gene}{string}{1}|$answers->{$sample}{$gene}{string}{2}";
				}
			} 
		}
	}		
	return(\%wrong);
}

sub checkTruth{
	my $gene = shift;
	my $truth = shift;
	my $a1_ref = shift;
	my $a2_ref = shift;
	my @a1 = keys %{$a1_ref};
	my @a2 = keys %{$a2_ref};
	my $match1 = 0;
	my $match2 = 0;
	#Initialize matches to 0	
	my $matches = 0;
	#Iterate over all combinations of allele1 and allele2 answers to find a match
	for (my $i = 0; $i<=$#a1; $i++){
		for(my $j=0; $j<=$#a2; $j++){
			if(exists $truth->{1}{$a1[$i]} && exists $truth->{2}{$a2[$j]}){
				$matches = 2;
				$match1 =1; 
				$match2 =1; 
			}elsif(exists $truth->{1}{$a2[$j]} && exists $truth->{2}{$a1[$i]}){
				$matches = 2;
				$match1 =1; 
				$match2 =1; 
			}elsif(exists $truth->{1}{$a1[$i]}){
				if ($matches < 1){
					$match1 = 1;
					$matches = 1;	
				}
			}elsif (exists $truth->{2}{$a1[$i]}){
				if ($matches < 1){
					$match2= 1;
					$matches = 1;	
				}
			}elsif(exists $truth->{1}{$a2[$j]}){
				if ($matches < 1){
					$match1 = 1;
					$matches = 1;	
				}
			}elsif (exists $truth->{2}{$a2[$j]}){
				if ($matches < 1){
					$match2 = 1;
					$matches = 1;	
				}
			}
		}
	}
	if($#a1 == -1){
		for (my $i = 0; $i<=$#a2; $i++){
			if(exists $truth->{1}{$a2[$i]} || exists $truth->{2}{$a2[$i]}){
				$matches = 1;
				$match2 = 1;
			}
		}
	}
	if($#a2 == -1){
		for (my $i = 0; $i<=$#a1; $i++){
			if(exists $truth->{1}{$a1[$i]} || exists $truth->{2}{$a1[$i]}){
				$matches = 1;
				$match1 = 1;
			}
		}
	}
	return($matches,$match1,$match2);	
}

sub summarizeStats{
	my $stats = shift;
	my @order = (2,4,6,8,"exact");
	print "Gene";
	foreach my $p (@order){
		print "\t${p}_correct\t${p}_possible\t${p}_percent";
	}
	print "\n";
	GENE:foreach my $gene (keys %{$stats}){
		print "$gene";
		foreach my $p (@order){
			my $count = $stats->{$gene}{$p}{correct} || 0;
			my $total = $stats->{$gene}{$p}{total} || 0;
			my $perc = "0.00";
			$perc = sprintf("%.2f", $count/$total*100) if($total != 0);
			print "\t$count\t$total\t$perc";
		}
		print "\n";
	}
}

sub summarizeSubsetStats{
	my $stats = shift;
	my @order = (2,4,6,8,"exact");
	print "Gene";
	foreach my $p (@order){
		print "\t${p}_correct\t${p}_possible\t${p}_percent";
	}
	print "\n";
	GENE:foreach my $gene (keys %{$stats}){
		print "$gene";
		foreach my $p (@order){
			my $count = $stats->{$gene}{subset}{$p}{correct} || 0;
			my $total = $stats->{$gene}{subset}{$p}{total} || 0;
			my $perc = "0.00";
			$perc = sprintf("%.2f", $count/$total*100) if($total != 0);
			print "\t$count\t$total\t$perc";
		}
		print "\n";
	}
}

sub printWrongAnswers{
	my $wrong = shift;
	my @order = (2,4,6,8,"exact", "subset");
	foreach my $p (@order){
		if(exists $wrong->{$p}){
			my %genes= %{$wrong->{$p}};
			foreach my $gene (keys %genes){
				my %hash = %{$wrong->{$p}{$gene}};
				foreach my $key (keys %hash){
					print $log "$key\t$hash{$key}\n";
				}
			}
		}
	}	
}
