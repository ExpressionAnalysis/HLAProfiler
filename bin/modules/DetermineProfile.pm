#!/usr/bin/env perl
package DetermineProfile;
use strict;
use warnings;
use Getopt::Long;
use Storable;
use Module::Load;

(my $SCRIPT_NAME = $0) =~ s/.*\///;
my $version = "1.0";
my $creation_date = "1 Nov 2016";
my $last_updated = "11 Jan 2017";

my $usage = "\n$SCRIPT_NAME v$version\n" .
	    "\nDESCRIPTION:\n" .
	    "A perl module created to identify gene profiles using sequence counts of reads simulated (assuming stranded protocol) from each allele. Sequence counts files are tab-delimeted file containing a sequence in column 1 and the number of times that sequence was observed in column 2. Determine profile iterates over the counts files in a specified directory and calculates the profiles associated of each allele of the specified gene. While designed to be called from within HLAProfiler, this script can be run from the commandline as well.\n" .
	    "\nUSAGE:\n" .
	    "perl $SCRIPT_NAME <options>\n" .
	    "\nRequired options:\n" .
	    "-gene|g\t\t\tThe gene from which the alleles originate. (requred)\n" .
	    "-simulation_dir|sim\tThe directory containing the sequence counts files for the simulated reads. (required)\n" .
	    "-prefix|p\t\tThe portion of the name of the sequence counts files preceding the period just before the allele name, i.e simulatedReads for simulatedReads.HLA00001.A_1.fq (required)\n" .
	    "\nGeneral options:\n" .
	    "-print_hash|ph\t\tPrints the profile hash to a storable file. (default:ON)\n" .
	    "-print_matrix|pm\tPrints the profile hash as a tab-delimeted matrix. (default:OFF)\n" .
	    "-output_dir|o\t\tThe directory where the profile will be stored. (default:'.')\n" .
	    "-scripts_dir|sd\t\tThe parent directory of module folder containing the HLAProfiler modules (default:'..')\n" . 
	    "-help|h\t\t\tDisplays this message\n" .
            "\nAUTHORS:\n" .
	    "Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	    "Chad Brown:chad.brown\@q2labsolutions.com\n" .
	    "\nCREATED:\n$creation_date\n" .
	    "\nLAST UPDATED:\n$last_updated\n" .
	    "\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
	    "\n";

sub runCommandline{
	my %opts = (help=>0, gene=>"", simulation_dir=>"", prefix=>"", output_dir=>".",scripts_dir=>"..", print_hash=>1, print_matrix=>0);
	GetOptions(\%opts, qw(help|h gene|g=s simulation_dir|sim=s prefix|p=s output_dir|o=s scripts_dir|sd=s print_hash|ph print_matrix|pm));

	if($opts{help}){
		print "$usage";
		exit;
	}elsif($opts{gene} eq ""){
		print STDERR "Please specify a gene. $usage\n";
		exit 1;	
	}elsif($opts{prefix} eq ""){
		print STDERR "Please specify a file prefix. $usage\n";	
		exit 1;
	}elsif($opts{simulation_dir} eq ""){
		print STDERR "Please specify a simulation directory. $usage\n";	
		exit 1;
	}else{
		createProfile($opts{gene}, $opts{simulation_dir}, $opts{prefix}, $opts{output_dir}, $opts{print_hash}, $opts{print_matrix}, $opts{scripts_dir});
	}
}

sub createProfile{
	my $gene = shift;
	if($gene eq "DetermineProfile"){
		$gene = shift;
	}
	my $simulation_directory = shift;
	my $prefix = shift;
	my $output_directory = shift;
	my $print_hash = shift;
	my $print_matrix = shift;
	my $scripts_dir = shift;
	load "$scripts_dir/modules/SequenceFunctions.pm";
	my $sig = determineProfile($gene, $simulation_directory, $prefix);
	if($print_hash){
		printProfileToHashFile($sig, "$output_directory/$gene.profile.ph");
	}
	if($print_matrix){
		printProfileToMatrixFile($sig, "$output_directory/$gene.profile.txt");
	}
}

sub printProfileToHashFile{
	my $sig_ref = shift;
	my $file = shift;
	store $sig_ref, $file;
}

sub printProfileToMatrixFile{
	my $sig_ref = shift;
	my $file = shift;
	my %s = %{$sig_ref->{sig}};
	my %t = %{$sig_ref->{totals}};
	my %n = %{$sig_ref->{names}};
	my %k = %{$sig_ref->{kmers}};
	
	open(OUT, ">$file");	
	my @kmers =sort keys %k;
	print OUT "Sample";
	foreach my $kmer (@kmers){
		print OUT "\t$kmer";
	}
	print OUT "\tTotal\n";
	foreach my $allele (sort keys %s){
		print OUT "$allele";
		foreach my $kmer (@kmers){
			my $count = $s{$allele}{$k{$kmer}} || 0;
			print OUT "\t$count";
		}
		print OUT "\t$t{$allele}\n";
	}
	close OUT;
}

sub determineProfile{
	my $gene = shift;
	if($gene eq "DetermineProfile"){
		$gene = shift;
	}
	my $simulation_directory = shift;
	my $prefix = shift;

	my %kmers = ();
	my %kmer_keys = ();
	my %allele_counts = ();
	my %names = ();

	$simulation_directory =~ s/\/$//;

	my $kmer_cnt = 0;
	#Read in counts files and store counts in profile hash
	opendir(my $dh, "$simulation_directory/");
	my @files = sort readdir $dh;
	foreach my $file (@files){
		if($file=~m/^$prefix\.(.*)\.${gene}_([12]).uniq.cnts/){
			my $allele = $1;
			my $read = $2;
			open(IN, "$simulation_directory/$file");
			while(<IN>){
				chomp;
				my @parts = split /\t/, $_;
				my $key = "";
				$parts[0] = SequenceFunctions->revcomp($parts[0]) if ($read == 1);
				if(defined $kmer_keys{$parts[0]}){
					$key = $kmer_keys{$parts[0]};
				}else{
					#The unique identifier of the profile is the kmer_cnt
					$key = $kmer_cnt;
					#Maps kmer to the key
					$kmer_keys{$parts[0]}=$key;
					#Maps key to kmer
					$names{$key}=$parts[0];
					$kmer_cnt++;
				}
				$kmers{$allele}{$key}+=$parts[1];
				$allele_counts{$allele}+=$parts[1];
			}
			close(IN);
		}		
	}
	close($dh);
	my %profile = ();
	##Hash of alleles as keys and another hash as a value.The value has has the kmer unique identifiers as keys and counts as values
	$profile{sig}=\%kmers;
	##Key: allele Value: allele_counts
	$profile{totals}=\%allele_counts;
	##Key:kmer Value:unique identifier
	$profile{names}=\%names;
	##Value:unique identifier Key:kmer
	$profile{kmers}=\%kmer_keys;
	return(\%profile);	
}

__PACKAGE__->runCommandline() unless caller;
