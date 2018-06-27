#!/usr/bin/env perl
package Tests;
use strict;
#use warnings;
use Test::More;
use Getopt::Long;
use Module::Load;
use File::Slurp;
use File::Compare;
use Storable;
use FindBin;
use Test::Trap qw/ :output(systemsafe) /;
use lib "$FindBin::Bin/../bin/modules";

my $creation_date = "1 Feb 2017";
my $last_updated = "14 Sep 2017";

my $usage = "$0\n" .
	     	"DESCRIPTION\n" .
	     	"\nTest suite for HLAProfiler.pl. This test suite will test all or parts of HLAProfiler.pl and the associated modules.\n" .
	     	"\nUSAGE:\n" .
	     	"perl $0 <options>\n" .
	     	"\nRequired options\n" .
		"-test|t\t\tDenotes the module to test\n" .
		"\t\tAvailable_tests:\n" .
		"\t\tall\n" .
		"\t\tSequenceFunctions\n" .
		"\t\tMergeDuplicates\n" .
		"\t\tHLATaxonomy\n" .
		"\t\tHLADistractome\n" .
		"\t\tTaxonomyDivisions\n" .
		"\t\tRunKraken\n" .
		"\t\tSimulateReads\n" .
		"\t\tReadCounter\n" .
		"\t\tDetermineProfile\n" .
		"\t\tPairPicker\n" .
		"\t\tAlleleRefiner\n" .
		#"\t\tHLAPredict\n" .
		"\nGeneral options:\n" .
		"-kraken_path|kp\tbase directory of kraken installation. (default:base directory of path returned by `which kraken`)\n" .
		"-directory|d\tlocation of test files. (default:;'.')\n" .
		"-output_directory|od\tlocation of temporary output files. (default:;'.')\n" .
		"-help|h\t\tprints this help prompt\n" .
		"\nAUTHORS:\n" .
		"Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .  
		"\nCREATED:\n$creation_date\n" .
		"\nLAST UPDATED:\n$last_updated\n" .
		"\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .   
		"\n"; 

my %opts = (test=>"all", kraken_path=>".", directory=>".", output_directory=>".");

my $reference_fasta;
my $reference_merged_fasta;
my $cwd_file;
my $transcripts_fasta;
my $transcripts_gtf;
my $exclude_bed;
my $names; 
my $nodes; 
my $divisions;
my $distractome;
my $kraken_fasta;
my $outdir;
my $indir;
my $kraken_path;
my %loaded_modules = ();


sub runCommandline{
	GetOptions(\%opts, qw(help|h test|t:s kraken_path|kp:s directory|d:s));
	if($opts{help}){
		print "$usage";
		exit;
	}else{
		runTests($opts{test}, $opts{kraken_path}, $opts{test_directory}, $opts{output_directory});
	}
}

sub runTests{
	my $test = shift;
	$test =shift if ($test eq "Tests");
	$kraken_path = shift;
	my $test_directory = shift || ".";
	my $output_directory = shift || ".";
	createDirectories($output_directory);
	setFilePaths($test_directory);
	if($test eq "all"){
		
		#Run Tests
		testSequenceFunctions();
		testMergeDuplicates();
		testHLATaxonomy();
		testHLADistractome();
		testTaxonomyDivisions();
		testRunKraken();
		testSimulateReads();
		testReadCounter();
		testDetermineProfile();
		testPairPicker();
		testAlleleRefiner();
		#testHLAPredict();
		#testHLAProfiler();
		done_testing(11);
	}elsif($test eq "SequenceFunctions"){
		testSequenceFunctions();
		done_testing(1);
	}elsif($test eq "MergeDuplicates"){
		testMergeDuplicates();
		done_testing(1);
	}elsif($test eq "HLATaxonomy"){
		testHLATaxonomy();
		done_testing(1);
	}elsif($test eq "HLADistractome"){
		testHLADistractome();
		done_testing(1);
	}elsif($test eq "TaxonomyDivisions"){
		testTaxonomyDivisions();
		done_testing(1);
	}elsif($test eq "RunKraken"){
		testRunKraken();
		done_testing(1);
	}elsif($test eq "SimulateReads"){
		testSimulateReads();
		done_testing(1);
	}elsif($test eq "ReadCounter"){
		testReadCounter();
		done_testing(1);
	}elsif($test eq "DetermineProfile"){
		testDetermineProfile();
		done_testing(1);
	}elsif($test eq "PairPicker"){
		testPairPicker();
		done_testing(1);
	}elsif($test eq "AlleleRefiner"){
		testAlleleRefiner();
		done_testing(1);
	}elsif($test eq "HLAPredict"){
		testHLAPredict();
		done_testing(1);
	}
}


sub createDirectories{
	my $dir = shift;
	$outdir = "$dir/test_output";	
	if(-e "$outdir" && -d "$outdir"){	
	}else{
		print "Directory $outdir does not exist. Creating directory.\n";
		`mkdir -p "$outdir"`;
	}
}

sub setFilePaths{
	my $dir = shift;
	$reference_fasta = "$dir/inputs/hla.ref.fa";
	$reference_merged_fasta = "$dir/outputs/hla.ref.merged.fa";
	$cwd_file = "$dir/inputs/hla.cwd.txt";
	$transcripts_fasta = "$dir/inputs/transcripts.fa.gz";
	$transcripts_gtf = "$dir/inputs/transcripts.gtf.gz";
	$exclude_bed = "$dir/inputs/exclude.bed";
	$names = "$dir/outputs/database/taxonomy/names.dmp"; 
	$nodes = "$dir/outputs/database/taxonomy/nodes.dmp"; 
	$divisions = "$dir/outputs/database/taxonomy/divisions.txt";
	$distractome = "$dir/outputs/database/data/reference/distractome.fa";
	$kraken_fasta = "$dir/outputs/database/data/reference/hla.ref.forKraken.fa";
	$indir = $dir;
}

sub compareFiles{
	my $expected_file = shift;
	my $observed_file = shift;
	my $test_name = shift;
	my $observed = read_file($observed_file);
	my $expected = read_file($expected_file);
	is($observed, $expected, $test_name);	
}

sub compareFileToString{
	my $expected = shift;
	my $observed_file = shift;
	my $test_name = shift;
	my $observed = read_file($observed_file);
	is($expected, $observed, $test_name);	
}

sub testMergeDuplicates{
	subtest 'MergeDuplicates.pm' => sub{
		plan tests=>11;
		## Cleaning merge outputs from test_output";
		if(-e "$outdir/hla.ref.function.merged.fa"){
			`rm $outdir/hla.ref.function.merged.fa`;
		}
		if(-e "$outdir/hla.ref.function.merged.allele_map.txt"){
			`rm $outdir/hla.ref.function.merged.allele_map.txt`;
		}
		if(-e "$outdir/hla.ref.commandline.merged.fa"){
			`rm $outdir/hla.ref.commandline.merged.fa`;
		}
		if(-e "$outdir/hla.ref.commandline.merged.allele_map.txt"){
			`rm $outdir/hla.ref.commandline.merged.allele_map.txt`;
		}
		if(-e "$outdir/MergeDuplicates.noreference.err"){
			`rm $outdir/MergeDuplicates.noreference.err`;
		}	
		if(-e "$outdir/MergeDuplicates.nocwd.err"){
			`rm $outdir/MergeDuplicates.nocwd.err`;
		}
		
		#Test 1 - Module load
		use_ok("MergeDuplicates");
		$loaded_modules{"MergeDuplicates"}=1;
		### Test mergeDuplicates function
		my @r = trap {MergeDuplicates::mergeDuplicates("test", "$outdir/hla.ref.function",$cwd_file)};
		#Test 2 - Invalid reference exit
		is($trap->exit, 1, "Function mergeDuplicates: No Reference Fasta Exit");
		#Test 3 - Invalid reference STDERR
		is($trap->stderr, "Reference fasta test does not exist.\n", "Function mergeDuplicates: No Reference Fasta error message"); 

		 
		@r = trap {MergeDuplicates::mergeDuplicates($reference_fasta, "$outdir/hla.ref.function","cwd")};
		#Test 4 - Invalid cwd allele exit
		is($trap->exit, 1, "Function mergeDuplicates: No CWD File Exit");
		#Test 5 - Invalid cwd allele STDERR
		is($trap->stderr, "CWD allele file cwd does not exist.\n", "Function mergeDuplicates: No CWD File error message"); 
		
		MergeDuplicates::mergeDuplicates($reference_fasta, "$outdir/hla.ref.function",$cwd_file);
		#Test 6 - Correct reference output
		compareFiles("$indir/outputs/hla.ref.merged.fa", "$outdir/hla.ref.function.merged.fa", "Function mergeDuplicates: Successful fasta output");	
		#Test 7 - Correct allele_map output
		compareFiles("$indir/outputs/hla.ref.merged.allele_map.txt", "$outdir/hla.ref.function.merged.allele_map.txt", "Function mergeDuplicates: Successful allele_map output");	
		
		### Test RunCommandline function
		`perl $FindBin::Bin/../bin/modules/MergeDuplicates.pm -c $cwd_file -o $outdir/hla.ref.commandline 2>&1`;

		#Test 8 - Invalid reference supplied exit
		is($?, 256, "Commandline mergeDuplicates: No Reference Fasta Exit");
		
		my $out = `perl $FindBin::Bin/../bin/modules/MergeDuplicates.pm -r $reference_fasta -o $outdir/hla.ref.commandline 2>&1`;
		#Test 9 - Invalid cwd allele file supplied exit
		is($?, 256, "Commandline mergeDuplicates: No CWD File Exit");

		$out = `perl $FindBin::Bin/../bin/modules/MergeDuplicates.pm -r $reference_fasta -o "$outdir/hla.ref.commandline" -c $cwd_file 2>&1`;
		#Test 10 - Correct reference output	
		compareFiles("$indir/outputs/hla.ref.merged.fa", "$outdir/hla.ref.commandline.merged.fa", "Commandline MergeDuplicates fasta output");	
		#Test 11 - Correct allele_map output	
		compareFiles("$indir/outputs/hla.ref.merged.allele_map.txt", "$outdir/hla.ref.commandline.merged.allele_map.txt", "Commandline MergeDuplicates allele_map output");
	};		
}

sub testHLATaxonomy{
	subtest 'HLATaxonomy.pm' => sub{
		plan tests=>16;

		## Creating test database directories";
		`mkdir -p $outdir/database_function`;
		`mkdir -p $outdir/database_commandline`;

		## Cleaning taxonomy outputs from test_output";
		if(-e "$outdir/database_function/taxonomy/nodes.dmp"){
			`rm $outdir/database_function/taxonomy/nodes.dmp`;
		}	
		if(-e "$outdir/database_function/taxonomy/names.dmp"){
			`rm $outdir/database_function/taxonomy/names.dmp`;
		}
		if(-e "$outdir/database_commandline/taxonomy/nodes.dmp"){
			`rm $outdir/database_commandline/taxonomy/nodes.dmp`;
		}	
		if(-e "$outdir/database_commandline/taxonomy/names.dmp"){
		`rm $outdir/database_commandline/taxonomy/names.dmp`;
		}
	
		#Test 1 - Module load
		use_ok("HLATaxonomy");
		$loaded_modules{"HLATaxonomy"}=1;
	
		my @r = trap {HLATaxonomy::createHLATaxonomy("","$outdir/database_function")};
		#Test 2 - Invalid cwd allele exit
		is($trap->exit, 1, "Function createHLATaxonomy: Invalid reference exit");
		#Test 3 - Invalid cwd allele STDERR
		is($trap->stderr, "Please enter a fastq file that exists and is not empty. Exiting\n", "Function createHLATaxonomy: Invalid reference error message"); 
	
		HLATaxonomy::createHLATaxonomy($reference_merged_fasta,"$outdir/database_function");
		#Test 4 - Data directory
		is(-d "$outdir/database_function/data",1, "Function createHLATaxonomy: data directory created");
		#Test 5 - Data directory
		is(-d "$outdir/database_function/data/reference",1, "Function createHLATaxonomy: data/reference directory created");
		#Test 6 - Data directory
		is(-d "$outdir/database_function/taxonomy",1, "Function createHLATaxonomy: taxonomy directory created");
		#Test 7 - Successful nodes.dmp
		compareFiles("$nodes", "$outdir/database_function/taxonomy/nodes.dmp", "Function createHLATaxonomy: Successful nodes.dmp");
		#Test 8 - Successful names.dmp
		compareFiles("$names", "$outdir/database_function/taxonomy/names.dmp", "Function createHLATaxonomy: Successful names.dmp");
		#Test 9 - Successful kraken fa
		compareFiles("$indir/outputs/database/data/reference/hla.ref.forKraken.fa", "$outdir/database_function/data/reference/hla.ref.forKraken.fa", "Function createHLATaxonomy: Successful kraken reference fasta");

		my $out = `perl $FindBin::Bin/../bin/modules/HLATaxonomy.pm -o "$outdir/database_commandline" 2>&1` ;
		#Test 10 - Invalid cwd allele exit
		is($?, 256, "Commandline HLATaxonomy: Invalid reference exit");
	
		`perl $FindBin::Bin/../bin/modules/HLATaxonomy.pm -r $reference_merged_fasta -o "$outdir/database_commandline"` ;
		#Test 11 - Data directory
		is(-d "$outdir/database_commandline/data",1, "Commandline HLATaxonomym: data directory created");
		#Test 12 - Data directory
		is(-d "$outdir/database_commandline/data/reference",1, "Commandline HLATaxonomy: data/reference directory created");
		#Test 13 - Data directory
		is(-d "$outdir/database_commandline/taxonomy",1, "Commandline HLATaxonomy: taxonomy directory created");
		#Test 14 - Successful nodes.dmp
		compareFiles("$nodes", "$outdir/database_commandline/taxonomy/nodes.dmp", "Commandline HLATaxonomy: Successful nodes.dmp");
		#Test 15 - Successful names.dmp
		compareFiles("$names", "$outdir/database_commandline/taxonomy/names.dmp", "Commandline HLATaxonomy: Successful names.dmp");
		#Test 16 - Successful kraken fa
		compareFiles("$indir/outputs/database/data/reference/hla.ref.forKraken.fa", "$outdir/database_commandline/data/reference/hla.ref.forKraken.fa", "Commandline HLATaxonomy: Successful kraken reference fasta");
		
	};
}

sub testHLADistractome{
	subtest 'HLADistractome.pm'=> sub{
		plan tests=>13;
		
		if(-e "$outdir/distractome.function.fa"){
			`rm $outdir/distractome.function.fa`;
		}	
		if(-e "$outdir/distractome.commandline.fa"){
			`rm $outdir/distractome.commandline.fa`;
		}	
		
		#Test 1 -Load module
		use_ok("HLADistractome");
		$loaded_modules{"HLADistractome"}=1;
		
		my @r = trap {HLADistractome::findExcludeGenes("gtf",$exclude_bed)};
		#Test 2 - Invalid gtf exit
		is($trap->exit, 1, "Function findExcludeGenes: invalid transcript gtf exit");
		#Test 3 - Invalid transcript exit		
		is($trap->stderr, "Transcripts GTF gtf does not exist. Exiting\n", "Function findExcludeGenes: invalid transcript gtf error message");

		@r = trap {HLADistractome::findExcludeGenes($transcripts_gtf,"exclude")};
		#Test 4 - Invalid gtf exit
		is($trap->exit, 1, "Function findExcludeGenes: invalid exclude bed exit");
		#Test 5 - Invalid transcript exit		
		is($trap->stderr, "Gene exclusion bed exclude does not exist. Exiting\n", "Function findExcludeGenes: invalid exclude bed error message");

		my %exclude_genes = ("HLA-A"=>1,LST1=>1,EHMT2=>1,HCG4B=>1);
		my $genes = HLADistractome::findExcludeGenes($transcripts_gtf,$exclude_bed);
		#Test 6 - Correct exclude hash
		is_deeply($genes,\%exclude_genes, "Function findExcludeGenes: correct exclude hash");

		@r = trap {HLADistractome::createDistractome("fasta", "$outdir/distractome.function.fa")};
		#Test 7 - Invalid transcripts fa exit
		is($trap->exit, 1, "Function createDistractome: invalid transcripts fasta exit");
		#Test 8 - Invalid transcripts fa error message 		
		is($trap->stderr, "Transcripts fa fasta does not exist. Exiting\n", "Function createDistractome: invalid transcripts fasta error message");
		
		HLADistractome::createDistractome($transcripts_fasta, "$outdir/distractome.function.fa");
		#Test 9 - Correct distractome fasta
		compareFiles("$indir/outputs/database/data/reference/distractome.fa","$outdir/distractome.function.fa", "Function HLADistractome: correct distractome fasta");		
		
		my $out = `perl $FindBin::Bin/../bin/modules/HLADistractome.pm -g $transcripts_gtf -e $exclude_bed -o "$outdir/distractome.commandline.fa" 2>&1`;
		#Test 10 - Missing transcripts fasta
		is($?, 256, "Commandline HLADistractome: missing transcripts fasta");
		
		$out = `perl $FindBin::Bin/../bin/modules/HLADistractome.pm -t $transcripts_fasta -e $exclude_bed -o "$outdir/distractome.commandline.fa" 2>&1`;
		#Test 11 - Missing transcripts gtf
		is($?, 256, "Commandline HLADistractome: missing transcripts gtf");
		
		$out = `perl $FindBin::Bin/../bin/modules/HLADistractome.pm -t $transcripts_fasta -g $transcripts_gtf -o "$outdir/distractome.commandline.fa" 2>&1`;
		#Test 12 - Missing exclude bed
		is($?, 256, "Commandline HLADistractome: missing exclude bed");
	
		$out = `perl $FindBin::Bin/../bin/modules/HLADistractome.pm -t $transcripts_fasta -g $transcripts_gtf -e $exclude_bed -o "$outdir/distractome.commandline.fa" 2>&1`;
		#Test 13 - Correct distractome fasta
		compareFiles("$indir/outputs/database/data/reference/distractome.fa","$outdir/distractome.commandline.fa", "Commandline HLADistractome: correct distractome fasta");		
	};
}

sub testTaxonomyDivisions{
	subtest 'TaxonomyDivisions.pm'=>sub {
		plan tests=> 13;
		
		if(-e "$outdir/divisions.txt"){
			`rm $outdir/divisions.txt`;
		}	
		
		#Test 1-Load Module
		use_ok("TaxonomyDivisions");
		$loaded_modules{"TaxonomyDivisions"}=1;

		my @r = trap {TaxonomyDivisions::createTaxonomyDivisions("names", "$nodes","$outdir/divisions.function.txt")};
		#Test 2- invalid names.dmp exit
		is($trap->exit, 1, "Function createTaxonomyDivisions: invalid names.dmp");
		#Test 3 -invalid names.dmp message
		is($trap->stderr, "Names file names does not exist or is empty. Please specify an existing and non-empty names.dmp file. Exiting.\n", "Function createTaxonomyDivisions: invalid names.dmp error message");
		
		@r = trap {TaxonomyDivisions::createTaxonomyDivisions("", "$nodes","$outdir/divisions.function.txt")};
		#Test 4- missing names.dmp exit
		is($trap->exit, 1, "Function createTaxonomyDivisions: missing names.dmp error");
		#Test 5 -missing names.dmp message
		is($trap->stderr, "Please specify a names.dmp file. Exiting.\n", "Function createTaxonomyDivisions: invalid names.dmp error message");
	
		@r = trap {TaxonomyDivisions::createTaxonomyDivisions("$names","nodes","$outdir/divisions.function.txt")};
		#Test 6- invalid nodes.dmp exit
		is($trap->exit, 1, "Function createTaxonomyDivisions: invalid nodes.dmp");
		#Test 7- invalid nodes.dmp message
		is($trap->stderr, "Nodes file nodes does not exist or is empty. Please specify an existing and non-empty nodes.dmp file. Exiting.\n", "Function createTaxonomyDivisions: invalid nodes.dmp error message");

		@r = trap {TaxonomyDivisions::createTaxonomyDivisions("$names","","$outdir/divisions.function.txt")};
		#Test 8- missing nodes.dmp exit
		is($trap->exit, 1, "Function createTaxonomyDivisions: missing nodes.dmp");
		#Test 9- invalid nodes.dmp message
		is($trap->stderr, "Please specify a nodes.dmp file. Exiting.\n", "Function createTaxonomyDivisions: missing nodes.dmp error message");

		TaxonomyDivisions::createTaxonomyDivisions("$names","$nodes","$outdir/divisions.function.txt");
		#Test 10 -correct divisions.txt
		compareFiles($divisions, "$outdir/divisions.function.txt", "Function createTaxonomyDivisions: correct divisions.txt");	
		
		my $out = `perl $FindBin::Bin/../bin/modules/TaxonomyDivisions.pm -t $nodes -o "$outdir/divisions.commandline.txt" 2>&1`;
		#Test 11 - Commandline missing names.txt
		is($?, 256, "Commandline TaxonomyDivisions: missing names.dmp");
		
		$out = `perl $FindBin::Bin/../bin/modules/TaxonomyDivisions.pm -n "$names" "$outdir/divisions.commandline.txt" 2>&1`;
		#Test 12 - Commandline missing nodes.txt
		is($?, 256, "Commandline TaxonomyDivisions: missing nodes.dmp");

		$out = `perl $FindBin::Bin/../bin/modules/TaxonomyDivisions.pm -n "$names" -t $nodes -o "$outdir/divisions.commandline.txt" 2>&1`;
		#Test 13 - Commandline correct divisions.txt
		compareFiles($divisions, "$outdir/divisions.commandline.txt", "Commandline createTaxonomyDivisions: correct divisions.txt");	
	};	
}

sub testRunKraken{
	subtest 'RunKraken.pm'=>sub {
		plan tests=>17;
			
		#Test 1 - Load module
		use_ok("RunKraken");
		$loaded_modules{"RunKraken"}=1;
		`mkdir -p $outdir/test_database`;
		`cp -r $indir/outputs/database/taxonomy $outdir/test_database`;
		if(-e "$outdir/test_database/library/added"){
			`rm $outdir/test_database/library/added/*.fna`;
		}
		if(-e "$outdir/test_database/lca.complete"){
			`rm $outdir/test_database/lca.complete`;
		}
		if(-e "$outdir/test_database/database.kdb"){
			`rm $outdir/test_database/database.kdb`;
		}
		if(-e "$outdir/test_database/database.jdb"){
			`rm $outdir/test_database/database.jdb`;
		}
		if(-e "$outdir/test_database/database.idx"){
			`rm $outdir/test_database/database.idx`;
		}
		if(-e "$outdir/HLA00001.commandline.A_1.fq"){
			`rm $outdir/HLA00001.commandline.A_1.fq`;
		}
		if(-e "$outdir/HLA00001.commandline.A_2.fq"){
			`rm $outdir/HLA00001.commandline.A_2.fq`;
		}
		if(-e "$outdir/HLA00001.function.A_1.fq"){
			`rm $outdir/HLA00001.function.A_1.fq`;
		}
		if(-e "$outdir/HLA00001.function.reads.txt"){
			`rm $outdir/HLA00001.function.reads.txt`;
		}
		if(-e "$outdir/HLA00001.function.A_2.fq"){
			`rm $outdir/HLA00001.function.A_2.fq`;
		}
		`echo ">" $outdir/test_database/taxonomy/gi_taxid_nucl.dmp`;
		
		RunKraken::addLibraryToKraken("$outdir","test_database",$distractome,$kraken_path,"$outdir/test.log");
		RunKraken::addLibraryToKraken("$outdir","test_database",$kraken_fasta,$kraken_path,"$outdir/test.log");
		my @list = `ls -t $outdir/test_database/library/added/*.fna`; 
		chomp($list[0]);
		chomp($list[1]);
		#Test 2 - Function addLibraryToKraken
		compareFiles($kraken_fasta,$list[0],"Function addLibraryToKraken: reference added");	
		#Test 3 - Function addLibraryToKraken
		compareFiles($distractome, $list[1],"Function addLibraryToKraken: distractome added");
		
		RunKraken::buildKraken("$outdir","test_database","--minimizer-len 3 --kmer-len 31 --jellyfish-hash-size 500000",1,$kraken_path,"$outdir/test.log");
		#Test 4 - Function buildKraken
		is(compare("$outdir/test_database/database.kdb","$indir/outputs/database/database.kdb"),0,"Function buildKraken: successful database.kdb");
		#Test 5 - Function buildKraken
		is(compare("$outdir/test_database/database.idx","$indir/outputs/database/database.idx"),0,"Function buildKraken: successful database.idx");
		
		RunKraken::filterReads("$outdir","test_database","$outdir/HLA00001.function",1,"$indir/inputs/HLA00001_1.fastq.gz","$indir/inputs/HLA00001_2.fastq.gz",$kraken_path,"$outdir/HLA00001.function.reads.txt","$outdir/test.log");
		#Test 6 - Function filterReads
		compareFiles("$indir/outputs/filtered/HLA00001.A_1.fq","$outdir/HLA00001.function.A_1.fq","Function filterReads: correct Read1");
		#Test 7 - Function filterReads
		compareFiles("$indir/outputs/filtered/HLA00001.A_2.fq","$outdir/HLA00001.function.A_2.fq","Function filterReads: correct Read2");
		#Test 8 - Function filterReads
		compareFiles("$indir/outputs/filtered/HLA00001.reads.txt","$outdir/HLA00001.function.reads.txt","Function filterReads: same classifications");
		
		my $out = `perl $FindBin::Bin/../bin/modules/RunKraken.pm 2>&1`;
		#Test 9 - Commandline missing -mode
		is($?, 256, "Commandline RunKraken: missing mode");
	
		$out = `perl $FindBin::Bin/../bin/modules/RunKraken.pm -m build -dd $outdir 2>&1`;
		#Test 10 - Commandline missing -database_name
		is($?, 256, "Commandline RunKraken: missing database_name");
		
		$out = `perl $FindBin::Bin/../bin/modules/RunKraken.pm -m filter  -dd $outdir -db test_database 2>&1`;
		#Test 11 - Commandline missing fastq1
		is($?, 256, "Commandline RunKraken filter: missing fastq1");
		
		$out = `perl $FindBin::Bin/../bin/modules/RunKraken.pm -m filter  -dd $outdir -db test_database -fq1 $indir/inputs/HLA00001_1.fastq.gz 2>&1`;
		#Test 12 - Commandline missing fastq2
		is($?, 256, "Commandline RunKraken filter: missing fastq2");
		
		$out = `perl $FindBin::Bin/../bin/modules/RunKraken.pm -m filter -dd $outdir -db test_database -fq1 $indir/inputs/HLA00001_1.fastq.gz -fq2 $indir/inputs/HLA00001_2.fastq.gz 2>&1`;
		#Test 13 - Commandline missing output_prefix
		is($?, 256, "Commandline RunKraken filter: missing output_prefix");
		
		$out = `perl $FindBin::Bin/../bin/modules/RunKraken.pm -m filter -dd $outdir -db test_database -fq1 $indir/inputs/HLA00001_1.fastq.gz -fq2 $indir/inputs/HLA00001_2.fastq.gz -output_prefix $outdir/HLA00001.commandline -p $outdir/HLA00001.commandline.reads.txt -c 1 -kp $kraken_path -l $outdir/test.log 2>&1`;
		#Test 14 - Commandline successful filter
		compareFiles("$indir/outputs/filtered/HLA00001.A_1.fq","$outdir/HLA00001.commandline.A_1.fq","Function filterReads: correct Read1");
		#Test 15 - Commandline successful filter
		compareFiles("$indir/outputs/filtered/HLA00001.A_2.fq","$outdir/HLA00001.commandline.A_2.fq","Function filterReads: correct Read2");
		#Test 16 - Commandline successful filter
		compareFiles("$indir/outputs/filtered/HLA00001.reads.txt","$outdir/HLA00001.commandline.reads.txt","Function filterReads: same classifications");
	
		$out = `perl $FindBin::Bin/../bin/modules/RunKraken.pm -m library 2>&1`;
		#Test17	- Commandline missing fa
		is($?, 256, "Commandline RunKraken library: missing fa");
		
	};
}

sub testSimulateReads{
	subtest 'SimulateReads.pm'=>sub{
		plan tests=>17;
	
		if(-e "$outdir/test_simulation.commandline.HLA00001_1.fastq"){
			`rm $outdir/test_simulation.commandline.HLA00001_1.fastq`;
		}
		if(-e "$outdir/test_simulation.commandline.HLA00001_2.fastq"){
			`rm $outdir/test_simulation.commandline.HLA00001_2.fastq`;
		}
		if(-e "$outdir/test_simulation.commandline.HLA00002_1.fastq"){
			`rm $outdir/test_simulation.commandline.HLA00002_1.fastq`;
		}
		if(-e "$outdir/test_simulation.commandline.HLA00002_2.fastq"){
			`rm $outdir/test_simulation.commandline.HLA00002_2.fastq`;
		}
		if(-e "$outdir/test_simulation.HLA00001_1.fastq"){
			`rm $outdir/test_simulation.HLA00001_1.fastq`;
		}
		if(-e "$outdir/test_simulation.HLA00001_2.fastq"){
			`rm $outdir/test_simulation.HLA00001_2.fastq`;
		}
		if(-e "$outdir/test_simulation.HLA00002_1.fastq"){
			`rm $outdir/test_simulation.HLA00002_1.fastq`;
		}
		if(-e "$outdir/test_simulation.HLA00002_2.fastq"){
			`rm $outdir/test_simulation.HLA00002_2.fastq`;
		}
		if(-e "$outdir/test_simulation.TestSimulation.commandline_1.fastq"){
			`rm $outdir/test_simulation.TestSimulation.commandline_1.fastq`;
		}
		if(-e "$outdir/test_simulation.TestSimulation.commandline_2.fastq"){
			`rm $outdir/test_simulation.TestSimulation.commandline_2.fastq`;
		}
		if(-e "$outdir/test_simulation.TestSimulation.function_1.fastq"){
			`rm $outdir/test_simulation.TestSimulation.function_1.fastq`;
		}
		if(-e "$outdir/test_simulation.TestSimulation.function_2.fastq"){
			`rm $outdir/test_simulation.TestSimulation.function_2.fastq`;
		}
		
			
		#Test 1 - Load Module
		use_ok("SimulateReads");
		$loaded_modules{"SimulateReads"}=1;

		my %expected = (reference=>"$indir/inputs/simulation.ref.fa", numReads=>100, readLength=>75, max_insert=>200, scale=>85, shape=>(1/.75),seed=>321, threads=>1, scripts_dir=>"$FindBin::Bin/../bin",sequence=>"", name=>"", outPrefix=>"$outdir/test_simulation"); 
		my $opts = SimulateReads::setSimulationOptions("$indir/inputs/simulation.ref.fa",100, 75,200, 85, .75,321, 1, "$FindBin::Bin/../bin", "$outdir/test_simulation"); 
		#Test 2 - Set simulation options
		is_deeply($opts, \%expected, "Function setSimulationOptions");

		SimulateReads::makeSim("TestSimulation.function","AGGATTTCGTGTTCCAGTTTAAGGGCATGTGCTACTTCACCAACGGGACGGAGCGCGTGCGTCTTGTGACCAGATACATCTATAACCGAGAGGAGTACGCGCGCTTCGACAGCGACGTGGGGGTGTACCGCGCGGTGACGCCGCAGGGGCGGCCTGATGCCGAGTACTGGAACAGCCAGAAGGAAGTCCTGGAGGGGACCCGGGCGGAGTTGGACACGATGTGCAGACACAACTACGAGGTGGCGTTCCGCGGGATCTTGCAGAGGAGAG", "$outdir/test_simulation");
		#Test 3 - function makeSim read1
		compareFiles("$indir/outputs/simulation.TestSimulation_1.fastq","$outdir/test_simulation.TestSimulation.function_1.fastq", "Function makeSim: correct read1");	
		#Test 4 - function makeSim read2
		compareFiles("$indir/outputs/simulation.TestSimulation_2.fastq","$outdir/test_simulation.TestSimulation.function_2.fastq", "Function makeSim: correct read2");	

		SimulateReads::simulateReads();
		#Test 5 - function simulateReads sequence1 fastq1
		compareFiles("$indir/outputs/simulation.HLA00001_1.fastq","$outdir/test_simulation.HLA00001_1.fastq", "Function makeSim: correct HLA00001 read1");	
		#Test 6 - function simulateReads sequence1 fastq2
		compareFiles("$indir/outputs/simulation.HLA00001_1.fastq","$outdir/test_simulation.HLA00001_1.fastq", "Function makeSim: correct HLA00001 read2");	
		
		#Test 7 - function simulateReads seuqence2 fastq1
		compareFiles("$indir/outputs/simulation.HLA00002_1.fastq","$outdir/test_simulation.HLA00002_1.fastq", "Function makeSim: correct HLA00002 read1");	
		#Test 8 - function simulateReads sequence2 fastq2
		compareFiles("$indir/outputs/simulation.HLA00002_1.fastq","$outdir/test_simulation.HLA00002_1.fastq", "Function makeSim: correct HLA00002 read2");	
		
		my $out = `perl $FindBin::Bin/../bin/modules/SimulateReads.pm 2>&1`;
		#Test 9 - Commandline no name, sequence or reference
		is($?, 256, "Commandline SimulateReads: missing reference, sequence and name");
				
		$out = `perl $FindBin::Bin/../bin/modules/SimulateReads.pm -name TestSequence 2>&1`;
		#Test 10 - Commandline no name, yes sequence 
		is($?, 256, "Commandline SimulateReads: name specified but missing sequence");
				
		$out = `perl $FindBin::Bin/../bin/modules/SimulateReads.pm -sequence ATCGCGA 2>&1`;
		#Test 11 - Commandline yes name, no sequence
		is($?, 256, "Commandline SimulateReads: sequence specified but missing name");

		$out = `perl $FindBin::Bin/../bin/modules/SimulateReads.pm -name TestSimulation.commandline -sequence AGGATTTCGTGTTCCAGTTTAAGGGCATGTGCTACTTCACCAACGGGACGGAGCGCGTGCGTCTTGTGACCAGATACATCTATAACCGAGAGGAGTACGCGCGCTTCGACAGCGACGTGGGGGTGTACCGCGCGGTGACGCCGCAGGGGCGGCCTGATGCCGAGTACTGGAACAGCCAGAAGGAAGTCCTGGAGGGGACCCGGGCGGAGTTGGACACGATGTGCAGACACAACTACGAGGTGGCGTTCCGCGGGATCTTGCAGAGGAGAG -outPrefix $outdir/test_simulation -numReads 100 -readLength 75 -max_insert 200 -scale 85 -shape .75 -seed 321 -threads 1 -scripts_dir "$FindBin::Bin/../bin" 2>&1`;
		#Test 12 - Commandline simulateReads sequence and name fastq1
		compareFiles("$indir/outputs/simulation.TestSimulation_1.fastq","$outdir/test_simulation.TestSimulation.commandline_1.fastq", "Function makeSim: correct read1");	
		#Test 13 - Commandline simulateReads sequence and name fastq2
		compareFiles("$indir/outputs/simulation.TestSimulation_2.fastq","$outdir/test_simulation.TestSimulation.commandline_2.fastq", "Function makeSim: correct read2");	
		
		$out = `perl $FindBin::Bin/../bin/modules/SimulateReads.pm -reference $indir/inputs/simulation.ref.fa -outPrefix $outdir/test_simulation.commandline -numReads 100 -readLength 75 -max_insert 200 -scale 85 -shape .75 -seed 321 -threads 1 -scripts_dir "$FindBin::Bin/../bin" 2>&1`;
		#Test 14 - Commandline simulateReads sequence1 fastq1
		compareFiles("$indir/outputs/simulation.HLA00001_1.fastq","$outdir/test_simulation.commandline.HLA00001_1.fastq", "Function makeSim: correct HLA00001 read1");	
		#Test 15 - Commandline simulateReads sequence1 fastq2
		compareFiles("$indir/outputs/simulation.HLA00001_1.fastq","$outdir/test_simulation.commandline.HLA00001_1.fastq", "Function makeSim: correct HLA00001 read2");	
		
		#Test 16 - Commandline simulateReads seuqence2 fastq1
		compareFiles("$indir/outputs/simulation.HLA00002_1.fastq","$outdir/test_simulation.commandline.HLA00002_1.fastq", "Function makeSim: correct HLA00002 read1");	
		#Test 17 - Commandline simulateReads sequence2 fastq2
		compareFiles("$indir/outputs/simulation.HLA00002_1.fastq","$outdir/test_simulation.commandline.HLA00002_1.fastq", "Function makeSim: correct HLA00002 read2");	
		
	};
}

sub testReadCounter{
	subtest 'ReadCounter.pm'=>sub{
		plan tests=>9;

		if(-e "$outdir/HLA00001_1.function.uniq.cnts"){
			`rm $outdir/HLA00001_1.function.uniq.cnts`;
		}
		if(-e "$outdir/HLA00001_2.function.uniq.cnts"){
			`rm $outdir/HLA00001_2.function.uniq.cnts`;
		}

		#Test 1 - Load module
		use_ok("ReadCounter");
		$loaded_modules{"ReadCounter"}=1;

		my @r = trap {ReadCounter::count_reads("input","$outdir/HLA00001_1.function.uniq.cnts")};
		#Test 2- invalid input
		is($trap->exit, 1, "Function count_reads: invalid input exit");
		#Test 3- invalid nodes.dmp message
		is($trap->stderr, "Unable to open input for reading\n", "Function count_reads: missing nodes.dmp error message");

		ReadCounter::count_reads("$indir/inputs/HLA00001_1.fastq.gz", "$outdir/HLA00001_1.function.uniq.cnts");
		#Test 4 - Function countReads correct output read1
		compareFiles("$indir/outputs/HLA00001_1.uniq.cnts","$outdir/HLA00001_1.function.uniq.cnts", "Function count_reads: correct read1 counts");
		
		ReadCounter::count_reads("$indir/inputs/HLA00001_2.fastq.gz", "$outdir/HLA00001_2.function.uniq.cnts");
		#Test 5 - Function countReads correct output read2
		compareFiles("$indir/outputs/HLA00001_2.uniq.cnts","$outdir/HLA00001_2.function.uniq.cnts", "Function count_reads: correct read2 counts");
		
		my $out = `$FindBin::Bin/../bin/modules/ReadCounter.pm -o $outdir/HLA00001_1.commandline.uniq.cnts 2>&1`;
		#Test 6 - Commandline no input
		is($?, 256, "Commandline ReadCounter: no input provided");
		
		$out = `$FindBin::Bin/../bin/modules/ReadCounter.pm -i $indir/inputs/HLA00001_1.fastq.gz 2>&1`;
		#Test 7 - Commandline no output
		is($?, 256, "Commandline ReadCounter: no output provided");
		
		$out = `$FindBin::Bin/../bin/modules/ReadCounter.pm -i $indir/inputs/HLA00001_1.fastq.gz -o $outdir/HLA00001_1.commandline.uniq.cnts`;
		##Test 8 - Commandline correct output read1
		compareFiles("$indir/outputs/HLA00001_1.uniq.cnts","$outdir/HLA00001_1.commandline.uniq.cnts", "Commandline ReadCounter: correct read1 counts");
		
		$out = `$FindBin::Bin/../bin/modules/ReadCounter.pm -i $indir/inputs/HLA00001_2.fastq.gz -o $outdir/HLA00001_2.commandline.uniq.cnts`;
		##Test 9 - Commandline correct output read2
		compareFiles("$indir/outputs/HLA00001_2.uniq.cnts","$outdir/HLA00001_2.commandline.uniq.cnts", "Commandline ReadCounter: correct read2 counts");
		
	};
}

sub testDetermineProfile{
	subtest 'DetermineProfile.pm'=>sub{
		plan tests=>9;
		#Test 1 - Load module
		use_ok("DetermineProfile");
		$loaded_modules{"DetermineProfile"}=1;
		
		`mkdir -p $outdir/profiles_function`;
		if(-e "$outdir/profiles_function/A.profile.ph"){
			`rm $outdir/profiles_function/A.profile.ph`;
		}	
		if(-e "$outdir/profiles_function/A.profile.txt"){
			`rm $outdir/profiles_function/A.profile.txt`;
		}	

		`mkdir -p $outdir/profiles_commandline`;
		if(-e "$outdir/profiles_commandline/A.profile.ph"){
			`rm $outdir/profiles_commandline/A.profile.ph`;
		}	
		if(-e "$outdir/profiles_commandline/A.profile.txt"){
			`rm $outdir/profiles_commandline/A.profile.txt`;
		}	

		my %sig = ("0"=>159,"1"=>158,"2"=>157,"3"=>156,"4"=>155,"5"=>155,"6"=>154,"7"=>153,"8"=>152,"9"=>152,"10"=>165,"11"=>163,"12"=>160,"13"=>159,"14"=>158,"15"=>157,"16"=>157,"17"=>156,"18"=>155,"19"=>154);
		my %profile = ();
		my %sigs =("HLA00001"=>\%sig);
		$profile{sig}=\%sigs;
		my %kmers = (	"CGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCC"=>10,
				"CATCCAGATAATGTATGGCTGCGACGTGGGGCCGGACGGGCGCTTCCTCC"=>11,
				"CCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGG"=>12,
				"AGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCA"=>13,
				"ACCCTGCGCGGCTACTACAACCAGAGCGAGGACGGTTCTCACACCATCCA"=>14,
				"CACGGAATATGAAGGCCCACTCACAGACTGACCGAGCGAACCTGGGGACC"=>15,
				"GGCCGGAGTATTGGGACCAGGAGACACGGAATATGAAGGCCCACTCACAG"=>16,
				"CAGACCTGGGCGGGCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTC"=>17,
				"TGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATAT"=>18,
				"CGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTT"=>19,
				"GATCACTGGAGCTGTGGTCGCTGCCGTGATGTGGAGGAGGAAGAGCTCAG"=>0,
				"GCTGTCTTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGG"=>1,
				"GAGGAAGAGCTCAGATAGAAAAGGAGGGAGTTACACTCAGGCTGCAAGCA"=>2,
				"TGACCCACCACCCCATCTCTGACCATGAGGCCACCCTGAGGTGCTGGGCC"=>3,
				"GGGAGGACCAGACCCAGGACACGGAGCTCGTGGAGACCAGGCCTGCAGGG"=>4,
				"TTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCC"=>5,
				"CAAGCCCCTCACCCTGAGATGGGAGCTGTCTTCCCAGCCCACCATCCCCA"=>6,
				"CCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCCTT"=>7,
				"ACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCCTTGGAGCTGT"=>8,
				"TGGGAGCTGTCTTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGG"=>9);
		$profile{kmers}= \%kmers;	
		my %names = (	"10"=>"CGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCC",
				"11"=>"CATCCAGATAATGTATGGCTGCGACGTGGGGCCGGACGGGCGCTTCCTCC",
				"12"=>"CCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGG",
				"13"=>"AGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCA",
				"14"=>"ACCCTGCGCGGCTACTACAACCAGAGCGAGGACGGTTCTCACACCATCCA",
				"15"=>"CACGGAATATGAAGGCCCACTCACAGACTGACCGAGCGAACCTGGGGACC",
				"16"=>"GGCCGGAGTATTGGGACCAGGAGACACGGAATATGAAGGCCCACTCACAG",
				"17"=>"CAGACCTGGGCGGGCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTC",
				"18"=>"TGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATAT",
				"19"=>"CGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTT",
				"0"=>"GATCACTGGAGCTGTGGTCGCTGCCGTGATGTGGAGGAGGAAGAGCTCAG",
				"1"=>"GCTGTCTTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGG",
				"2"=>"GAGGAAGAGCTCAGATAGAAAAGGAGGGAGTTACACTCAGGCTGCAAGCA",
				"3"=>"TGACCCACCACCCCATCTCTGACCATGAGGCCACCCTGAGGTGCTGGGCC",
				"4"=>"GGGAGGACCAGACCCAGGACACGGAGCTCGTGGAGACCAGGCCTGCAGGG",
				"5"=>"TTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCC",
				"6"=>"CAAGCCCCTCACCCTGAGATGGGAGCTGTCTTCCCAGCCCACCATCCCCA",
				"7"=>"CCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCCTT",
				"8"=>"ACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCCTTGGAGCTGT",
				"9"=>"TGGGAGCTGTCTTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGG");
		$profile{names}=\%names;
		my %cnts = ("HLA00001"=>3135);
		$profile{totals}=\%cnts;
		
		DetermineProfile::createProfile("A", "$indir/outputs/filtered", "sim", "$outdir/profiles_function", 1,1,"$FindBin::Bin/../bin/");	
		my $loaded_profile = retrieve("$outdir/profiles_function/A.profile.ph");
		#Test 2 - determineProfile correct profile hash file
		is_deeply($loaded_profile, \%profile, "Function createProfile: correct profile hash file");
		#Test 3 - determineProfile correct profile matrix
		compareFiles("$indir/outputs/filtered/A.profile.txt", "$outdir/profiles_function/A.profile.txt", "Function createProfile: correct profile matrix file");
		
		my $observed_hash = DetermineProfile::determineProfile("A", "$indir/outputs/filtered", "sim");
		#Test 4 - determineProfile correct profiles
		is_deeply($observed_hash, \%profile, "Function determineProfile: correct profile");
	
		my $out = `perl $FindBin::Bin/../bin/modules/DetermineProfile.pm -prefix sim -simulation_dir$indir/ outputs/filtered -output_dir $outdir/profiles_commandline -sd $FindBin::Bin/../bin/ -print_hash -print_matrix 2>&1`;
		#Test 5 - Commandline no gene
		is($?, 256, "Commandline DetermineProfile: missing gene");	
	
		$out = `perl $FindBin::Bin/../bin/modules/DetermineProfile.pm -gene A -simulation_dir $indir/outputs/filtered -output_dir $outdir/profiles_commandline -sd $FindBin::Bin/../bin/ -print_hash -print_matrix 2>&1`;
		#Test 6 - Commandline no prefix
		is($?, 256, "Commandline DetermineProfile: missing prefix");	
	
		$out = `perl $FindBin::Bin/../bin/modules/DetermineProfile.pm -gene A -prefix sim -output_dir $outdir/profiles_commandline -sd $FindBin::Bin/../bin/ -print_hash -print_matrix 2>&1`;
		#Test 7 - Commandline no simulation dir
		is($?, 256, "Commandline DetermineProfile: missing simulation directory");	
	
		$out = `perl $FindBin::Bin/../bin/modules/DetermineProfile.pm -gene A -prefix sim -simulation_dir $indir/outputs/filtered -output_dir $outdir/profiles_commandline/ -sd $FindBin::Bin/../bin/ -print_hash -print_matrix 2>&1`;
		#Test 8 - Commandline correct profile hash file
		$loaded_profile = retrieve("$outdir/profiles_commandline/A.profile.ph");
		is_deeply($loaded_profile, \%profile, "Commandline DetermineProfile: correct profile hash file");
		#Test 9 - Commandline correct profile matrix
		compareFiles("$indir/outputs/filtered/A.profile.txt", "$outdir/profiles_commandline/A.profile.txt", "Commandline DetermineProfile: correct profile matrix file");
		
	};
}

sub testSequenceFunctions{
	subtest 'SequenceFunctions.pm'=>sub{
		plan tests=>16;
		
		#Test 1- Load modules
		use_ok("SequenceFunctions");
		$loaded_modules{"SequenceFunctions"}=1;
		
		#Test 2 - Test empty string revcomp
		is(SequenceFunctions::revcomp(""), "", "Function revcomp: empty string");
		
		#Test 3 - Test DNA string revcomp
		is(SequenceFunctions::revcomp("ATCGAGAGATAGTACCCATAGATAGAACGATAGA"), "TCTATCGTTCTATCTATGGGTACTATCTCTCGAT", "Function revcomp: DNA");
		
		#Test 4 - Test RNA string revcomp
		is(SequenceFunctions::revcomp("AGGCCCGGUUCUGACCUGA"), "TCAGGTCAGAACCGGGCCT", "Function revcomp: RNA contains U");
		
		#Test 5 - Test RNA string revcomp 
		is(SequenceFunctions::revcomp("AGGcccGGUUCUGAccUGA"), "TCAggTCAGAACCgggCCT", "Function revcomp: RNA contains U lowercase");
		
		#Test 6 - Test DNA with N string revcomp
		is(SequenceFunctions::revcomp("ATCGANGAGATAGTACCCNATAGATAGAACGATAGA"), "TCTATCGTTCTATCTATNGGGTACTATCTCNTCGAT", "Function revcomp: DNA contains N");
		
		#Test 7 - Test random string string revcomp
		is(SequenceFunctions::revcomp("INCONCEIVABLE"), "ELBTVIEGNOGNI", "Function revcomp: random string");
		
		#Test 8 - Test empty string rev
		is(SequenceFunctions::rev(""), "", "Function rev: empty string");
		
		#Test 9 - Test rev DNA
		is(SequenceFunctions::rev("TATCGAGACGAGGACGATAGACGAGAGATA"), "ATAGAGAGCAGATAGCAGGAGCAGAGCTAT", "Function rev: DNA");
		
		#Test 10 - Test rev Quality string
		is(SequenceFunctions::rev("IIBBII^^^III"), "III^^^IIBBII", "Function rev: Quality string");
		
		#Test 11 - Random string 
		is(SequenceFunctions::rev("May the force be with you"), "uoy htiw eb ecrof eht yaM", "Function rev: random string");

		my $out = `perl $FindBin::Bin/../bin/modules/SequenceFunctions.pm -rv -rc 2>&1`;
		#Test 12 - Commandline no sequence
		is($?, 256, "Commandline SequenceFunctions: missing sequence");	
		
		$out = `perl $FindBin::Bin/../bin/modules/SequenceFunctions.pm -seq ATCGANGAGATAGTACCCNATAGATAGAACGATAGA 2>&1`;
		#Test 13 - Commandline no operation
		is($?, 256, "Commandline SequenceFunctions: missing operations");

		$out = `perl $FindBin::Bin/../bin/modules/SequenceFunctions.pm -rc -seq ATCGANGAGATAGTACCCNATAGATAGAACGATAGA 2>&1`;
		#Test 14 - Commandline revcomp
		is($out, "TCTATCGTTCTATCTATNGGGTACTATCTCNTCGAT\n", "Commandline SequenceFunctions: revcomp");	

		$out = `perl $FindBin::Bin/../bin/modules/SequenceFunctions.pm -rv -seq ATCGANGAGATAGTACCCNATAGATAGAACGATAGA 2>&1`;
		#Test 15 - Commandline rev
		is($out, "AGATAGCAAGATAGATANCCCATGATAGAGNAGCTA\n", "Commandline SequenceFunctions: rev");	

		$out = `perl $FindBin::Bin/../bin/modules/SequenceFunctions.pm -rc -rv -seq ATCGANGAGATAGTACCCNATAGATAGAACGATAGA 2>&1`;
		#Test 16 - Commandline no sequence
		is($out, "TCTATCGTTCTATCTATNGGGTACTATCTCNTCGAT\nAGATAGCAAGATAGATANCCCATGATAGAGNAGCTA\n", "Commandline SequenceFunctions: revcomp and rev");	
			
	};
}

sub testPairPicker{
	subtest 'PairPicker.pm' =>sub{
		plan tests=>14;
	
		#Test- 1 load module
		use_ok("PairPicker");	
		$loaded_modules{"PairPicker"}=1;
	
			
		my @r = trap{PairPicker::runModule("prefix", "$indir/outputs/hla.ref.merged.fa", "HLA00001|HLA00001,HLA00002|HLA00002", "A", 20, 0, "$FindBin::Bin/../bin/", 0)};	
		#Test 2 - invalid prefix
		is($trap->exit, 1, "Function runModule: invalid prefix exit");
		#Test 3 - invalid prefix message
		is($trap->stderr, "Cannot find fastq file with prefix prefix.A\n", "Function runModule: invalid prefix error message");
		
		@r = trap{PairPicker::runModule("$indir/inputs/PairPicker", "$indir/outputs/hla.ref.merged.fa", "", "A", 20, 0, "$FindBin::Bin/../bin/", 0)};	
		#Test 3 - empty candidates
		is($trap->exit, 1, "Function runModule: empty candidates exit");
		#Test 5 - empty candidates  message
		is($trap->stderr, "Candidate list is empty!\n", "Function runModule: empty candidates error message");
	
		my @cc1 = ([0,2.86666666666667],[3.33333333333333,0]);
		my @cp1 = ("HLA00001|HLA00001","HLA00002|HLA00002");
		my $urc1 = 6.2;
		my @pp1 = (\@cc1,\@cp1, $urc1);
		my @ref = PairPicker::runModule("$indir/inputs/PairPicker", "$indir/outputs/hla.ref.merged.fa", "HLA00001|HLA00001,HLA00002|HLA00002", "A", 20, 0, "$FindBin::Bin/../bin/", 0);	
		#Test 6 - Success
		is_deeply(\@ref, \@pp1, "Function runModule: correct matrix");
		
		my @cc2 = ([0,2],[3,0]);
		my @cp2 = ("HLA00001|HLA00001","HLA00002|HLA00002");
		my $urc2 = 5;
		my @pp2 = (\@cc2,\@cp2, $urc2);
		@ref = PairPicker::runModule("$indir/inputs/PairPicker", "$indir/outputs/hla.ref.merged.fa", "HLA00001|HLA00001,HLA00002|HLA00002", "A", 35, 0, "$FindBin::Bin/../bin/", 0);	
		#Test 7 - Success min quality
		is_deeply(\@ref, \@pp2, "Function runModule: min quality correct matrix");
		
		my @cc3 = ([0,7.86666666666667],[8.33333333333333,0]);
		my @cp3 = ("HLA00001|HLA00001","HLA00002|HLA00002");
		my $urc3 = 6.2;
		my @pp3 = (\@cc3,\@cp3, $urc3);
		@ref = PairPicker::runModule("$indir/inputs/PairPicker", "$indir/outputs/hla.ref.merged.fa", "HLA00001|HLA00001,HLA00002|HLA00002", "A", 20, 1, "$FindBin::Bin/../bin/", 0);	
		#Test 8 - Success using scale 1
		is_deeply(\@ref, \@pp3, "Function runModule: scale 1 correct matrix");
		
		my @cc4 = ([0,5.36666666666667],[5.83333333333333,0]);
		my @cp4 = ("HLA00001|HLA00001","HLA00002|HLA00002");
		my $urc4 = 6.2;
		my @pp4 = (\@cc4,\@cp4, $urc4);
		@ref = PairPicker::runModule("$indir/inputs/PairPicker", "$indir/outputs/hla.ref.merged.fa", "HLA00001|HLA00001,HLA00002|HLA00002", "A", 20, .5, "$FindBin::Bin/../bin/", 0);	
		#Test 9 - Success using scale .5
		is_deeply(\@ref, \@pp4, "Function runModule: scale .5 correct matrix");
		
		my @cc5 = ([0,5.73333333333333],[6.66666666666667,0]);
		my @cp5 = ("HLA00001|HLA00001","HLA00002|HLA00002");
		my $urc5 = 12.4;
		my @pp5 = (\@cc5,\@cp5, $urc5);
		@ref = PairPicker::runModule("$indir/inputs/PairPicker", "$indir/outputs/hla.ref.merged.fa", "HLA00001|HLA00001,HLA00002|HLA00002", "A", 20, 0, "$FindBin::Bin/../bin/", 1);	
		#Test 10 - Success using ClassI
		is_deeply(\@ref, \@pp5, "Function runModule: use ClassI reads correct matrix");

		my $out = `perl $FindBin::Bin/../bin/modules/PairPicker.pm -g A -o "$outdir/PairPicker.matrix.txt" -c "HLA00001|HLA00001,HLA00002|HLA00002" -r $indir/outputs/hla.ref.merged.fa -sd $FindBin::Bin/../bin/ 2>&1`;
		#Test 11 - Commandline no input
		is($?, 256, "Commandline PairPicker: missing input");
		
		$out = `perl $FindBin::Bin/../bin/modules/PairPicker.pm -i $indir/inputs/PairPicker -o "$outdir/PairPicker.matrix.txt" -c "HLA00001|HLA00001,HLA00002|HLA00002" -r $indir/outputs/hla.ref.merged.fa -sd $FindBin::Bin/../bin/ 2>&1`;
		#Test 12 - Commandline no operation
		is($?, 256, "Commandline PairPicker: missing gene");
		
		$out = `perl $FindBin::Bin/../bin/modules/PairPicker.pm -i $indir/inputs/PairPicker -g A -o "$outdir/PairPicker.matrix.txt" -c "HLA00001|HLA00001,HLA00002|HLA00002" -sd $FindBin::Bin/../bin/ 2>&1`;
		#Test 13 - Commandline no operation
		is($?, 256, "Commandline PairPicker: missing reference");
		
		$out = `perl $FindBin::Bin/../bin/modules/PairPicker.pm -i $indir/inputs/PairPicker -g A -o "$outdir/PairPicker.matrix.txt" -c "HLA00001|HLA00001,HLA00002|HLA00002" -r $indir/outputs/hla.ref.merged.fa -sd $FindBin::Bin/../bin/ 2>&1`;
		#Test 14 - Commandline success
		compareFiles("$indir/outputs/pair_picker.matrix.txt", "$outdir/PairPicker.matrix.txt", "Commandline PairPicker: correct matrix");
	};
}

sub testAlleleRefiner{
	subtest 'AlleleRefiner.pm' =>sub{
		plan tests=>56;

		if(-e "$outdir/HLA00001_HLA00002.coverage.txt"){
			`rm $outdir/HLA00001_HLA00002.coverage.txt`;
		}
		if(-e "$outdir/HLA00001_HLA00002.mismatches.txt"){
			`rm $outdir/HLA00001_HLA00002.mismatches.txt`;
		}
		if(-e "$outdir/HLA00001_HLA00002.extra_kmers.txt"){
			`rm $outdir/HLA00001_HLA00002.extra_kmers.txt`;
		}
		if(-e "$outdir/HLA00001N_HLA00002.coverage.txt"){
			`rm $outdir/HLA00001N_HLA00002.coverage.txt`;
		}
		if(-e "$outdir/HLA00001N_HLA00002.mismatches.txt"){
			`rm $outdir/HLA00001N_HLA00002.mismatches.txt`;
		}
		if(-e "$outdir/HLA00001N_HLA00002.extra_kmers.txt"){
			`rm $outdir/HLA00001N_HLA00002.extra_kmers.txt`;
		}
		if(-e "$outdir/HLA00001_HLA00002.nf.coverage.txt"){
			`rm $outdir/HLA00001_HLA00002.nf.coverage.txt`;
		}
		if(-e "$outdir/HLA00001_HLA00002.nf.mismatches.txt"){
			`rm $outdir/HLA00001_HLA00002.nf.mismatches.txt`;
		}
		if(-e "$outdir/HLA00001_HLA00002.nf.extra_kmers.txt"){
			`rm $outdir/HLA00001_HLA00002.nf.extra_kmers.txt`;
		}
		
		#Test- 1 load module
		use_ok("AlleleRefiner");	
		$loaded_modules{"AlleleRefiner"}=1;
		#Test- 2 load module
		SKIP:{
			skip "Sequence Functions already loaded", 1 unless (!defined $loaded_modules{"SequenceFunctions"});
			use_ok("SequenceFunctions");
			$loaded_modules{"SequenceFunctions"}=1;
		}
	
		my %candidate_truth_hash = ("HLA00001"=>"ACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGT",
				      	    "HLA00002"=>"ACTCCATGAGGTATTTCTCCACATCCGTGTCCCGGCCCGGCAGTGGAGAGCCCCGCTTCATCGCAGT");

		my %reference_truth_hash = ();
		$reference_truth_hash{GAGCCCCGCTTCATCGCAGT}{"HLA00002::48"}=1;
		$reference_truth_hash{AGTGGAGAGCCCCGCTTCAT}{"HLA00002::42"}=1;
		$reference_truth_hash{AGGTATTTCTCCACATCCGT}{"HLA00002::9"}=1;
		$reference_truth_hash{CAGTGGAGAGCCCCGCTTCA}{"HLA00002::41"}=1;
		$reference_truth_hash{GGAGCCCCGCTTCATCGCCG}{"HLA00001::47"}=1;
		$reference_truth_hash{ATCCGTGTCCCGGCCCGGCC}{"HLA00001::23"}=1;
		$reference_truth_hash{GTCCCGGCCCGGCAGTGGAG}{"HLA00002::29"}=1;
		$reference_truth_hash{TCCACATCCGTGTCCCGGCC}{"HLA00002::18"}=1;
		$reference_truth_hash{CGGCCCGGCCGCGGGGAGCC}{"HLA00001::33"}=1;
		$reference_truth_hash{TTTCTCCACATCCGTGTCCC}{"HLA00002::14"}=1;
		$reference_truth_hash{GCCGCGGGGAGCCCCGCTTC}{"HLA00001::40"}=1;
		$reference_truth_hash{CGGGGAGCCCCGCTTCATCG}{"HLA00001::44"}=1;
		$reference_truth_hash{CCGCGGGGAGCCCCGCTTCA}{"HLA00001::41"}=1;
		$reference_truth_hash{TCTTCACATCCGTGTCCCGG}{"HLA00001::16"}=1;
		$reference_truth_hash{CGGCCGCGGGGAGCCCCGCT}{"HLA00001::38"}=1;
		$reference_truth_hash{GGCCCGGCCGCGGGGAGCCC}{"HLA00001::34"}=1;
		$reference_truth_hash{CTTCACATCCGTGTCCCGGC}{"HLA00001::17"}=1;
		$reference_truth_hash{TCTCCACATCCGTGTCCCGG}{"HLA00002::16"}=1;
		$reference_truth_hash{TGTCCCGGCCCGGCAGTGGA}{"HLA00002::28"}=1;
		$reference_truth_hash{CGGCAGTGGAGAGCCCCGCT}{"HLA00002::38"}=1;
		$reference_truth_hash{TCCCGGCCCGGCCGCGGGGA}{"HLA00001::30"}=1;
		$reference_truth_hash{GGAGAGCCCCGCTTCATCGC}{"HLA00002::45"}=1;
		$reference_truth_hash{ATGAGGTATTTCTTCACATC}{"HLA00001::6"}=1;
		$reference_truth_hash{GTATTTCTCCACATCCGTGT}{"HLA00002::11"}=1;
		$reference_truth_hash{TCCATGAGGTATTTCTTCAC}{"HLA00001::3"}=1;
		$reference_truth_hash{GGTATTTCTCCACATCCGTG}{"HLA00002::10"}=1;
		$reference_truth_hash{ATGAGGTATTTCTCCACATC}{"HLA00002::6"}=1;
		$reference_truth_hash{GCGGGGAGCCCCGCTTCATC}{"HLA00001::43"}=1;
		$reference_truth_hash{CCGGCCGCGGGGAGCCCCGC}{"HLA00001::37"}=1;
		$reference_truth_hash{TATTTCTTCACATCCGTGTC}{"HLA00001::12"}=1;
		$reference_truth_hash{CTCCATGAGGTATTTCTCCA}{"HLA00002::2"}=1;
		$reference_truth_hash{AGGTATTTCTTCACATCCGT}{"HLA00001::9"}=1;
		$reference_truth_hash{GCAGTGGAGAGCCCCGCTTC}{"HLA00002::40"}=1;
		$reference_truth_hash{GAGGTATTTCTTCACATCCG}{"HLA00001::8"}=1;
		$reference_truth_hash{TGTCCCGGCCCGGCCGCGGG}{"HLA00001::28"}=1;
		$reference_truth_hash{TATTTCTCCACATCCGTGTC}{"HLA00002::12"}=1;
		$reference_truth_hash{GGCAGTGGAGAGCCCCGCTT}{"HLA00002::39"}=1;
		$reference_truth_hash{CCCGGCCCGGCAGTGGAGAG}{"HLA00002::31"}=1;
		$reference_truth_hash{CACATCCGTGTCCCGGCCCG}{"HLA00001::20"}=1;
		$reference_truth_hash{CACATCCGTGTCCCGGCCCG}{"HLA00002::20"}=1;
		$reference_truth_hash{CCGGCCCGGCAGTGGAGAGC}{"HLA00002::32"}=1;
		$reference_truth_hash{TCCGTGTCCCGGCCCGGCAG}{"HLA00002::24"}=1;
		$reference_truth_hash{GTGGAGAGCCCCGCTTCATC}{"HLA00002::43"}=1;
		$reference_truth_hash{GTCCCGGCCCGGCCGCGGGG}{"HLA00001::29"}=1;
		$reference_truth_hash{TGAGGTATTTCTCCACATCC}{"HLA00002::7"}=1;
		$reference_truth_hash{CCGTGTCCCGGCCCGGCCGC}{"HLA00001::25"}=1;
		$reference_truth_hash{CTCCACATCCGTGTCCCGGC}{"HLA00002::17"}=1;
		$reference_truth_hash{TCCGTGTCCCGGCCCGGCCG}{"HLA00001::24"}=1;
		$reference_truth_hash{GGTATTTCTTCACATCCGTG}{"HLA00001::10"}=1;
		$reference_truth_hash{GAGAGCCCCGCTTCATCGCA}{"HLA00002::46"}=1;
		$reference_truth_hash{ATCCGTGTCCCGGCCCGGCA}{"HLA00002::23"}=1;
		$reference_truth_hash{ATTTCTTCACATCCGTGTCC}{"HLA00001::13"}=1;
		$reference_truth_hash{AGAGCCCCGCTTCATCGCAG}{"HLA00002::47"}=1;
		$reference_truth_hash{CCGGCAGTGGAGAGCCCCGC}{"HLA00002::37"}=1;
		$reference_truth_hash{GGCCCGGCAGTGGAGAGCCC}{"HLA00002::34"}=1;
		$reference_truth_hash{CCGTGTCCCGGCCCGGCAGT}{"HLA00002::25"}=1;
		$reference_truth_hash{TGAGGTATTTCTTCACATCC}{"HLA00001::7"}=1;
		$reference_truth_hash{ACATCCGTGTCCCGGCCCGG}{"HLA00001::21"}=1;
		$reference_truth_hash{ACATCCGTGTCCCGGCCCGG}{"HLA00002::21"}=1;
		$reference_truth_hash{CATGAGGTATTTCTCCACAT}{"HLA00002::5"}=1;
		$reference_truth_hash{GCCCGGCCGCGGGGAGCCCC}{"HLA00001::35"}=1;
		$reference_truth_hash{ATTTCTCCACATCCGTGTCC}{"HLA00002::13"}=1;
		$reference_truth_hash{CCCGGCCCGGCCGCGGGGAG}{"HLA00001::31"}=1;
		$reference_truth_hash{GTGTCCCGGCCCGGCCGCGG}{"HLA00001::27"}=1;
		$reference_truth_hash{CCATGAGGTATTTCTCCACA}{"HLA00002::4"}=1;
		$reference_truth_hash{CCCGGCAGTGGAGAGCCCCG}{"HLA00002::36"}=1;
		$reference_truth_hash{CCCGGCCGCGGGGAGCCCCG}{"HLA00001::36"}=1;
		$reference_truth_hash{ACTCCATGAGGTATTTCTCC}{"HLA00002::1"}=1;
		$reference_truth_hash{CTCCATGAGGTATTTCTTCA}{"HLA00001::2"}=1;
		$reference_truth_hash{GTGTCCCGGCCCGGCAGTGG}{"HLA00002::27"}=1;
		$reference_truth_hash{TTCACATCCGTGTCCCGGCC}{"HLA00001::18"}=1;
		$reference_truth_hash{GTATTTCTTCACATCCGTGT}{"HLA00001::11"}=1;
		$reference_truth_hash{TTCTTCACATCCGTGTCCCG}{"HLA00001::15"}=1;
		$reference_truth_hash{GGGGAGCCCCGCTTCATCGC}{"HLA00001::45"}=1;
		$reference_truth_hash{GCCCGGCAGTGGAGAGCCCC}{"HLA00002::35"}=1;
		$reference_truth_hash{CCGGCCCGGCCGCGGGGAGC}{"HLA00001::32"}=1;
		$reference_truth_hash{GGGAGCCCCGCTTCATCGCC}{"HLA00001::46"}=1;
		$reference_truth_hash{TCCCGGCCCGGCAGTGGAGA}{"HLA00002::30"}=1;
		$reference_truth_hash{CATCCGTGTCCCGGCCCGGC}{"HLA00001::22"}=1;
		$reference_truth_hash{CATCCGTGTCCCGGCCCGGC}{"HLA00002::22"}=1;
		$reference_truth_hash{CATGAGGTATTTCTTCACAT}{"HLA00001::5"}=1;
		$reference_truth_hash{CGCGGGGAGCCCCGCTTCAT}{"HLA00001::42"}=1;
		$reference_truth_hash{GAGCCCCGCTTCATCGCCGT}{"HLA00001::48"}=1;
		$reference_truth_hash{GAGGTATTTCTCCACATCCG}{"HLA00002::8"}=1;
		$reference_truth_hash{CGGCCCGGCAGTGGAGAGCC}{"HLA00002::33"}=1;
		$reference_truth_hash{TCACATCCGTGTCCCGGCCC}{"HLA00001::19"}=1;
		$reference_truth_hash{ACTCCATGAGGTATTTCTTC}{"HLA00001::1"}=1;
		$reference_truth_hash{GGCCGCGGGGAGCCCCGCTT}{"HLA00001::39"}=1;
		$reference_truth_hash{CGTGTCCCGGCCCGGCAGTG}{"HLA00002::26"}=1;
		$reference_truth_hash{CGTGTCCCGGCCCGGCCGCG}{"HLA00001::26"}=1;
		$reference_truth_hash{TTCTCCACATCCGTGTCCCG}{"HLA00002::15"}=1;
		$reference_truth_hash{CCACATCCGTGTCCCGGCCC}{"HLA00002::19"}=1;
		$reference_truth_hash{CCATGAGGTATTTCTTCACA}{"HLA00001::4"}=1;
		$reference_truth_hash{TCCATGAGGTATTTCTCCAC}{"HLA00002::3"}=1;
		$reference_truth_hash{TGGAGAGCCCCGCTTCATCG}{"HLA00002::44"}=1;
		$reference_truth_hash{TTTCTTCACATCCGTGTCCC}{"HLA00001::14"}=1;
		my %reference_hash = ();		
		my %candidate_hash = ();
		
		AlleleRefiner->addAlleleToReference("HLA00001","ACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGT",20,\%reference_hash, \%candidate_hash);
		AlleleRefiner->addAlleleToReference("HLA00002","ACTCCATGAGGTATTTCTCCACATCCGTGTCCCGGCCCGGCAGTGGAGAGCCCCGCTTCATCGCAGT",20,\%reference_hash, \%candidate_hash);
		

		##Test 3 - Add allele reference correct reference hash	
		is_deeply(\%reference_hash, \%reference_truth_hash, "Function addAlleleToReference: correct reference hash");	
		##Test 4 - Add allele reference corerct candidate reference hash	
		is_deeply(\%candidate_hash, \%candidate_truth_hash, "Function addAlleleToReference: correct candidate reference hash");
		my %alignment1_truth_hash = (0=>{alignments=>{"HLA00001"=>{12=>{}}},"counts"=>{"HLA00001"=>1},ref_counts=>{"HLA00001"=>{12=>{12=>1,13=>1,14=>1,15=>1,16=>1,17=>1,18=>1,19=>1,20=>1,21=>1,22=>1,23=>1,24=>1,25=>1,26=>1,27=>1,28=>1,29=>1,30=>1,31=>1,32=>1,33=>1,34=>1,35=>1,36=>1}}},scores=>{"HLA00001"=>{12=>{12=>1,13=>1,14=>1,15=>1,16=>1,17=>1,18=>1,19=>1,20=>1,21=>1,22=>1,23=>1,24=>1,25=>1,26=>1,27=>1,28=>1,29=>1,30=>1,31=>1,32=>1,33=>1,34=>1,35=>1,36=>1}}}});
		my %alignments1 = ();
		my %pos1_hash=(12=>1);
		my $mm =  AlleleRefiner->identifyMismatches("HLA00001", 25, "TATTTCTTCACATCCGTGTCCCGGC", \%pos1_hash,\%alignments1, 1, "IIIIIIIIIIIIIIIIIIIIIIIII", \%candidate_truth_hash, 1,20);	
		##Test 5 - identifyMismatches correct alignment hash	
		is_deeply(\%alignments1, \%alignment1_truth_hash, "Function identifyMismatches: 0 mismatch correct alignment hash");
		##Test 6 - identifyMismatches correct minimum mismatch
		is($mm,0,"Function identifyMismatches: 0 mismatch correct minimum mismatch");
		
		my %alignment2_truth_hash = (1=>{alignments=>{"HLA00001"=>{12=>{34=>'A'}}},"counts"=>{"HLA00001"=>1},ref_counts=>{"HLA00001"=>{12=>{12=>1,13=>1,14=>1,15=>1,16=>1,17=>1,18=>1,19=>1,20=>1,21=>1,22=>1,23=>1,24=>1,25=>1,26=>1,27=>1,28=>1,29=>1,30=>1,31=>1,32=>1,33=>1,35=>1,36=>1}}},scores=>{"HLA00001"=>{12=>{12=>1,13=>1,14=>1,15=>1,16=>1,17=>1,18=>1,19=>1,20=>1,21=>1,22=>1,23=>1,24=>1,25=>1,26=>1,27=>1,28=>1,29=>1,30=>1,31=>1,32=>1,33=>1,34=>1,35=>1,36=>1}}}});
		my %alignments2 = ();
		$mm =  AlleleRefiner->identifyMismatches("HLA00001", 25, "TATTTCTTCACATCCGTGTCCCAGC", \%pos1_hash,\%alignments2, 1, "IIIIIIIIIIIIIIIIIIIIIIIII", \%candidate_truth_hash, 1,20);	
		##Test 7 - identifyMismatches correct alignment hash	
		is_deeply(\%alignments2, \%alignment2_truth_hash, "Function identifyMismatches: 1 mismatch correct alignment hash");
		##Test 8 - identifyMismatches correct minimum mismatch
		is($mm,1,"Function identifyMismatches: 1 mismatch correct minimum mismatch");
		
		my %alignment3_truth_hash = (1=>{alignments=>{"HLA00001"=>{12=>{34=>'A'}}},"counts"=>{"HLA00001"=>1},ref_counts=>{"HLA00001"=>{12=>{12=>1,13=>1,14=>1,15=>1,16=>1,17=>1,18=>1,19=>1,20=>1,21=>1,22=>1,23=>1,24=>1,25=>1,26=>1,27=>1,28=>1,29=>1,30=>1,31=>1,32=>1,33=>1,35=>1,36=>1}}},scores=>{"HLA00001"=>{12=>{12=>1,13=>1,14=>1,15=>1,16=>1,17=>1,18=>1,19=>1,20=>1,21=>1,22=>1,23=>1,24=>1,25=>1,26=>1,27=>1,28=>1,29=>1,30=>1,31=>1,32=>1,33=>1,34=>0,35=>1,36=>1}}}});
		my %alignments3 = ();
		$mm =  AlleleRefiner->identifyMismatches("HLA00001", 25, "TATTTCTTCACATCCGTGTCCCAGC", \%pos1_hash,\%alignments3, 1, "IIIIIIIIIIIIIIIIIIIIII+II", \%candidate_truth_hash, 1,20);	
		##Test 9 - identifyMismatches correct alignment hash	
		is_deeply(\%alignments3, \%alignment3_truth_hash, "Function identifyMismatches: 1 mismatch Q10 correct alignment hash");
		##Test 10 - identifyMismatches correct minimum mismatch
		is($mm,1,"Function identifyMismatches: 1 mismatch Q10 correct minimum mismatch");
	
		my %alignment4_truth_hash = (1=>{alignments=>{"HLA00001"=>{12=>{34=>'A'}}},"counts"=>{"HLA00001"=>1},ref_counts=>{"HLA00001"=>{12=>{12=>1,13=>1,14=>1,15=>1,16=>1,17=>1,18=>1,19=>1,20=>1,21=>1,22=>1,23=>1,24=>1,25=>1,26=>1,27=>1,28=>1,29=>1,30=>1,31=>1,32=>1,33=>1,35=>1,36=>1}}},scores=>{"HLA00001"=>{12=>{12=>1,13=>1,14=>1,15=>1,16=>1,17=>1,18=>1,19=>1,20=>1,21=>1,22=>1,23=>1,24=>1,25=>1,26=>1,27=>1,28=>1,29=>1,30=>1,31=>1,32=>1,33=>1,34=>.666666666666667,35=>1,36=>1}}}});
		my %alignments4 = ();
		$mm =  AlleleRefiner->identifyMismatches("HLA00001", 25, "TATTTCTTCACATCCGTGTCCCAGC", \%pos1_hash,\%alignments4, 1, "IIIIIIIIIIIIIIIIIIIIII?II", \%candidate_truth_hash, 1,20);	
		##Test 11 - identifyMismatches correct alignment hash	
		is_deeply(\%alignments4, \%alignment4_truth_hash, "Function identifyMismatches: 1 mismatch Q30 correct alignment hash");
		##Test 12 - identifyMismatches correct minimum mismatch
		is($mm,1,"Function identifyMismatches: 1 mismatch Q30 correct minimum mismatch");
		
		my %pos5_hash=(-6=>1);
		my %alignment5_truth_hash = (1=>{alignments=>{"HLA00001"=>{1=>{1=>'TCCCAGCA'}}},"counts"=>{"HLA00001"=>1},ref_counts=>{"HLA00001"=>{1=>{12=>1,13=>1,14=>1,15=>1,16=>1,17=>1,18=>1,10=>1,11=>1,2=>1,3=>1,4=>1,5=>1,6=>1,7=>1,8=>1,9=>1}}},scores=>{"HLA00001"=>{1=>{12=>1,13=>1,14=>1,15=>1,16=>1,17=>1,18=>1,11=>1,10=>1,9=>1,8=>1,7=>1,6=>1,5=>1,4=>1,3=>1,2=>1,1=>1}}}});
		my %alignments5 = ();
		$mm =  AlleleRefiner->identifyMismatches("HLA00001", 25, "TCCCAGCACTCCATGAGGTATTTCT", \%pos5_hash,\%alignments5, 1, "IIIIIIIIIIIIIIIIIIIIIIIII", \%candidate_truth_hash, 1,20);	
		##Test 13 - identifyMismatches correct alignment hash	
		is_deeply(\%alignments5, \%alignment5_truth_hash, "Function identifyMismatches: Beginning overhang mismatch correct alignment hash");
		##Test 14 - identifyMismatches correct minimum mismatch
		is($mm,1,"Function identifyMismatches: Beginning overhang correct minimum mismatch");
		
		%alignment5_truth_hash = (1=>{alignments=>{"HLA00001"=>{1=>{1=>'TCCCAGCA'}}},"counts"=>{"HLA00001"=>1},ref_counts=>{"HLA00001"=>{1=>{12=>1,13=>1,14=>1,15=>1,16=>1,17=>1,18=>1,10=>1,11=>1,2=>1,3=>1,4=>1,5=>1,6=>1,7=>1,8=>1,9=>1}}},scores=>{"HLA00001"=>{1=>{12=>1,13=>1,14=>1,15=>1,16=>1,17=>1,18=>1,11=>1,10=>1,9=>1,8=>1,7=>1,6=>1,5=>1,4=>1,3=>1,2=>1,1=>0}}}});
		%alignments5 = ();
		$mm =  AlleleRefiner->identifyMismatches("HLA00001", 25, "TCCCAGCACTCCATGAGGTATTTCT", \%pos5_hash,\%alignments5, 1, "II+IIIIIIIIIIIIIIIIIIIIII", \%candidate_truth_hash, 1,20);	
		##Test 15 - identifyMismatches correct alignment hash	
		is_deeply(\%alignments5, \%alignment5_truth_hash, "Function identifyMismatches: Beginning overhang Q10 mismatch correct alignment hash");
		##Test 16 - identifyMismatches correct minimum mismatch
		is($mm,1,"Function identifyMismatches: Beginning overhang Q10 correct minimum mismatch");
		
		my %pos6_hash=(60=>1);
		my %alignment6_truth_hash = (1=>{alignments=>{"HLA00001"=>{60=>{67=>'TTCCCAGCACTCCATGAG'}}},"counts"=>{"HLA00001"=>1},ref_counts=>{"HLA00001"=>{60=>{60=>1,61=>1,62=>1,63=>1,64=>1,65=>1,66=>1}}},scores=>{"HLA00001"=>{60=>{60=>1,61=>1,62=>1,63=>1,64=>1,65=>1,66=>1,67=>1}}}});
		my %alignments6 = ();
		$mm =  AlleleRefiner->identifyMismatches("HLA00001", 25, "ATCGCCGTTCCCAGCACTCCATGAG", \%pos6_hash,\%alignments6, 1, "IIIIIIIIIIIIIIIIIIIIIIIII", \%candidate_truth_hash, 1,20);	
		##Test 17 - identifyMismatches correct alignment hash	
		is_deeply(\%alignments6, \%alignment6_truth_hash, "Function identifyMismatches: Ending overhang Q10 correct alignment hash");
		##Test 18 - identifyMismatches correct minimum mismatch
		is($mm,1,"Function identifyMismatches: Ending overhang Q10 correct minimum mismatch");
		
		%alignment6_truth_hash = (1=>{alignments=>{"HLA00001"=>{60=>{67=>'TTCCCAGCACTCCATGAG'}}},"counts"=>{"HLA00001"=>1},ref_counts=>{"HLA00001"=>{60=>{60=>1,61=>1,62=>1,63=>1,64=>1,65=>1,66=>1}}},scores=>{"HLA00001"=>{60=>{60=>1,61=>1,62=>1,63=>1,64=>1,65=>1,66=>1,67=>0}}}});
		%alignments6 = ();
		$mm =  AlleleRefiner->identifyMismatches("HLA00001", 25, "ATCGCCGTTCCCAGCACTCCATGAG", \%pos6_hash,\%alignments6, 1, "IIIIIIIIIIIIIIIIIIIII+III", \%candidate_truth_hash, 1,20);	
		##Test 19 - identifyMismatches correct alignment hash	
		is_deeply(\%alignments6, \%alignment6_truth_hash, "Function identifyMismatches: Ending overhang correct alignment hash");
		##Test 20 - identifyMismatches correct minimum mismatch
		is($mm,1,"Function identifyMismatches: Ending overhang correct minimum mismatch");
		
		my @candidates  = ("HLA00001","HLA00002");
		my %mismatches1_truth_hash = ();
		my %reference_counts1_truth_hash = ("HLA0001"=>{8=>1,9=>1,10=>1,11=>1,12=>1,13=>1,14=>1,15=>1,16=>1,17=>1,18=>1,19=>1,20=>1,21=>1,22=>1,23=>1,24=>1,25=>1,26=>1,27=>1,28=>1,29=>1,30=>1,31=>1,32=>1,33=>1,42=>1,43=>1,44=>1,45=>1,46=>1,47=>1,48=>1,49=>1,50=>1,51=>1,52=>1,53=>1,54=>1,55=>1,56=>1,57=>1,58=>1,59=>1,60=>1,61=>1,62=>1,63=>1,64=>1,65=>1,66=>1,67=>1});
		my %mismatches1 = ();
		my %reference_counts1 = ();
		AlleleRefiner->alignRead("GGGACACGGATGTGAAGAAATACCTC","CGCGGGGAGCCCCGCTTCATCGCCGT",1,"IIIIIIIIIIIIIIIIIIIIIIIIII","IIIIIIIIIIIIIIIIIIIIIIIIII",20,1,20,\%reference_truth_hash,\%candidate_truth_hash,\%mismatches1,\%reference_counts1,\@candidates);
		##Test 21 - alignRead k1r_k2 correct mismatches hash 
		is_deeply(\%mismatches1,\%mismatches1_truth_hash,"Function alignRead k1r_k2 correct mismatches hash");
		##Test 22 - alignRead k1r_k2 correct mismatches hash 
		is_deeply(\%mismatches1,\%mismatches1_truth_hash,"Function alignRead k1r_k2 correct reference_counts hash");

		my %mismatches2_truth_hash = ();
		my %reference_counts2_truth_hash = ("HLA00002"=>{8=>1,9=>1,10=>1,11=>1,12=>1,13=>1,14=>1,15=>1,16=>1,17=>1,18=>1,19=>1,20=>1,21=>1,22=>1,23=>1,24=>1,25=>1,26=>1,27=>1,28=>1,29=>1,30=>1,31=>1,32=>1,33=>1,42=>1,43=>1,44=>1,45=>1,46=>1,47=>1,48=>1,49=>1,50=>1,51=>1,52=>1,53=>1,54=>1,55=>1,56=>1,57=>1,58=>1,59=>1,60=>1,61=>1,62=>1,63=>1,64=>1,65=>1,66=>1,67=>1});
		my %mismatches2 = ();
		my %reference_counts2 = ();
		AlleleRefiner->alignRead("GAGGTATTTCTCCACATCCGTGTCCC","ACTGCGATGAAGCGGGGCTCTCCACT",1,"IIIIIIIIIIIIIIIIIIIIIIIIII","IIIIIIIIIIIIIIIIIIIIIIIIII",20,1,20,\%reference_truth_hash,\%candidate_truth_hash,\%mismatches2,\%reference_counts2,\@candidates);
		##Test 23 - alignRead k1r_k2 correct mismatches hash 
		is_deeply(\%mismatches2,\%mismatches2_truth_hash,"Function alignRead k1_k2r correct mismatches hash");
		##Test 24 - alignRead k1r_k2 correct mismatches hash 
		is_deeply(\%reference_counts2,\%reference_counts2_truth_hash,"Function alignRead k1_k2r correct reference_counts hash");
		
		my %mismatches3_truth_hash = ("HLA00002"=>{10=>{'C'=>1}});
		my %reference_counts3_truth_hash = ("HLA00002"=>{8=>1,9=>1,11=>1,12=>1,13=>1,14=>1,15=>1,16=>1,17=>1,18=>1,19=>1,20=>1,21=>1,22=>1,23=>1,24=>1,25=>1,26=>1,27=>1,28=>1,29=>1,30=>1,31=>1,32=>1,33=>1,42=>1,43=>1,44=>1,45=>1,46=>1,47=>1,48=>1,49=>1,50=>1,51=>1,52=>1,53=>1,54=>1,55=>1,56=>1,57=>1,58=>1,59=>1,60=>1,61=>1,62=>1,63=>1,64=>1,65=>1,66=>1,67=>1});
		my %mismatches3 = ();
		my %reference_counts3 = ();
		AlleleRefiner->alignRead("GACGTATTTCTCCACATCCGTGTCCC","ACTGCGATGAAGCGGGGCTCTCCACT",1,"IIIIIIIIIIIIIIIIIIIIIIIIII","IIIIIIIIIIIIIIIIIIIIIIIIII",20,1,20,\%reference_truth_hash,\%candidate_truth_hash,\%mismatches3,\%reference_counts3,\@candidates);
		##Test 25 - alignRead 1 mismatch correct mismatches hash 
		is_deeply(\%mismatches3,\%mismatches3_truth_hash,"Function alignRead k1_k2r 1 mismatch correct mismatches hash");
		##Test 26 - alignRead 1 mismatch correct mismatches hash 
		is_deeply(\%reference_counts3,\%reference_counts3_truth_hash,"Function alignRead k1_k2r 1 mismatch correct reference_counts hash");
		
		my %mismatches4_truth_hash = ("HLA00002"=>{10=>{'C'=>0}});
		my %reference_counts4_truth_hash = ("HLA00002"=>{8=>1,9=>1,14=>1,11=>1,12=>1,13=>1,15=>1,16=>1,17=>1,18=>1,19=>1,20=>1,21=>1,22=>1,23=>1,24=>1,25=>1,26=>1,27=>1,28=>1,29=>1,30=>1,31=>1,32=>1,33=>1,42=>1,43=>1,44=>1,45=>1,46=>1,47=>1,48=>1,49=>1,50=>1,51=>1,52=>1,53=>1,54=>1,55=>1,56=>1,57=>1,58=>1,59=>1,60=>1,61=>1,62=>1,63=>1,64=>1,65=>1,66=>1,67=>1});
		my %mismatches4 = ();
		my %reference_counts4 = ();
		AlleleRefiner->alignRead("GACGTATTTCTCCACATCCGTGTCCC","ACTGCGATGAAGCGGGGCTCTCCACT",1,"II+IIIIIIIIIIIIIIIIIIIIIIII","IIIIIIIIIIIIIIIIIIIIIIIIII",20,1,20,\%reference_truth_hash,\%candidate_truth_hash,\%mismatches4,\%reference_counts4,\@candidates);
		##Test 27 - alignRead 1 mismatch correct mismatches hash 
		is_deeply(\%mismatches4,\%mismatches4_truth_hash,"Function alignRead k1_k2r 1 mismatch Q10 correct mismatches hash");
		##Test 28 - alignRead 1 mismatch correct mismatches hash 
		is_deeply(\%reference_counts4,\%reference_counts4_truth_hash,"Function alignRead k1_k2r 1 mismatch Q10 correct reference_counts hash");
		
		my %mismatches5_truth_hash = ("HLA00002"=>{10=>{'C'=>.666666666666667}});
		my %reference_counts5_truth_hash = ("HLA00002"=>{8=>1,9=>1,14=>1,11=>.666666666666667,12=>1,13=>1,15=>1,16=>1,17=>1,18=>1,19=>1,20=>1,21=>1,22=>1,23=>1,24=>1,25=>1,26=>1,27=>1,28=>1,29=>1,30=>1,31=>1,32=>1,33=>1,42=>1,43=>1,44=>1,45=>1,46=>1,47=>1,48=>1,49=>1,50=>1,51=>1,52=>1,53=>1,54=>1,55=>1,56=>1,57=>1,58=>1,59=>1,60=>1,61=>1,62=>1,63=>1,64=>1,65=>1,66=>1,67=>1});
		my %mismatches5 = ();
		my %reference_counts5 = ();
		AlleleRefiner->alignRead("GACGTATTTCTCCACATCCGTGTCCC","ACTGCGATGAAGCGGGGCTCTCCACT",1,"II??IIIIIIIIIIIIIIIIIIIIIII","IIIIIIIIIIIIIIIIIIIIIIIIII",20,1,20,\%reference_truth_hash,\%candidate_truth_hash,\%mismatches5,\%reference_counts5,\@candidates);
		##Test 29 - alignRead 1 mismatch correct mismatches hash 
		is_deeply(\%mismatches5,\%mismatches5_truth_hash,"Function alignRead k1_k2r 1 mismatch Q30 correct mismatches hash");
		##Test 30 - alignRead 1 mismatch correct mismatches hash 
		is_deeply(\%reference_counts5,\%reference_counts5_truth_hash,"Function alignRead k1_k2r 1 mismatch Q30 correct reference_counts hash");
		
		my %mismatches6_truth_hash = ();
		my %reference_counts6_truth_hash = ("HLA00001"=>{0=>.5,1=>.5,2=>.5,3=>.5,4=>.5,5=>.5,6=>.5,7=>.5,8=>.5,9=>.5,10=>.5,20=>.5,21=>.5,22=>.5,23=>.5,24=>.5,25=>.5,26=>.5,27=>.5,28=>.5,29=>.5,40=>.5,31=>.5,32=>.5,33=>.5,34=>.5,35=>.5,36=>.5,37=>.5,38=>.5,39=>.5,40=>.5,11=>.5,12=>.5,13=>.5,14=>.5,15=>.5,16=>.5,17=>.5,18=>.5,19=>.5},"HLA00002"=>{0=>.5,1=>.5,2=>.5,3=>.5,4=>.5,5=>.5,6=>.5,7=>.5,8=>.5,9=>.5,10=>.5,20=>.5,21=>.5,22=>.5,23=>.5,24=>.5,25=>.5,26=>.5,27=>.5,28=>.5,29=>.5,40=>.5,31=>.5,32=>.5,33=>.5,34=>.5,35=>.5,36=>.5,37=>.5,38=>.5,39=>.5,40=>.5,11=>.5,12=>.5,13=>.5,14=>.5,15=>.5,16=>.5,17=>.5,18=>.5,19=>.5});
		my %mismatches6 = ();
		my %reference_counts6 = ();
		AlleleRefiner->alignRead("GGAGAAATACCTCATGGAGT","ACATCCGTGTCCCGGCCCGG",1,"IIIIIIIIIIIIIIIIIIII","IIIIIIIIIIIIIIIIIIII",20,1,20,\%reference_truth_hash,\%candidate_truth_hash,\%mismatches6,\%reference_counts6,\@candidates);
		##Test 31 - alignRead 1 mismatch correct mismatches hash 
		is_deeply(\%mismatches6,\%mismatches6_truth_hash,"Function alignRead k1r_k2r matches both references correct mismatches hash");
		
		my %complex_input_hash = ('TCCGTGTCCCGGCCCGGCAGTGGAGAGCC'=>1,'GGCAGTGGAGAGCC'=>1,'TCCCGGCCCGGCAGTGGAGAGCC'=>1,'TGTCCCGGCCCGGCAGTGGAGAGCC'=>1,'GAGCC'=>1,'CGTGTCCCGGCCCGGCCGCGGAGCC'=>1,'GAGCC'=>1,'GCCGCGGAGCC'=>1,'TGTCCCGGCCCGGCCGCGGAGCC'=>1);
		my %complex_truth_hash1 = ('TCCGTGTCCCGGCCCGGCAGTGGAGAGCC'=>4.5, 'CGTGTCCCGGCCCGGCCGCGGAGCC'=>3.5);
		my $complex_hash1  = AlleleRefiner->compileComplexMutations(\%complex_input_hash);
		##Test 32 - compileComplexMutations correct hash
		is_deeply($complex_hash1, \%complex_truth_hash1, "Function compileComplexMutations: correct output hash");
		
		my %gene_hash = ("HLA00001"=>{seq=>"ACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGT",name=>'A*01:01:01'},"HLA00003"=>{seq=>"ACTCCATGACGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGT",name=>'A*01:01:02'},"HLA00002"=>{seq=>"ACTCCATGAGGTATTTCTCCACATCCGTGTCCCGGCCCGGCAGTGGAGAGCCCCGCTTCATCGCAGT",name=>"A*02:01"},"HLA00004"=>{seq=>"CCTACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGT",name=>'A*01:01:03'});
		my %novel_allele_mutref1=(10=>{base=>'C',freq=>20});
		my @novel_truth1 = ("HLA00003U","A*01:01:02_updated","ACTCCATGACGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGT",0); 
		my @novel_answer1 = AlleleRefiner->identifyNovelAllele("HLA00001",\%novel_allele_mutref1,\%candidate_truth_hash,\%gene_hash);
		##Test 33 - identifyNovelAlleles updated reference correct output
		is_deeply(\@novel_answer1,\@novel_truth1, "Function identifyNovelAlleles: updated reference mismatch correct output");
		
		my %novel_allele_mutref2=(1=>{base=>'TGGGCCTA',freq=>20});
		my @novel_truth2 = ("HLA00004U","A*01:01:03_updated","TGGGCCTACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGT",0); 
		my @novel_answer2 = AlleleRefiner->identifyNovelAllele("HLA00001",\%novel_allele_mutref2,\%candidate_truth_hash,\%gene_hash);
		##Test 34 - identifyNovelAlleles updated reference correct output
		is_deeply(\@novel_answer2,\@novel_truth2, "Function identifyNovelAlleles: updated reference overhang correct output");
		
		my %novel_allele_mutref3=(60=>{base=>'T',freq=>20});
		my @novel_truth3 = ("HLA00002N","A*02:01_novel","ACTCCATGAGGTATTTCTCCACATCCGTGTCCCGGCCCGGCAGTGGAGAGCCCCGCTTCTTCGCAGT",1); 
		my @novel_answer3 = AlleleRefiner->identifyNovelAllele("HLA00002",\%novel_allele_mutref3,\%candidate_truth_hash,\%gene_hash);
		##Test 35 - identifyNovelAlleles updated reference correct output
		is_deeply(\@novel_answer3,\@novel_truth3, "Function identifyNovelAlleles: novel reference mismatch correct output");
		
		my %novel_allele_mutref4=(3=>{base=>'C',freq=>"20"},67=>{base=>'TGCTTGGCAA',freq=>20});
		my @novel_truth4 = ("HLA00002N","A*02:01_novel","ACCCCATGAGGTATTTCTCCACATCCGTGTCCCGGCCCGGCAGTGGAGAGCCCCGCTTCATCGCAGTGCTTGGCAA",1); 
		my @novel_answer4 = AlleleRefiner->identifyNovelAllele("HLA00002",\%novel_allele_mutref4,\%candidate_truth_hash,\%gene_hash);
		##Test 36 - identifyNovelAlleles updated reference correct output
		is_deeply(\@novel_answer4,\@novel_truth4, "Function identifyNovelAlleles: novel reference overhang correct output");

		my %allele_hash_truth1 = ("HLA00001"=>{replace=>0},"HLA00002"=>{replace=>0});
		my %coverage_hash_truth1= ("HLA00001"=>{"mcon"=>54, "mcov"=>108, "pcon"=>4.92,"pcov"=>9.84},"HLA00002"=>{"mcon"=>302, "mcov"=>302, "pcon"=>'27.50',"pcov"=>'27.50'});
		my($allele_hash1, $coverage_hash1)= AlleleRefiner->runModule("A","$indir/inputs/HLA00001_HLA00002.A","HLA00001|HLA00002","$indir/inputs/hla.allele_refinement.fa",20,10,20,1,.75,"$FindBin::Bin/../bin/","",0,0,0);
		##Test 37 - runModule correct allele hash
		is_deeply($allele_hash1, \%allele_hash_truth1, "Fuction runModule: no mismatches correct allele hash");
		##Test 38 - runModule correct coverage_hash
		is_deeply($coverage_hash1, \%coverage_hash_truth1, "Fuction runModule: no mismatches correct coverage hash");
		
		my %allele_hash_truth2 = ("HLA00001"=>{replace=>1,novel=>1,allele=>"HLA00001N",name=>"A*01:01:01:01_novel",seq=>"ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCGGGCTCCCACTCCATGAGGTGTTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATATGAAGGCCCACTCACAGACTGACCGAGCGAACC"},"HLA00002"=>{replace=>0});
		my %coverage_hash_truth2= ("HLA00001"=>{"mcon"=>212, "mcov"=>301, "pcon"=>19.31,"pcov"=>27.41},"HLA00002"=>{"mcon"=>302, "mcov"=>302, "pcon"=>'27.50',"pcov"=>'27.50'});
		my($allele_hash2, $coverage_hash2)= AlleleRefiner->runModule("A","$indir/inputs/HLA00001N_HLA00002.A","HLA00001|HLA00002","$indir/inputs/hla.allele_refinement.fa",20,10,20,1,.75,"$FindBin::Bin/../bin/","",0,0,0);
		##Test 39 - runModule correct allele hash
		is_deeply($allele_hash2, \%allele_hash_truth2, "Fuction runModule: novel allele correct allele hash");
		##Test 40 - runModule correct allele hash
		is_deeply($coverage_hash2, \%coverage_hash_truth2, "Fuction runModul:e novel allele correct coverage hash");
		
		my @r = trap{AlleleRefiner->runModule("A","$indir/inputs/HLA00001N_HLA00002.A","","$indir/inputs/hla.allele_refinement.fa",20,10,20,1,.75,"$FindBin::Bin/../bin/","",0,0,0)};
		##Test 41 - runModule no candidates
		is($trap->exit, 1, "Function runModule no candidates exit");
		##Test 42 -  runModule no candidates error message
		is($trap->stderr, "Candidate list is empty!\n", "Function runModule: no candidates error message");

		my $out = `perl $FindBin::Bin/../bin/modules/AlleleRefiner.pm -gene A -r $indir/inputs/hla.allele_refinement.fa -c "HLA00001|HLA00002" -m 1 -pc -pm -pk -o $outdir/HLA00001_HLA00002 -f .75  -sd $FindBin::Bin/../bin/ 2>&1`;
		##Test 43 - Commandline no input
		is($?, 256, "Commandline AlleleRefiner: no input");

		$out = `perl $FindBin::Bin/../bin/modules/AlleleRefiner.pm -gene A -input $indir/inputs/HLA00001_HLA00002.A -r $indir/inputs/hla.allele_refinement.fa -c "HLA00001|HLA00002" -m 1 -pc -pm -pk -f .75  -sd $FindBin::Bin/../bin/ 2>&1`;
		##Test 44 - Commandline no output
		is($?, 256, "Commandline AlleleRefiner: no output");

		$out = `perl $FindBin::Bin/../bin/modules/AlleleRefiner.pm -input $indir/inputs/HLA00001_HLA00002.A -r $indir/inputs/hla.allele_refinement.fa -c "HLA00001|HLA00002" -m 1 -pc -pm -pk -o $outdir/HLA00001_HLA00002 -f .75  -sd $FindBin::Bin/../bin/ 2>&1`;
		##Test 45 - Commandline no gene
		is($?, 256, "Commandline AlleleRefiner: no gene");

		$out = `perl $FindBin::Bin/../bin/modules/AlleleRefiner.pm -gene A -input $indir/inputs/HLA00001_HLA00002.A -c "HLA00001|HLA00002" -m 1 -pc -pm -pk -o $outdir/HLA00001_HLA00002 -f .75  -sd $FindBin::Bin/../bin/ 2>&1`;
		##Test 46 - Commandline no reference fa
		is($?, 256, "Commandline AlleleRefiner: no reference fa");

		$out = `perl $FindBin::Bin/../bin/modules/AlleleRefiner.pm -gene A -input $indir/inputs/HLA00001_HLA00002.A -r $indir/inputs/hla.allele_refinement.fa -m 1 -pc -pm -pk -o $outdir/HLA00001_HLA00002 -f .75  -sd $FindBin::Bin/../bin/ 2>&1`;
		##Test 47 - Commandline no candidates
		is($?, 256, "Commandline AlleleRefiner: no candidates");

		$out = `perl $FindBin::Bin/../bin/modules/AlleleRefiner.pm -gene A -input $indir/inputs/HLA00001_HLA00002.A -r $indir/inputs/hla.allele_refinement.fa -c "HLA00001|HLA00002" -m 1 -pc -pm -pk -o $outdir/HLA00001_HLA00002 -f .75  -sd $FindBin::Bin/../bin/ 2>&1`;
		
		##Test 48 - Commandline no mismatch correct mismatch file
		compareFiles("$indir/outputs/HLA00001_HLA00002.mismatches.txt", "$outdir/HLA00001_HLA00002.mismatches.txt", "Commandline AlleleRefiner: no mismatch correct mismatch file");
		##Test 49 - Commandline no mismatch correct coverage file
		compareFiles("$indir/outputs/HLA00001_HLA00002.coverage.txt", "$outdir/HLA00001_HLA00002.coverage.txt", "Commandline AlleleRefiner: no mismatch correct coverage file");
		##Test 50 - Commandline no mismatch correct kmer file
		compareFiles("$indir/outputs/HLA00001_HLA00002.extra_kmers.txt", "$outdir/HLA00001_HLA00002.extra_kmers.txt", "Commandline AlleleRefiner: no mismatch correct kmer file");
		##Test 51 - Commandline no mismatch correct output
		is($out, "HLA00001\t54\t4.92\t108\t9.84\nHLA00002\t302\t27.50\t302\t27.50\n", "Commandline AlleleRefiner: no mismatches correct output");

		my $out_truth = "HLA00001\t212\t19.31\t301\t27.41\n".
				"HLA00002\t302\t27.50\t302\t27.50\n". 
				"Based on the observed data HLA00001 is really predicted to be HLA00001N A*01:01:01:01_novel.\n".
                                ">HLA00001N A*01:01:01:01_novel 304 bp\n".
				"ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCGGGCTCCCACTCCATGAGGTGTTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATATGAAGGCCCACTCACAGACTGACCGAGCGAACC\n";

		$out = `perl $FindBin::Bin/../bin/modules/AlleleRefiner.pm -gene A -input $indir/inputs/HLA00001N_HLA00002.A -r $indir/inputs/hla.allele_refinement.fa -c "HLA00001|HLA00002" -m 1 -pc -pm -pk -o $outdir/HLA00001N_HLA00002 -f .75  -sd $FindBin::Bin/../bin/ 2>&1`;
		
		##Test 52 - Commandline 1 mismatch correct mismatch file
		compareFiles("$indir/outputs/HLA00001N_HLA00002.mismatches.txt", "$outdir/HLA00001N_HLA00002.mismatches.txt", "Commandline AlleleRefiner: 1 mismatch correct mismatch file");
		##Test 53 - Commandline 1 mismatch correct coverage file
		compareFiles("$indir/outputs/HLA00001N_HLA00002.coverage.txt", "$outdir/HLA00001N_HLA00002.coverage.txt", "Commandline AlleleRefiner: 1 mismatch correct coverage file");
		##Test 54 - Commandline 1 mismatch correct kmer file
		compareFiles("$indir/outputs/HLA00001N_HLA00002.extra_kmers.txt", "$outdir/HLA00001N_HLA00002.extra_kmers.txt", "Commandline AlleleRefiner: 1 mismatch correct kmer file");
		##Test 55 - Commandline 1 mismatch correct output
		is($out, $out_truth, "Commandline AlleleRefiner: 1 mismatch correct output");

		$out = `perl $FindBin::Bin/../bin/modules/AlleleRefiner.pm -gene A -input $indir/inputs/HLA00001_HLA00002.A -r $indir/inputs/hla.allele_refinement.fa -c "HLA00001|HLA00002" -m 1 -o $outdir/HLA00001_HLA00002.nf -f .75  -sd $FindBin::Bin/../bin/ 2>&1`;
		
		##Test 56 - Commandline no output files 
		is(-e "$outdir/HLA00001_HLA00002.nf.extra_kmers.txt" || -e "$outdir/HLA00001_HLA00002.nf.coverage.txt" || -e "$outdir/HLA00001_HLA00002.nf.mismatches.txt" , undef, "Commandline AlleleRefiner: no output files");
		
	};
}

sub testHLAPredict{
	subtest 'HLAPredict.pm' =>sub{
		plan tests=>27;
		
		if(-e "$outdir/refinement_files/" && -d "$outdir/refinement_files/"){
			`rm -r $outdir/refinement_files/`;
		}
		if(-e "$outdir/function1" && -d "$outdir/predict_function1"){
			`rm -r $outdir/predict_function1/*`;
		}else{
			`mkdir -p $outdir/predict_function1/`;
		}	
		if(-e "$outdir/function2" && -d "$outdir/predict_function2"){
			`rm -r $outdir/predict_function2/*`;
		}else{
			`mkdir -p $outdir/predict_function2/`;
		}	
		if(-e "$outdir/predict_commandline" && -d "$outdir/predict_commandline"){
			`rm -r $outdir/predict_commandline/*`;
		}else{
			`mkdir -p $outdir/predict_commandline/`;
		}	
		if(-e "$outdir/A.refinement.log"){
			`rm $outdir/A.refinement.log`
		}

		#Test 1 - load module
		use_ok("HLAPredict");

		SKIP:{
			skip "DetermineProfile already tested", 1 unless (!defined $loaded_modules{"DetermineProfile"});
			#Test 2 - Test Determine Profile
			testDetermineProfile();
			$loaded_modules{"DetermineProfile"}=1;
		}
		
		my $truth_profile = DetermineProfile->determineProfile("A", "$indir/outputs/predict_profile_alleles", "sim");
		$truth_profile->{sig}{cnt}={"ATCGCAGTGGGCTACGTGGACGACACGCAG"=>2,"TGAAGGCCCACTCACAGACTGACCGAGCGA"=>6,"CATCGCAGTGGGCTACGTGGACGACACGCA"=>2,"GCCCCGCTTCATCGCAGTGGGCTACGTGGA"=>4,"CGCGAGCCAGAAGATGGAGCCGCGGGCGCC"=>6,"CTACGTGGACGACACGCAGTTCGTGCGGTT"=>6,"TGACCCAGACCTGGGCGGGCTCCCACTCCA"=>2,"ATTGGGACCAGGAGACACGGAATATGAAGG"=>8,"GGCGCCGTGGATAGAGCAGGAGGGGCCGGA"=>14,"CCGGCCCGGCAGTGGAGAGCCCCGCTTCAT"=>4,"CATCCGTGTCCCGGCCCGGCAGTGGAGAGC"=>2,"CCGCGGGCGCCGTGGATAGAGCAGGAGGGG"=>14,"AGACCTGGGCGGGCTCCCACTCCATGAGGT"=>2,"AGTGGAGAGCCCCGCTTCATCGCAGTGGGC"=>2,"GCCGCGAGCCAGAAGATGGAGCCGCGGGCG"=>6,"GCGACGCCGCGAGCCAGAAGATGGAGCCGC"=>6,"CGTGTCCCGGCCCGGCAGTGGAGAGCCCCG"=>2,"GGCCGGAGTATTGGGACCAGGAGACACGGA"=>12,"CAGTGGGCTACGTGGACGACACGCAGTTCG"=>2,"CCGCTTCATCGCAGTGGGCTACGTGGACGA"=>2,"TCATCGCAGTGGGCTACGTGGACGACACGC"=>2,"GCAGTTCGTGCGGTTCGACAGCGACGCCGC"=>6,"TGGACGACACGCAGTTCGTGCGGTTCGACA"=>6,"GCGGGCTCCCACTCCATGAGGTATTTCTCC"=>2,"TGGGCGGGCTCCCACTCCATGAGGTATTTC"=>2,"GCCCTGACCCAGACCTGGGCGGGCTCCCAC"=>2,"GGAGTATTGGGACCAGGAGACACGGAATAT"=>10,"CCGTGGATAGAGCAGGAGGGGCCGGAGTAT"=>8,"CGCCGCGAGCCAGAAGATGGAGCCGCGGGC"=>6,"GCTACGTGGACGACACGCAGTTCGTGCGGT"=>2,"AGTATTGGGACCAGGAGACACGGAATATGA"=>8,"AGAAGATGGAGCCGCGGGCGCCGTGGATAG"=>14,"CGACGCCGCGAGCCAGAAGATGGAGCCGCG"=>6,"GGGCGGGCTCCCACTCCATGAGGTATTTCT"=>2,"CCTGGGCGGGCTCCCACTCCATGAGGTATT"=>2,"AGATGGAGCCGCGGGCGCCGTGGATAGAGC"=>14,"CTTCATCGCAGTGGGCTACGTGGACGACAC"=>2,"AGACACGGAATATGAAGGCCCACTCACAGA"=>6,"GCAGTGGGCTACGTGGACGACACGCAGTTC"=>2,"TCGTGCGGTTCGACAGCGACGCCGCGAGCC"=>6,"AGCAGGAGGGGCCGGAGTATTGGGACCAGG"=>16,"AATATGAAGGCCCACTCACAGACTGACCGA"=>6,"CAGTTCGTGCGGTTCGACAGCGACGCCGCG"=>6,"GCCGGAGTATTGGGACCAGGAGACACGGAA"=>12,"AGTGGGCTACGTGGACGACACGCAGTTCGT"=>2,"CGCAGTGGGCTACGTGGACGACACGCAGTT"=>2,"GTCCCGGCCCGGCAGTGGAGAGCCCCGCTT"=>2,"CGGGCGCCGTGGATAGAGCAGGAGGGGCCG"=>14,"GCGCCGTGGATAGAGCAGGAGGGGCCGGAG"=>14,"ACCTGGGCGGGCTCCCACTCCATGAGGTAT"=>2,"GGAGACACGGAATATGAAGGCCCACTCACA"=>6,"GCCCGGCAGTGGAGAGCCCCGCTTCATCGC"=>4,"CGGCAGTGGAGAGCCCCGCTTCATCGCAGT"=>4,"GACACGGAATATGAAGGCCCACTCACAGAC"=>6,"GCGAGCCAGAAGATGGAGCCGCGGGCGCCG"=>6,"CCCTGACCCAGACCTGGGCGGGCTCCCACT"=>2,"GATGGAGCCGCGGGCGCCGTGGATAGAGCA"=>14,"CAGGAGACACGGAATATGAAGGCCCACTCA"=>6,"GCCGTGGATAGAGCAGGAGGGGCCGGAGTA"=>8,"TGGGACCAGGAGACACGGAATATGAAGGCC"=>2,"GGAATATGAAGGCCCACTCACAGACTGACC"=>6,"CGCCGTGGATAGAGCAGGAGGGGCCGGAGT"=>14,"GAGACACGGAATATGAAGGCCCACTCACAG"=>6,"GTGGAGAGCCCCGCTTCATCGCAGTGGGCT"=>2,"CGAGCCAGAAGATGGAGCCGCGGGCGCCGT"=>6,"ATATGAAGGCCCACTCACAGACTGACCGAG"=>6,"CCAGGAGACACGGAATATGAAGGCCCACTC"=>6,"AAGGCCCACTCACAGACTGACCGAGCGAAC"=>6,"GACCCAGACCTGGGCGGGCTCCCACTCCAT"=>2,"CCCCGCTTCATCGCAGTGGGCTACGTGGAC"=>4,"ACACGCAGTTCGTGCGGTTCGACAGCGACG"=>6,"CAGGAGGGGCCGGAGTATTGGGACCAGGAG"=>10,"TGGGCTACGTGGACGACACGCAGTTCGTGC"=>2,"GGAGGGGCCGGAGTATTGGGACCAGGAGAC"=>12,"GACACGCAGTTCGTGCGGTTCGACAGCGAC"=>6,"GTGGACGACACGCAGTTCGTGCGGTTCGAC"=>6,"CGCAGTTCGTGCGGTTCGACAGCGACGCCG"=>6,"CACGCAGTTCGTGCGGTTCGACAGCGACGC"=>6,"GCGGGCGCCGTGGATAGAGCAGGAGGGGCC"=>14,"CGCTTCATCGCAGTGGGCTACGTGGACGAC"=>2,"GTATTGGGACCAGGAGACACGGAATATGAA"=>8,"GGGCTACGTGGACGACACGCAGTTCGTGCG"=>2,
					  "AGGGGCCGGAGTATTGGGACCAGGAGACAC"=>12,"GGGGCCGGAGTATTGGGACCAGGAGACACG"=>12,"GAATATGAAGGCCCACTCACAGACTGACCG"=>6,"CACGGAATATGAAGGCCCACTCACAGACTG"=>6,"GGATAGAGCAGGAGGGGCCGGAGTATTGGG"=>10,"GAGCCCCGCTTCATCGCAGTGGGCTACGTG"=>2,"AGCCAGAAGATGGAGCCGCGGGCGCCGTGG"=>6,"GCAGGAGGGGCCGGAGTATTGGGACCAGGA"=>10,"CGTGGACGACACGCAGTTCGTGCGGTTCGA"=>6,"ATAGAGCAGGAGGGGCCGGAGTATTGGGAC"=>10,"GAGCAGGAGGGGCCGGAGTATTGGGACCAG"=>16,"ACGTGGACGACACGCAGTTCGTGCGGTTCG"=>6,"ACCAGGAGACACGGAATATGAAGGCCCACT"=>8,"CAGAAGATGGAGCCGCGGGCGCCGTGGATA"=>12,"AGAGCAGGAGGGGCCGGAGTATTGGGACCA"=>16,"TTCATCGCAGTGGGCTACGTGGACGACACG"=>2,"CCAGACCTGGGCGGGCTCCCACTCCATGAG"=>2,"GGGCGCCGTGGATAGAGCAGGAGGGGCCGG"=>14,"GAAGATGGAGCCGCGGGCGCCGTGGATAGA"=>14,"CCGCGAGCCAGAAGATGGAGCCGCGGGCGC"=>6,"ACGGAATATGAAGGCCCACTCACAGACTGA"=>6,"ATGAAGGCCCACTCACAGACTGACCGAGCG"=>6,"GACCAGGAGACACGGAATATGAAGGCCCAC"=>2,"GGAGCCGCGGGCGCCGTGGATAGAGCAGGA"=>8,"CAGACCTGGGCGGGCTCCCACTCCATGAGG"=>2,"AGAGCCCCGCTTCATCGCAGTGGGCTACGT"=>2,"GGCCCGGCAGTGGAGAGCCCCGCTTCATCG"=>4,"GTTCGTGCGGTTCGACAGCGACGCCGCGAG"=>6,"CGGAATATGAAGGCCCACTCACAGACTGAC"=>6,"GAGAGCCCCGCTTCATCGCAGTGGGCTACG"=>2,"TCCGTGTCCCGGCCCGGCAGTGGAGAGCCC"=>2,"CGGCCCGGCAGTGGAGAGCCCCGCTTCATC"=>4,"CGTGGATAGAGCAGGAGGGGCCGGAGTATT"=>8,"GCCAGAAGATGGAGCCGCGGGCGCCGTGGA"=>6,"CCGGCAGTGGAGAGCCCCGCTTCATCGCAG"=>4,"CCCGGCCCGGCAGTGGAGAGCCCCGCTTCA"=>2,"GGACGACACGCAGTTCGTGCGGTTCGACAG"=>6,"AGGAGACACGGAATATGAAGGCCCACTCAC"=>6,"TTCGTGCGGTTCGACAGCGACGCCGCGAGC"=>6,"CGACACGCAGTTCGTGCGGTTCGACAGCGA"=>6,"ACCCAGACCTGGGCGGGCTCCCACTCCATG"=>2,"TTGGGACCAGGAGACACGGAATATGAAGGC"=>2,"GGCGGGCTCCCACTCCATGAGGTATTTCTC"=>2,"GCAGTGGAGAGCCCCGCTTCATCGCAGTGG"=>4,"AAGATGGAGCCGCGGGCGCCGTGGATAGAG"=>14,"AGCCGCGGGCGCCGTGGATAGAGCAGGAGG"=>8,"CGGGCTCCCACTCCATGAGGTATTTCTCCA"=>2,"AGCCCCGCTTCATCGCAGTGGGCTACGTGG"=>2,"ACGACACGCAGTTCGTGCGGTTCGACAGCG"=>6,"GGCTACGTGGACGACACGCAGTTCGTGCGG"=>2,"GACGCCGCGAGCCAGAAGATGGAGCCGCGG"=>6,"GACGACACGCAGTTCGTGCGGTTCGACAGC"=>6,"GTGGGCTACGTGGACGACACGCAGTTCGTG"=>2,"GAGGGGCCGGAGTATTGGGACCAGGAGACA"=>12,"TAGAGCAGGAGGGGCCGGAGTATTGGGACC"=>10,"ACGCCGCGAGCCAGAAGATGGAGCCGCGGG"=>6,"GGAGAGCCCCGCTTCATCGCAGTGGGCTAC"=>2,"CAGTGGAGAGCCCCGCTTCATCGCAGTGGG"=>2,"CCCGGCAGTGGAGAGCCCCGCTTCATCGCA"=>4,"GATAGAGCAGGAGGGGCCGGAGTATTGGGA"=>10,"CACATCCGTGTCCCGGCCCGGCAGTGGAGA"=>2,"TCGCAGTGGGCTACGTGGACGACACGCAGT"=>2,"CGCGGGCGCCGTGGATAGAGCAGGAGGGGC"=>14,"ACATCCGTGTCCCGGCCCGGCAGTGGAGAG"=>2,"TGTCCCGGCCCGGCAGTGGAGAGCCCCGCT"=>2,"GAAGGCCCACTCACAGACTGACCGAGCGAA"=>6,"GAGCCGCGGGCGCCGTGGATAGAGCAGGAG"=>8,"ATCCGTGTCCCGGCCCGGCAGTGGAGAGCC"=>2,"GCCGCGGGCGCCGTGGATAGAGCAGGAGGG"=>14,"AGTTCGTGCGGTTCGACAGCGACGCCGCGA"=>6,"CCCAGACCTGGGCGGGCTCCCACTCCATGA"=>2,"CCGTGTCCCGGCCCGGCAGTGGAGAGCCCC"=>2,"GGGACCAGGAGACACGGAATATGAAGGCCC"=>2,"GGACCAGGAGACACGGAATATGAAGGCCCA"=>2,"GTGGATAGAGCAGGAGGGGCCGGAGTATTG"=>8,"CCCGCTTCATCGCAGTGGGCTACGTGGACG"=>4,"ATGGAGCCGCGGGCGCCGTGGATAGAGCAG"=>14,"TCCCGGCCCGGCAGTGGAGAGCCCCGCTTC"=>2,"ACACGGAATATGAAGGCCCACTCACAGACT"=>6,"TATTGGGACCAGGAGACACGGAATATGAAG"=>8,"TGGAGCCGCGGGCGCCGTGGATAGAGCAGG"=>8,"TGGATAGAGCAGGAGGGGCCGGAGTATTGG"=>10,
					  "GGCAGTGGAGAGCCCCGCTTCATCGCAGTG"=>4,"GACCTGGGCGGGCTCCCACTCCATGAGGTA"=>2,"CGGAGTATTGGGACCAGGAGACACGGAATA"=>10,"GGCCCTGACCCAGACCTGGGCGGGCTCCCA"=>2,"GTGTCCCGGCCCGGCAGTGGAGAGCCCCGC"=>2,"TGGAGAGCCCCGCTTCATCGCAGTGGGCTA"=>2,"GAGTATTGGGACCAGGAGACACGGAATATG"=>8,"CCTGACCCAGACCTGGGCGGGCTCCCACTC"=>2,"CTGGGCGGGCTCCCACTCCATGAGGTATTT"=>2,"AGGAGGGGCCGGAGTATTGGGACCAGGAGA"=>10,"TATGAAGGCCCACTCACAGACTGACCGAGC"=>6,"CCAGAAGATGGAGCCGCGGGCGCCGTGGAT"=>12,"GGGCCGGAGTATTGGGACCAGGAGACACGG"=>12,"GAGCCAGAAGATGGAGCCGCGGGCGCCGTG"=>6,"GCTTCATCGCAGTGGGCTACGTGGACGACA"=>2,"CCGGAGTATTGGGACCAGGAGACACGGAAT"=>12,"ACGCAGTTCGTGCGGTTCGACAGCGACGCC"=>6,"CTGACCCAGACCTGGGCGGGCTCCCACTCC"=>2,"TACGTGGACGACACGCAGTTCGTGCGGTTC"=>6};
	
		my $profile = DetermineProfile->determineProfile("A", "$indir/outputs/predict_profile_alleles", "sim");
			
		my ($cnts,$total_cnts) = HLAPredict->readUnique12("$indir/inputs/HLA00001_HLA00002.A", 30, $profile->{kmers}, $profile->{sig});
		#Test 3 - Function readUnique12 profile hash
		is_deeply($profile, $truth_profile, "Function readUnique12: correct profile hash");
		#Test 4 - Function readUnique12 read count
		is($cnts, 26, "Function readUnique12: correct read count");
		#Test 5 - Function readUnique12 read count
		is($total_cnts, 1092, "Function readUnique12: correct total kmer count");
		
		my @topAlleles = ("HLA00001","HLA00001N", "HLA00002");
		my %truth_combined_kmers = ("GCGGGGAGCCCCGCTTCATCGCCGTGGGCT"=>7,"ATCGCAGTGGGCTACGTGGACGACACGCAG"=>"","TGAAGGCCCACTCACAGACTGACCGAGCGA"=>"","CATCGCAGTGGGCTACGTGGACGACACGCA"=>"","CGCGAGCCAGAAGATGGAGCCGCGGGCGCC"=>"","GCCCCGCTTCATCGCAGTGGGCTACGTGGA"=>107,"TCCCGGCCCGGCCGCGGGGAGCCCCGCTTC"=>79,"CTACGTGGACGACACGCAGTTCGTGCGGTT"=>"","TGACCCAGACCTGGGCGGGCTCCCACTCCA"=>"","ATTGGGACCAGGAGACACGGAATATGAAGG"=>2,"GGCGCCGTGGATAGAGCAGGAGGGGCCGGA"=>36,"CCGGCCCGGCAGTGGAGAGCCCCGCTTCAT"=>106,"CATCCGTGTCCCGGCCCGGCAGTGGAGAGC"=>112,"CCGCGGGCGCCGTGGATAGAGCAGGAGGGG"=>"","TTCACATCCGTGTCCCGGCCCGGCCGCGGG"=>85,"AGACCTGGGCGGGCTCCCACTCCATGAGGT"=>"","AGTGGAGAGCCCCGCTTCATCGCAGTGGGC"=>"","GCCGCGAGCCAGAAGATGGAGCCGCGGGCG"=>"","TCGCCGTGGGCTACGTGGACGACACGCAGT"=>82,"GCGACGCCGCGAGCCAGAAGATGGAGCCGC"=>26,"GGGGCCCTGGCCCTGACCCAGACCTGGGCG"=>75,"CGTGTCCCGGCCCGGCAGTGGAGAGCCCCG"=>"","CCGTGGGCTACGTGGACGACACGCAGTTCG"=>16,"GGCCGGAGTATTGGGACCAGGAGACACGGA"=>37,"CAGTGGGCTACGTGGACGACACGCAGTTCG"=>96,"CCGCTTCATCGCAGTGGGCTACGTGGACGA"=>103,"GGTTCGACAGCGACGCCGCGAGCCAGAAGA"=>39,"GCAGTTCGTGCGGTTCGACAGCGACGCCGC"=>"","TCATCGCAGTGGGCTACGTGGACGACACGC"=>100,"GTATTTCTTCACATCCGTGTCCCGGCCCGG"=>77,"TGGACGACACGCAGTTCGTGCGGTTCGACA"=>49,"GCGGGCTCCCACTCCATGAGGTATTTCTCC"=>114,"TGGGCGGGCTCCCACTCCATGAGGTATTTC"=>"","ACTCCATGAGGTATTTCTTCACATCCGTGT"=>51,"GCCCTGACCCAGACCTGGGCGGGCTCCCAC"=>"","GGAGTATTGGGACCAGGAGACACGGAATAT"=>10,"CCGTGGATAGAGCAGGAGGGGCCGGAGTAT"=>9,"CGCCGCGAGCCAGAAGATGGAGCCGCGGGC"=>24,"GCTACGTGGACGACACGCAGTTCGTGCGGT"=>"","AGTATTGGGACCAGGAGACACGGAATATGA"=>"","AGAAGATGGAGCCGCGGGCGCCGTGGATAG"=>"","GCGGGCTCCCACTCCATGAGGTATTTCTTC"=>67,"CGACGCCGCGAGCCAGAAGATGGAGCCGCG"=>"","GGGCGGGCTCCCACTCCATGAGGTATTTCT"=>72,"CCTGGGCGGGCTCCCACTCCATGAGGTATT"=>"","AGATGGAGCCGCGGGCGCCGTGGATAGAGC"=>"","CTTCATCGCAGTGGGCTACGTGGACGACAC"=>113,"CCGCTTCATCGCCGTGGGCTACGTGGACGA"=>38,"AGACACGGAATATGAAGGCCCACTCACAGA"=>"","CTGGGCGGGCTCCCACTCCATGAGGTGTTT"=>88,"GCAGTGGGCTACGTGGACGACACGCAGTTC"=>"","CCCGCTTCATCGCCGTGGGCTACGTGGACG"=>19,"TCGTGCGGTTCGACAGCGACGCCGCGAGCC"=>4,"AGCAGGAGGGGCCGGAGTATTGGGACCAGG"=>1,"AATATGAAGGCCCACTCACAGACTGACCGA"=>"","CAGTTCGTGCGGTTCGACAGCGACGCCGCG"=>"","GCCGGAGTATTGGGACCAGGAGACACGGAA"=>"","GGTATTTCTTCACATCCGTGTCCCGGCCCG"=>76,"AGTGGGCTACGTGGACGACACGCAGTTCGT"=>"","CGCAGTGGGCTACGTGGACGACACGCAGTT"=>"","GTCCCGGCCCGGCAGTGGAGAGCCCCGCTT"=>"","TGAGGTATTTCTTCACATCCGTGTCCCGGC"=>48,"GGCCGTCATGGCGCCCCGAACCCTCCTCCT"=>69,"GGGCTCCCACTCCATGAGGTATTTCTCCAC"=>116,"GCGCCGTGGATAGAGCAGGAGGGGCCGGAG"=>"","CGGGCGCCGTGGATAGAGCAGGAGGGGCCG"=>"","GGTGTTTCTTCACATCCGTGTCCCGGCCCG"=>92,"CCGAACCCTCCTCCTGCTACTCTCGGGGGC"=>55,"ACCTGGGCGGGCTCCCACTCCATGAGGTAT"=>"","GGAGACACGGAATATGAAGGCCCACTCACA"=>"","AGGCCCACTCACAGACTGACCGAGCGAACC"=>31,"CACATCCGTGTCCCGGCCCGGCCGCGGGGA"=>44,"GCCCGGCAGTGGAGAGCCCCGCTTCATCGC"=>"","CGGCAGTGGAGAGCCCCGCTTCATCGCAGT"=>"","GACACGGAATATGAAGGCCCACTCACAGAC"=>"","GCGAGCCAGAAGATGGAGCCGCGGGCGCCG"=>"","CCCTGACCCAGACCTGGGCGGGCTCCCACT"=>"","GATGGAGCCGCGGGCGCCGTGGATAGAGCA"=>"","CAGGAGACACGGAATATGAAGGCCCACTCA"=>"","GCCGTGGATAGAGCAGGAGGGGCCGGAGTA"=>"",
			"TGGGACCAGGAGACACGGAATATGAAGGCC"=>"","GGAATATGAAGGCCCACTCACAGACTGACC"=>"","GTATTTCTCCACATCCGTGTCCCGGCCCGG"=>118,"CGCCGTGGATAGAGCAGGAGGGGCCGGAGT"=>0,"GAGACACGGAATATGAAGGCCCACTCACAG"=>22,"CGAGCCAGAAGATGGAGCCGCGGGCGCCGT"=>"","GTGGAGAGCCCCGCTTCATCGCAGTGGGCT"=>95,"ATATGAAGGCCCACTCACAGACTGACCGAG"=>"","CCAGGAGACACGGAATATGAAGGCCCACTC"=>"","AAGGCCCACTCACAGACTGACCGAGCGAAC"=>5,"GACCCAGACCTGGGCGGGCTCCCACTCCAT"=>64,"CCCCGCTTCATCGCAGTGGGCTACGTGGAC"=>101,"TCCTGCTACTCTCGGGGGCCCTGGCCCTGA"=>81,"ACACGCAGTTCGTGCGGTTCGACAGCGACG"=>"","CAGGAGGGGCCGGAGTATTGGGACCAGGAG"=>"","TGGGCTACGTGGACGACACGCAGTTCGTGC"=>"","TCCACATCCGTGTCCCGGCCCGGCAGTGGA"=>119,"GGAGGGGCCGGAGTATTGGGACCAGGAGAC"=>"","CCCGAACCCTCCTCCTGCTACTCTCGGGGG"=>54,"GACACGCAGTTCGTGCGGTTCGACAGCGAC"=>"","GTGGACGACACGCAGTTCGTGCGGTTCGAC"=>"","CGCAGTTCGTGCGGTTCGACAGCGACGCCG"=>"","TCCTCCTGCTACTCTCGGGGGCCCTGGCCC"=>80,"GCGCCCCGAACCCTCCTCCTGCTACTCTCG"=>66,"GCGGGCGCCGTGGATAGAGCAGGAGGGGCC"=>"","CTTCATCGCCGTGGGCTACGTGGACGACAC"=>63,"CACGCAGTTCGTGCGGTTCGACAGCGACGC"=>27,"GGGGAGCCCCGCTTCATCGCCGTGGGCTAC"=>74,"ACTCCATGAGGTATTTCTCCACATCCGTGT"=>109,"CGCTTCATCGCAGTGGGCTACGTGGACGAC"=>102,"GTCATGGCGCCCCGAACCCTCCTCCTGCTA"=>78,"CCCCGCTTCATCGCCGTGGGCTACGTGGAC"=>32,"GTATTGGGACCAGGAGACACGGAATATGAA"=>"","GGGCTACGTGGACGACACGCAGTTCGTGCG"=>"","AGGGGCCGGAGTATTGGGACCAGGAGACAC"=>"","CATCCGTGTCCCGGCCCGGCCGCGGGGAGC"=>53,"GGGGCCGGAGTATTGGGACCAGGAGACACG"=>"","GAATATGAAGGCCCACTCACAGACTGACCG"=>"","CACGGAATATGAAGGCCCACTCACAGACTG"=>"","GGATAGAGCAGGAGGGGCCGGAGTATTGGG"=>"","TCATCGCCGTGGGCTACGTGGACGACACGC"=>28,"CGCTTCATCGCCGTGGGCTACGTGGACGAC"=>33,"GAGCCCCGCTTCATCGCAGTGGGCTACGTG"=>"","AGCCAGAAGATGGAGCCGCGGGCGCCGTGG"=>"","CGTGGACGACACGCAGTTCGTGCGGTTCGA"=>"","GCAGGAGGGGCCGGAGTATTGGGACCAGGA"=>"","ATAGAGCAGGAGGGGCCGGAGTATTGGGAC"=>"","GAGCAGGAGGGGCCGGAGTATTGGGACCAG"=>"","ACGTGGACGACACGCAGTTCGTGCGGTTCG"=>50,"ACCAGGAGACACGGAATATGAAGGCCCACT"=>8,"CAGAAGATGGAGCCGCGGGCGCCGTGGATA"=>34,"AGAGCAGGAGGGGCCGGAGTATTGGGACCA"=>41,"TTCATCGCAGTGGGCTACGTGGACGACACG"=>99,"CCAGACCTGGGCGGGCTCCCACTCCATGAG"=>"","GGGCGCCGTGGATAGAGCAGGAGGGGCCGG"=>"","GACCAGGAGACACGGAATATGAAGGCCCAC"=>"","ACGGAATATGAAGGCCCACTCACAGACTGA"=>"","CCGCGAGCCAGAAGATGGAGCCGCGGGCGC"=>25,"ATGAAGGCCCACTCACAGACTGACCGAGCG"=>17,"GAAGATGGAGCCGCGGGCGCCGTGGATAGA"=>"","GGAGCCGCGGGCGCCGTGGATAGAGCAGGA"=>"","CAGACCTGGGCGGGCTCCCACTCCATGAGG"=>"","AGAGCCCCGCTTCATCGCAGTGGGCTACGT"=>94,"GGCCCGGCAGTGGAGAGCCCCGCTTCATCG"=>"","CGGAATATGAAGGCCCACTCACAGACTGAC"=>"","GTTCGTGCGGTTCGACAGCGACGCCGCGAG"=>"","GGGCTCCCACTCCATGAGGTATTTCTTCAC"=>73,"GAGAGCCCCGCTTCATCGCAGTGGGCTACG"=>97,"ACTCCATGAGGTGTTTCTTCACATCCGTGT"=>87,"TCCGTGTCCCGGCCCGGCAGTGGAGAGCCC"=>"","CGGCCCGGCAGTGGAGAGCCCCGCTTCATC"=>"","GCGGGCTCCCACTCCATGAGGTGTTTCTTC"=>89,"TGAGGTGTTTCTTCACATCCGTGTCCCGGC"=>86,"GCGGTTCGACAGCGACGCCGCGAGCCAGAA"=>43,"GCCAGAAGATGGAGCCGCGGGCGCCGTGGA"=>"","CGTGGATAGAGCAGGAGGGGCCGGAGTATT"=>"","GGAGCCCCGCTTCATCGCCGTGGGCTACGT"=>6,"CCGGCAGTGGAGAGCCCCGCTTCATCGCAG"=>"","CCCGGCCCGGCAGTGGAGAGCCCCGCTTCA"=>"","GGACGACACGCAGTTCGTGCGGTTCGACAG"=>68,
			"TTCATCGCCGTGGGCTACGTGGACGACACG"=>20,"AGGAGACACGGAATATGAAGGCCCACTCAC"=>"","TTCGTGCGGTTCGACAGCGACGCCGCGAGC"=>"","CGACACGCAGTTCGTGCGGTTCGACAGCGA"=>"","ACCCAGACCTGGGCGGGCTCCCACTCCATG"=>"","TTGGGACCAGGAGACACGGAATATGAAGGC"=>"","GGCGGGCTCCCACTCCATGAGGTATTTCTC"=>"","GCAGTGGAGAGCCCCGCTTCATCGCAGTGG"=>"","GCCCCGAACCCTCCTCCTGCTACTCTCGGG"=>65,"AGCCGCGGGCGCCGTGGATAGAGCAGGAGG"=>15,"AAGATGGAGCCGCGGGCGCCGTGGATAGAG"=>"","CGGGCTCCCACTCCATGAGGTATTTCTCCA"=>"","ACGACACGCAGTTCGTGCGGTTCGACAGCG"=>"","GGTATTTCTCCACATCCGTGTCCCGGCCCG"=>117,"AGCCCCGCTTCATCGCAGTGGGCTACGTGG"=>110,"GACGCCGCGAGCCAGAAGATGGAGCCGCGG"=>14,"GGCTACGTGGACGACACGCAGTTCGTGCGG"=>13,"GACGACACGCAGTTCGTGCGGTTCGACAGC"=>29,"GAGGGGCCGGAGTATTGGGACCAGGAGACA"=>"","GTGGGCTACGTGGACGACACGCAGTTCGTG"=>12,"ACGCCGCGAGCCAGAAGATGGAGCCGCGGG"=>"","TAGAGCAGGAGGGGCCGGAGTATTGGGACC"=>"","GGAGAGCCCCGCTTCATCGCAGTGGGCTAC"=>115,"CCCGGCAGTGGAGAGCCCCGCTTCATCGCA"=>"","CAGTGGAGAGCCCCGCTTCATCGCAGTGGG"=>111,"GATAGAGCAGGAGGGGCCGGAGTATTGGGA"=>"","CACATCCGTGTCCCGGCCCGGCAGTGGAGA"=>105,"TCGCAGTGGGCTACGTGGACGACACGCAGT"=>121,"CGCGGGCGCCGTGGATAGAGCAGGAGGGGC"=>"","ACATCCGTGTCCCGGCCCGGCAGTGGAGAG"=>"","TGTCCCGGCCCGGCAGTGGAGAGCCCCGCT"=>"","GGGAGCCCCGCTTCATCGCCGTGGGCTACG"=>18,"GAAGGCCCACTCACAGACTGACCGAGCGAA"=>"","CCTCCTCCTGCTACTCTCGGGGGCCCTGGC"=>57,"ATCCGTGTCCCGGCCCGGCAGTGGAGAGCC"=>"","GAGCCGCGGGCGCCGTGGATAGAGCAGGAG"=>21,"CCGCGGGGAGCCCCGCTTCATCGCCGTGGG"=>56,"AGTTCGTGCGGTTCGACAGCGACGCCGCGA"=>"","GCCGCGGGCGCCGTGGATAGAGCAGGAGGG"=>"","CCCAGACCTGGGCGGGCTCCCACTCCATGA"=>"","CCGTGTCCCGGCCCGGCAGTGGAGAGCCCC"=>"","GGGACCAGGAGACACGGAATATGAAGGCCC"=>"","CGCCCCGAACCCTCCTCCTGCTACTCTCGG"=>58,"CGTCATGGCGCCCCGAACCCTCCTCCTGCT"=>60,"GTGTTTCTTCACATCCGTGTCCCGGCCCGG"=>93,"GGACCAGGAGACACGGAATATGAAGGCCCA"=>40,"GTGGATAGAGCAGGAGGGGCCGGAGTATTG"=>"","CCCGCTTCATCGCAGTGGGCTACGTGGACG"=>98,"TCCCGGCCCGGCAGTGGAGAGCCCCGCTTC"=>120,"ATGGAGCCGCGGGCGCCGTGGATAGAGCAG"=>3,"ACACGGAATATGAAGGCCCACTCACAGACT"=>"","TGGATAGAGCAGGAGGGGCCGGAGTATTGG"=>"","TGGAGCCGCGGGCGCCGTGGATAGAGCAGG"=>"","TATTGGGACCAGGAGACACGGAATATGAAG"=>23,"GGCAGTGGAGAGCCCCGCTTCATCGCAGTG"=>"","GACCTGGGCGGGCTCCCACTCCATGAGGTA"=>"","CGGAGTATTGGGACCAGGAGACACGGAATA"=>35,"GTGTCCCGGCCCGGCAGTGGAGAGCCCCGC"=>"","GCCCCGCTTCATCGCCGTGGGCTACGTGGA"=>46,"GGCCCTGACCCAGACCTGGGCGGGCTCCCA"=>47,"GAGTATTGGGACCAGGAGACACGGAATATG"=>"","TGGAGAGCCCCGCTTCATCGCAGTGGGCTA"=>122,"CCTGACCCAGACCTGGGCGGGCTCCCACTC"=>"","GCTTCATCGCCGTGGGCTACGTGGACGACA"=>42,"CTGGGCGGGCTCCCACTCCATGAGGTATTT"=>62,"TATGAAGGCCCACTCACAGACTGACCGAGC"=>"","AGGAGGGGCCGGAGTATTGGGACCAGGAGA"=>"","GGGCGGGCTCCCACTCCATGAGGTGTTTCT"=>90,"TCTCGGGGGCCCTGGCCCTGACCCAGACCT"=>83,"CCAGAAGATGGAGCCGCGGGCGCCGTGGAT"=>"","GGGCCGGAGTATTGGGACCAGGAGACACGG"=>"","TGGCGCCCCGAACCCTCCTCCTGCTACTCT"=>84,"GGGCCCTGGCCCTGACCCAGACCTGGGCGG"=>71,
			"GAGCCAGAAGATGGAGCCGCGGGCGCCGTG"=>"","CGGGGAGCCCCGCTTCATCGCCGTGGGCTA"=>59,"GCTTCATCGCAGTGGGCTACGTGGACGACA"=>104,"CCGGAGTATTGGGACCAGGAGACACGGAAT"=>11,"AGCCCCGCTTCATCGCCGTGGGCTACGTGG"=>52,"ACGCAGTTCGTGCGGTTCGACAGCGACGCC"=>30,"GGCGCCCCGAACCCTCCTCCTGCTACTCTC"=>70,"CCGGCCCGGCCGCGGGGAGCCCCGCTTCAT"=>45,"CTGACCCAGACCTGGGCGGGCTCCCACTCC"=>61,"TGAGGTATTTCTCCACATCCGTGTCCCGGC"=>108,"TACGTGGACGACACGCAGTTCGTGCGGTTC"=>"","GGGCTCCCACTCCATGAGGTGTTTCTTCAC"=>91);

		my $combined_kmers = HLAPredict->findCombinedKmers($truth_profile->{sig}, $truth_profile->{names},\@topAlleles, 0);
		#Test 6 - Find combined kmers
		is_deeply($combined_kmers, \%truth_combined_kmers, "Function findCombinedKmers: correct combined hash");

		my ($cor, $propReads, $propSig) = HLAPredict->getSSE("HLA00001", "HLA00002", $total_cnts, $truth_profile->{totals}{HLA00001} + $truth_profile->{totals}{HLA00002}, 1/($cnts*2), $truth_profile->{sig},2);
		#Test 7 - getSSE correct correlation
		is_deeply($cor, -0.176606658649094, "Function getSSE: correct correlation");
		#Test 8 - getSSE correct proportion of signal
		is_deeply($propReads, 0.304029304029304, "Function getSSE: correct proportion of signal");
		#Test 9 - getSSE correct proportion of reads
		is_deeply($propSig, 0.530434782608696, "Function getSSE: correct proportion of reads");

		my @pp_results = ([0,10,5],[1,0,2],[3,5,0]);
		my @finalists = ("HLA00001|HLA00002","HLA00001|HLA00001N","HLA00001N|HLA00002");
		my %truth_matrix = ("HLA00001|HLA00002"=>{row=>15.02,col=>4.02},"HLA00001|HLA00001N"=>{row=>3.02,col=>15.02},"HLA00001N|HLA00002"=>{row=>8.02,col=>7.02});
		
		my $matrix_sums = HLAPredict->sumMatrix(\@pp_results, \@finalists, .01);
		#Test 10 - sumMatrix correct hash
		is_deeply($matrix_sums,\%truth_matrix, "Function sumMatrix: correct matrix sums");

		my %options_truth = (allele_map=>"",cntThresh=>2, numFinalists=>5, numSemifinalists=>10, PairPickerPower=>1,minQuality=>20,threads=>1, pp_scale=>0, out_dir=>"$outdir",scripts_dir=>"$FindBin::Bin/../bin/", allele_refinement=>"all", kraken_db=>"$indir/outputs/database", kraken_path=>$kraken_path, sim_num_reads=>50,sim_read_length=>50, sim_max_insert=>50, sim_seed=>1234,sim_scale=>80, sim_shape=>0.7, unique1=>"",fastq_dir=>"", log=>"",minReads=>100,reference=>"",profile=>"", output=>"");
		my %sim_opts = (num_reads=>50,read_length=>50, max_insert=>50, seed=>1234,scale=>80, shape=>0.7);
		my $options = HLAPredict->setOptions(2,5,10,1,20,1,0,"$outdir","$FindBin::Bin/../bin/","all",\%sim_opts,"$indir/outputs/database",$kraken_path,100);
		#Test 11 - setOptions
		is_deeply($options, \%options_truth, "Function setOptions: correct options");


		my %add_hash =(kmers=>{"CGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGGGG"=>0},names=>{0=>"CGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGGGG"},totals=>{"HLA00000"=>10},sig=>{"HLA00000"=>{0=>10}} );

		my %add_hash_truth =(kmers=>{"CGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGGGG"=>0,"CGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCC"=>1,"CATCCAGATAATGTATGGCTGCGACGTGGGGCCGGACGGGCGCTTCCTCC"=>2,"CCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGG"=>3,"AGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCA"=>4,"ACCCTGCGCGGCTACTACAACCAGAGCGAGGACGGTTCTCACACCATCCA"=>5,"CACGGAATATGAAGGCCCACTCACAGACTGACCGAGCGAACCTGGGGACC"=>6,"GGCCGGAGTATTGGGACCAGGAGACACGGAATATGAAGGCCCACTCACAG"=>7,"CAGACCTGGGCGGGCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTC"=>8,"TGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATAT"=>9,"CGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTT"=>10},
				     names=>{0=>"CGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGGGG",1=>"CGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCC",,2=>"CATCCAGATAATGTATGGCTGCGACGTGGGGCCGGACGGGCGCTTCCTCC",3=>"CCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGG",4=>"AGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCA",5=>"ACCCTGCGCGGCTACTACAACCAGAGCGAGGACGGTTCTCACACCATCCA",6=>"CACGGAATATGAAGGCCCACTCACAGACTGACCGAGCGAACCTGGGGACC",7=>"GGCCGGAGTATTGGGACCAGGAGACACGGAATATGAAGGCCCACTCACAG",8=>"CAGACCTGGGCGGGCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTC",9=>"TGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATAT",10=>"CGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTT"},
				     totals=>{"HLA00000"=>10,"HLA00001"=>1584},sig=>{"HLA00000"=>{0=>10},"HLA00001"=>{1=>165,2=>163,3=>160,4=>159,5=>158,6=>157,7=>157,8=>156,9=>155,10=>154}} );
		my %add_combined_kmer_truth = ("CGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGGGG"=>0,"CGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCC"=>1,"CATCCAGATAATGTATGGCTGCGACGTGGGGCCGGACGGGCGCTTCCTCC"=>2,"CCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGG"=>3,"AGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCA"=>4,"ACCCTGCGCGGCTACTACAACCAGAGCGAGGACGGTTCTCACACCATCCA"=>5,"CACGGAATATGAAGGCCCACTCACAGACTGACCGAGCGAACCTGGGGACC"=>6,"GGCCGGAGTATTGGGACCAGGAGACACGGAATATGAAGGCCCACTCACAG"=>7,"CAGACCTGGGCGGGCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTC"=>8,"TGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATAT"=>9,"CGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTT"=>10);
		my %add_combined_kmers = ("CGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGGGG"=>0);
		HLAPredict->addReadsToProfile("$indir/outputs/filtered/sim.HLA00001.A_2.uniq.cnts","HLA00001",0,$add_hash{kmers},$add_hash{names},\%add_combined_kmers,$add_hash{sig},$add_hash{totals});	
		
		#Test 12 - addReadsToProfile correct profile
		is_deeply(\%add_hash, \%add_hash_truth, "Function addReadsToProfile: correct profile"); 
		#Test 13 - addReadsToProfile correct combined kmers
		is_deeply(\%add_combined_kmers, \%add_combined_kmer_truth, "Function addReadsToProfile: correct combined kmers"); 

		my $new_allele = "HLA00001N";
		my $seq1 = "ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCGGGCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATATGAAGGCCCACTCACAGACTGACCGAGCGAACC";
		my $seq2 = "ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCGGGCTCCCACTCCATGAGGTATTTCTCCACATCCGTGTCCCGGCCCGGCAGTGGAGAGCCCCGCTTCATCGCAGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATATGAAGGCCCACTCACAGACTGACCGAGCGAACC";
		my $truth_profile1 = DetermineProfile->determineProfile("A", "$indir/outputs/predict_profile_alleles", "profiletest");
		my ($cnts1,$total_cnts1) = HLAPredict->readUnique12("$indir/inputs/HLA00001_HLA00002.A", 50, $truth_profile1->{kmers}, $truth_profile1->{sig});
		my %base_profile1=(sig=>{},kmers=>{},totals=>{},names=>{});	
		my ($cnts2,$total_cnts2) = HLAPredict->readUnique12("$indir/inputs/HLA00001_HLA00002.A", 50, $base_profile1{kmers}, $base_profile1{sig});
			
		my @topAlleles1 = ("HLA00001", "HLA00002");	
		my $combined_kmers1 = HLAPredict->findCombinedKmers($base_profile1{sig}, $base_profile1{names},\@topAlleles1,0);
		my $combined_kmers_truth1 = retrieve "$indir/outputs/predict_profile_alleles/profiletest.combined_kmers.ph";		
		my ($ncor, $npropReads, $npropSig, $refined_reference)= HLAPredict->createProfileAndCalculateSignal("A","HLA00001","HLA00002","A*01:01:01:01","A*01:02",1,1,$seq1,$seq2,52,(1/52),$base_profile1{sig}, $base_profile1{names},$base_profile1{kmers},$base_profile1{totals},"",2);
		my $tmp = 0;
	
		#Test 14 - addProfileAndCalculateSignal correct profile
		is_deeply(\%base_profile1, $truth_profile1, "Function addProfileAndCalculateSignal: correct profile");
		#Test 15 - addProfileAndCalculateSignal combined kmers
		is_deeply($combined_kmers1, $combined_kmers_truth1, "Function addProfileAndCalculateSignal: correct combined kmer list");
		#Test 16 - addProfileAndCalculateSignal correct correlation
		is_deeply($ncor, 0.631216028002741, "Function addProfileAndCalculateSignal: correct correlation");
		#Test 17 - addProfileAndCalculateSignal correct proportion of signal
		is_deeply($npropReads, 1, "Function addProfileAndCalculateSignal: correct proportion of signal");
		#Test 18 - addProfileAndCalculateSignal correct proportion of reads
		is_deeply($npropSig, 0.114754098360656, "Function addProfileAndCalculateSignal: correct proportion of reads");
		
		my $answer_truth1="HLA00002\tHLA00001\tA*01:02\tA*01:01:01:01\t1\t0.306550218340611\t0.341718511967625\t1.35173126969176\t1681\t1243.59037753364\t\t\n".
			"HLA00001\tHLA00091\tA*01:01:01:01\tA*30:03\t1\t0.230314960629921\t0.274559999791732\t1.49512503957835\t1.31089225341761\t0.876777673248858\t\t\n".
			"HLA00001\tHLA00089\tA*01:01:01:01\tA*30:01:01\t1\t0.227626459143969\t0.259689152402247\t1.51268438845378\t0.903835705748126\t0.597504484509156\t\t\n".
			"HLA00002\tHLA10377\tA*01:02\tA*11:160\t0.91375\t0.213515456506111\t0.230189928340654\t1.64254461515324\t0.473639599932192\t0.288357220597022\t\t\n".
			"HLA00002\tHLA00120\tA*01:02\tA*68:04\t0.91375\t0.196039603960396\t0.202576547592082\t1.68763384844752\t0.377849683777026\t0.22389316505154\t\t\n";
		my $answer1 = HLAPredict->predictHLAType("$indir/inputs/full.HLA00001_HLA00002.A_1.uniq.cnts", "$indir/inputs/", "$indir/outputs/database/data/kmer_profiles/A.profile.ph","$indir/inputs/hla.ref.fa","$indir/outputs/hla.ref.merged.allele_map.txt","$outdir/predict.log","$outdir/predict_function1");
		#Test 19 - predictHLAType
		is($answer1,$answer_truth1,"Function predictHLAType: no novel allele correct answer");		

		my $answer_truth2 = "HLA00002\tHLA00001N\tA*01:02\tA*01:01:01:01_novel\t0.78705\t0.24952741020794\t0.299177261657758\t1.6642453281343\t1842.55\t1107.13845419989\t\t Predicted novel allele based on differences between top allele HLA00001 and observed data.\n".
			"HLA00002\tHLA00001\tA*01:02\tA*01:01:01:01\t0.92105\t0.268122270742358\t0.343920267864929\t1.46690746139271\t6.62712397734424\t4.51775190444005\t\t\n".
			"HLA00001\tHLA00091\tA*01:01:01:01\tA*30:03\t0.92105\t0.201443569553806\t0.264947373742645\t1.61255905670355\t0.548125867653864\t0.339910569709219\t\t\n".
			"HLA00002\tHLA10377\tA*01:02\tA*11:160\t0.902125\t0.213515456506111\t0.26681725879586\t1.61754228469803\t0.494031422803476\t0.305421025142291\t\t\n".
			"HLA00001\tHLA00089\tA*01:01:01:01\tA*30:01:01\t0.92105\t0.199092088197147\t0.255606619379559\t1.62425129242329\t0.391547749725576\t0.241063529733388\t\t\n".
			"HLA00002\tHLA00120\tA*01:02\tA*68:04\t0.902125\t0.196039603960396\t0.238338616264372\t1.66349677977523\t0.38202495307497\t0.229651753895543\t\t\n";
		my $answer2 = HLAPredict->predictHLAType("$indir/inputs/full.HLA00001N_HLA00002.A_1.uniq.cnts", "$indir/inputs/", "$indir/outputs/database/data/kmer_profiles/A.profile.ph","$indir/inputs/hla.ref.fa","$indir/outputs/hla.ref.merged.allele_map.txt","$outdir/predict.log","$outdir/predict_function2");
		#Test 20 - predictHLAType
		is($answer2,$answer_truth2,"Function predictHLAType: novel allele correct answer");		
	
		my $out = `perl $FindBin::Bin/../bin/modules/HLAPredict.pm -unique1 $indir/inputs/full.HLA00001N_HLA00002.A_1.uniq.cnt -fastq_dir $indir/inputs -profile $indir/outputs/database/data/kmer_profiles/A.profile.ph -allele_map$indir/ outputs/hla.ref.merged.allele_map.txt -out_dir test_dir/predict_commandline 2>&1`;
		#Test 21 - Commandline HLAPredict no reference
		is($?, 256, "Commandline HLAPredict: no reference"); 
		
		$out = `perl $FindBin::Bin/../bin/modules/HLAPredict.pm -fastq_dir $indir/inputs -reference $indir/inputs/hla.ref.fa -profile $indir/outputs/database/data/kmer_profiles/A.profile.ph -allele_map $indir/outputs/hla.ref.merged.allele_map.txt -out_dir test_dir/predict_commandline 2>&1`; 
		#Test 22 - Commandline HLAPredict no unique1
		is($?, 256, "Commandline HLAPredict: no unique1"); 
		
		$out = `perl $FindBin::Bin/../bin/modules/HLAPredict.pm -unique1 $indir/inputs/full.HLA00001N_HLA00002.A_1.uniq.cnt -reference $indir/inputs/hla.ref.fa -profile $indir/outputs/database/data/kmer_profiles/A.profile.ph -allele_map $indir/outputs/hla.ref.merged.allele_map.txt -out_dir test_dir/predict_commandline 2>&1`; 
		#Test 23 - Commandline HLAPredict no fastq directory
		is($?, 256, "Commandline HLAPredict: no fastq directory");
 
		$out = `perl $FindBin::Bin/../bin/modules/HLAPredict.pm -unique1 $indir/inputs/full.HLA00001N_HLA00002.A_1.uniq.cnt -fastq_dir $indir/inputs -reference $indir/inputs/hla.ref.fa -allele_map $indir/outputs/hla.ref.merged.allele_map.txt -out_dir test_dir/predict_commandline 2>&1`; 
		#Test 24 - Commandline HLAPredict no profile
		is($?, 256, "Commandline HLAPredict: no profile");
		
		$out = `perl $FindBin::Bin/../bin/modules/HLAPredict.pm -unique1 $indir/inputs/full.HLA00001N_HLA00002.A_1.uniq.cnt -fastq_dir $indir/inputs -reference $indir/inputs/hla.ref.fa -profile $indir/outputs/database/data/kmer_profiles/A.profile.ph  -out_dir test_dir/predict_commandline 2>&1`; 
		#Test 25 - Commandline HLAPredict no allele_map
		is($?, 256, "Commandline HLAPredict: no allele_map");
 
		$out = `perl $FindBin::Bin/../bin/modules/HLAPredict.pm -unique1 $indir/inputs/full.HLA00001N_HLA00002.A_1.uniq.cnts -fastq_dir $indir/inputs -reference $indir/inputs/hla.ref.fa -profile $indir/outputs/database/data/kmer_profiles/A.profile.ph  -out_dir test_dir/predict_commandline -allele_map "a" 2>&1`; 
		#Test 26 - Commandline HLAPredict invalid allele_map
		is($?, 256, "Commandline HLAPredict: invalid allele_map");
		
		$out = `perl $FindBin::Bin/../bin/modules/HLAPredict.pm -unique1 $indir/inputs/full.HLA00001N_HLA00002.A_1.uniq.cnts -fastq_dir $indir/inputs -reference $indir/inputs/hla.ref.fa -profile $indir/outputs/database/data/kmer_profiles/A.profile.ph -allele_map $indir/outputs/hla.ref.merged.allele_map.txt -out_dir $outdir/predict_commandline -sd $FindBin::Bin/../bin/ -kraken_path ~/opt/pkg-big/kraken_ea -kraken_db $indir/outputs/database -numFinalists 5 -numSemifinalists 10 -PairPickerPower 1 -pp_scale 0 -allele_refinement all -sim_max_insert 50 -sim_read_length 50 -sim_num_reads 50`; 
		#Test 27 - Commandline HLAPredict correct predictions
		is($out, $answer_truth2, "Commandline HLAPredict: correct predictions");
	};
}

sub testHLAProfiler{
	subtest 'HLAProfiler.pl'=>sub{
		plan tests=>8;
		testHLAProfilerBuild();
		testHLAProfilerPredict();
		testHLAProfilerCreateTaxonomy();
		testHLAProfilerBuildTaxonomy();
		testHLAProfilerCreateProfiles();
		testHLAProfilerFilter();
		testHLAProfilerCountReads();
		testHLAProfilerPredictOnly();
	};
}

sub testHLAProfilerBuild{
	subtest 'HLAProfiler.pl build' =>sub{
		plan tests=>10;
	};
}


sub testHLAProfilerPredict{
	subtest 'HLAProfiler.pl predict' =>sub{
		plan tests=>10;
	};
}


sub testHLAProfilerCreateTaxonomy{
	subtest 'HLAProfiler.pl create_taxonomy' =>sub{
		plan tests=>10;
	};
}

sub testHLAProfilerBuildTaxonomy{
	subtest 'HLAProfiler.pl build_taxonomy' =>sub{
		plan tests=>10;
	};
}


sub testHLAProfilerCreateProfiles{
	subtest 'HLAProfiler.pl create_profiles' =>sub{
		plan tests=>10;
	};
}

sub testHLAProfilerFilter{
	subtest 'HLAProfiler.pl filter' =>sub{
		plan tests=>10;
	};
}

sub testHLAProfilerCountReads{
	subtest 'HLAProfiler.pl count_reads' =>sub{
		plan tests=>10;
	};
}

sub testHLAProfilerPredictOnly{
	subtest 'HLAProfiler.pl predict_only' =>sub{
		plan tests=>10;
	};
}
__PACKAGE__->runCommandline() unless caller;
