#!/usr/bin/env perl

(my $SCRIPTS_DIR = $0) =~ s/HLAProfiler.pl//;
(my $SCRIPT_NAME = $0) =~ s/.*\///;

use warnings;
use strict;
use Parallel::ForkManager;
use Getopt::Long;
use Module::Load;
use File::Copy;

my $version = "1.0";
my $creation_date = "1 Oct 2016";
my $last_updated = "13 Jan 2017";

my $usage = "\n$SCRIPT_NAME v$version\n" .
	    "\nDESCRIPTION\n" .
	    "A tool for predicting HLA types using NGS Paired-end sequencing data.\n" .
	    "\nUSAGE:\n" . 
	    "perl $SCRIPT_NAME <module>\n" .
	    "\nGlobal Modules. These module will run the process from start to finish.\n" .
	    "  build\t\t\tBuild the HLAProfiler reference using a reference FASTA\n" .
	    "  predict\t\tPredict the HLA types using paired end sequencing data.\n" .
	    "  test_modules\t\tThis modules simply tests whether HLAProfiler.pl can access the required accessory perl modules\n" .
	    "\nBuild Module Components. These module complete individual steps in the build module.\n" .
	    "  create_taxonomy\tCreates a taxonomy needed to create a HLA database\n" .
	    "  build_taxonomy\tBuilds the taxonomy into a database\n" .
	    "  create_profiles\tSimulated reads from the HLA reference and creates kmer profiles\n" . 
	    "\nPredict Module Components. These module complete individual steps in the predict module.\n" .	    
	    "  filter\t\tFilters paired fastq files into HLA genes using and HLA database\n" .
	    "  count_reads\t\tCounts the reads in fastq files to create a read count file\n" .
	    "  predict_only\t\tOnly predicts HLA type using previously filtered and counted reads\n" .  
	    "\nFor module specific instructions run each individual module with the help flag.\n" .
	    "i.e. perl $SCRIPT_NAME  build -h\n" .
	            "\nAUTHORS:\n" . 
	    "Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	    "Chad Brown\n" .
	    "\nCREATED:\n$creation_date\n" .
	    "\nLAST UPDATED:\n$last_updated\n" .
	    "\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
	    "\n";

my $module = shift || "";
my $kraken_path =`which kraken`;
$kraken_path=~s/scripts\/kraken//;
chomp($kraken_path);

my $log;
if($module eq "test_modules"){
	load "$SCRIPTS_DIR/modules/RunKraken.pm";
	load "$SCRIPTS_DIR/modules/MergeDuplicates.pm";
	load "$SCRIPTS_DIR/modules/HLATaxonomy.pm";
	load "$SCRIPTS_DIR/modules/TaxonomyDivisions.pm";
	load "$SCRIPTS_DIR/modules/RunKraken.pm";
	load "$SCRIPTS_DIR/modules/HLADistractome.pm";
	load "$SCRIPTS_DIR/modules/SimulateReads.pm";
	load "$SCRIPTS_DIR/modules/DetermineProfile.pm";
	load "$SCRIPTS_DIR/modules/RunKraken.pm";
	load "$SCRIPTS_DIR/modules/ReadCounter.pm";
	load "$SCRIPTS_DIR/modules/HLAPredict.pm";
	print "Modules loaded successfully\n";	
}elsif($module eq "build"){
	my $build_usage = "\n$SCRIPT_NAME build\n" .
			"\nDESCRIPTION\n" .
			"A tool for building the HLAProfiler reference from a fasta file containing the HLA reference and a fasta file containing GENCODE transcripts.\n" .
			"\nUSAGE:\n" .
			"perl $SCRIPT_NAME build <options>\n" .
			"\nRequired Options:\n" .
			"-transcripts|t\t\tlocation of fasta file containing transcripts. Currently only GENCODE transcripts are supported.(required)\n" .
			"-transcript_gtf|g\tlocation of gtf file containing transcripts corresponding to the -transcripts option. Currently only GENCODE transcripts are supported.(required)\n" .
			"-exclusion_bed|e\tlocation of bed file containing the coordinated any regions to be excluded from the distractome. i.e. HLA region.(required)\n" .
			"-reference|r\t\tlocation of fasta file containing HLA reference. IPD-IMGT/HLA reference recommended.(required)\n" .
			"-cwd\t\t\tFile containing the names of common and well-documented alleles. This file can be blank but must be specified.(required)\n" .
			"\nOutput Options:\n" .
			"-output_dir|o\t\tlocation of output directory(default:\".\")\n" .
			"-database_name|db\tname of the HLA database to be created(default:hla)\n" .
			"-kraken_path|kp\tbase directory of kraken installation. (default:base directory of path returned by `which kraken`)\n" . 
			"\nHLA database creation\n" .
			"-k_mer|k\t\tsize of the k-mer used to create database.(default:31)\n" .
			"-minimizer|mi\t\tsize of the k-mer minimizer used to crate database.(default:13)\n" .
			"\nK-mer profile creation\n" .
			"-num_reads|nr\t\tnumber of reads to simulated per reference allele for k-mer profile creations.(default:500000)\n" .
			"-read_length|rl\t\tlength of reads simulated for k-mer profile. Same as the length of the k-mers in the profile.(default:50)\n" .
			"-max_insert|m\t\tmaximum size of insert (default:1000)\n" .
			"-scale|sc\t\tscale of pareto distribution to determine insert size (default:80)\n" .
			"-shape|sh\t\tshape of pareto distribution to determine insert size (default:0.7)\n" .
			"-seed|sd\t\tseed of random number generator for simulation (default:1234)\n" .
			"-filter_reads|f\t\ttoggle whether or not to filter reads using in the HLA database when building the k-mer profile.It is STRONGLY recommended to use the default for this setting. Possibile values 0 or 1. (default:1)\n" .
			"-intermediate_files|if\ttoggles flag to keep intermediate files (default:off)\n" .
			"\nGeneral options:\n" .
			"-threads|c\tnumber of threads to uses for processing.(default:1)\n" .
			"-help|h\t\tprints this help prompt\n" .
	            		"\nAUTHORS:\n" . 
	    		"Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	    		"Chad Brown:chad.brown\@q2labsolutions.com\n" .
	    		"\nCREATED:\n$creation_date\n" .
	    		"\nLAST UPDATED:\n$last_updated\n" .
	    		"\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
	    		"\n";
			  
	my %build_opts = (transcripts=>"", transcript_gtf=>"", reference=>"", output_dir=>".",cwd=>"", database_name=>"hla", num_reads=>500000,read_length=>50,exclusion_bed=>"", threads=>1, filter_reads=>1,k_mer=>31, minimizer=>13, kraken_path=>$kraken_path, max_insert=>1000, seed=>1234, scale=>80, shape=>0.7);
	GetOptions(\%build_opts, qw(reference|r=s transcript_gtf|g=s exclusion_bed|e=s output_dir|o=s cwd=s database_name|db=s num_reads|nr=s intermediate_files|if read_length|l=s transcripts|t=s threads|c=s filter_reads|f=s help|h minimizer|m=s k_mer|k=s kraken_path|kp=s max_insert|mi=s scale|sc=s shape|sh=s seed|sd=s));
	if($build_opts{help}){
		print "$build_usage\n";
		exit;
	}elsif($build_opts{reference} eq ""){
		print "Missing required option -reference. Please enter the location of the reference fasta.\n$build_usage\n";
		exit;
	}elsif($build_opts{cwd} eq ""){
		print "Missing required option -cwd. Please enter the location of the cwd file.\n$build_usage\n";
		exit;
	}elsif (! -e $build_opts{reference}){
		print "File $build_opts{reference} does not exist. Please enter a valid location for the reference fasta.\n$build_usage\n";
		exit;
	}elsif ($build_opts{transcripts} eq ""){
		print "Missing required option -transcripts. Please enter a valid location for the transcripts fasta.\n$build_usage\n";
		exit;
	}elsif (! -e $build_opts{transcripts}){
		print "File $build_opts{transcripts} does not exist. Please enter a valid location for the transcripts fasta.\n$build_usage\n";
		exit;
	}elsif ($build_opts{transcript_gtf} eq ""){
		print "Missing required option -transcript_gtf. Please enter a valid location for the transcripts gtf.\n$build_usage\n";
		exit;
	}elsif (! -e $build_opts{transcript_gtf}){
		print "File $build_opts{transcript_gtf} does not exist. Please enter a valid location for the transcripts gtf.\n$build_usage\n";
		exit;
	}elsif ($build_opts{exclusion_bed} eq ""){
		print "Missing required option -exclusion_bed. Please enter a valid location for the exclusion bed.\n$build_usage\n";
		exit;
	}elsif (! -e $build_opts{exclusion_bed}){
		print "File $build_opts{exclusion_bed} does not exist. Please enter a valid location for the exclusion bed.\n$build_usage\n";
		exit;
	}else{
		run_build(\%build_opts);
	}
}elsif($module eq "predict"){
	my $predict_usage = "\n$SCRIPT_NAME predict\n" .
			"\nDESCRIPTION\n" .
			"A tool for predicting the HLA type of Paired-end NGS data.\n" .
			"\nUSAGE:\n" .
			"perl $SCRIPT_NAME predict <options>\n" .
			"\nRequired Options:\n" .
			"-fastq1|fq1\t\tlocation of read1 fastq (required)\n" .
			"-fastq2|fq2\t\tlocation of read2 fastq (required)\n" .
			"-database_name|db\tname of HLA database (required)\n" .
			"-directory_dir|dd\tname of parent directory of database (required)\n" .
			"-reference|r\treference fa used to create the database (required)\n" .
			"\nAllele Refinement Options:\n" .
			"-allele_refinement|ar\tSpecifies the level to which the predicted alleles are to be refined based on the observed reads (default:all)\n" .
			"     Possible values:\n" .
			"\trefine_only\tRefines the allelle call by looking predicting the true allele sequence using observed reads and looking for a better match in the reference\n" .
			"\tpredict_only\tReports if the observe reads support a novel allele sequence not found in the reference\n" .
			"\trefineAndPredict\tRefines the allele call (-refine_only) and report novel alleles (-novel_only)\n" . 
			"\tall\t\tRefines the allele call (-refine_only) and report novel alleles (-novel_only), creates a profile for the refined/novel allele sequence and calculates prediction metrics.\n" . 
			"\tnone\t\tTurns off refinement and novel allele prediction.\n" . 
			"\n-num_reads|nr\t\tnumber of reads to simulated per reference allele for k-mer profile creations.(default:500000)\n" .
			"-read_length|rl\t\tlength of reads simulated for k-mer profile. Same as the length of the k-mers in the profile.(default:50)\n" .
			"-max_insert|m\t\tmaximum size of insert (default:1000)\n" .
			"-scale|sc\t\tscale of pareto distribution to determine insert size (default:80)\n" .
			"-shape|sh\t\tshape of pareto distribution to determine insert size (default:0.7)\n" .
			"-seed|sd\t\tseed of random number generator for simulation (default:1234)\n" .
			"\nGeneral Options:\n" .
			"-intermediate_files|if\ttoggles flag to keep intermediate files (default:off)\n" .
			"-minimum_reads|min\tminimum number of reads from a gene before attempting to call HLA types.(default:100)\n" .
			"-threads|c\t\tnumber of threads (default:1)\n" .
			"-output_dir|od\t\toutput directory (default:\" .\")\n" .
			"-kraken_path|kp\t\tbase directory of kraken installation. (default:base directory of path returned by `which kraken`)\n" . 
			"-log|l\t\t\tname of the prediction log file\n" .
			"-help|h\t\t\tprints this help prompt\n" .
	            		"\nAUTHORS:\n" . 
	    		"Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	    		"Chad Brown:chad.brown\@q2labsolutions.com\n" .
	    		"\nCREATED:\n$creation_date\n" .
	    		"\nLAST UPDATED:\n$last_updated\n" .
	    		"\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
	    		"\n";
	my %predict_opts=(kraken_path=>$kraken_path,database_dir=>"", fastq1=>"", fastq2=>"",database_name=>"",output_dir=>".", threads=>1, log=>"", allele_refinement=>"all", num_reads=>500000,read_length=>50, max_insert=>1000, seed=>1234, scale=>80, shape=>0.7, minimum_reads=>100);
	GetOptions(\%predict_opts, qw(allele_refinement|ar=s intermediate_files|if fastq1|fq1=s fastq2|fq2=s threads|c=s reference|r=s database_name|db=s database_dir|dd=s num_reads|nr=s read_length|rl=s max_insert|m=s scale|sc=s shape|sh=s seed|sd=s output_dir|od=s kraken_path|kp=s log|l=s threads|c=s help|h minimum_reads|min=s));
	if($predict_opts{help}){
		print "$predict_usage\n";
		exit;
	}elsif(! (-e $predict_opts{fastq1} && -e $predict_opts{fastq2})){#
		print "Please enter valid paired-end fastq files.\n$predict_usage\n";
		exit; 
	}elsif($predict_opts{database_dir} eq ""){
		print "Please include a database_directory using -database_dir.\n$predict_usage\n";
		exit;
	}elsif($predict_opts{database_name} eq ""){
		print "Please include a database_name using -database_name.\n$predict_usage\n";
		exit;
	}elsif(! -e $predict_opts{reference}){
		print "Please include a valid reference fasta using -reference.\n$predict_usage\n";
		exit;
	}else{
		my $allele_map = $predict_opts{reference};
		$allele_map =~ s/\.fa/\.allele_map.txt/;
		if(! -e $allele_map){
			print "Could not find the allele map $allele_map in the reference directory.\n";
			exit;
		}
		run_predict(\%predict_opts, $allele_map);
	}
}elsif($module eq "create_taxonomy"){
	my $create_usage = "\n$SCRIPT_NAME create_taxonomy\n" .
			"\nDESCRIPTION\n" .
			"A tool for creating a kraken-compatable taxonomy for HLA alleles from a reference fasta.\n" .
			"\nUSAGE:\n" .
			"perl $SCRIPT_NAME create_taxonomy <options>\n" .
			"\nRequired Options\n" .
			"-reference|r\tHLA reference fasta (required)\n" .
			"-cwd\t\tFile containing the names of common and well-documented alleles. This file can be blank but must be specified.(required)\n" .
			"\nOptions:\n" .
			"-output_dir|od\tparent directory of taxonomy (default:\".\")\n" .  
			"-help|h\ti\tprints this help prompt\n" .
			"-log|l\t\tname of the prediction log file\n" .
	            		"\nAUTHORS:\n" . 
	    		"Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	    		"Chad Brown:chad.brown\@q2labsolutions.com\n" .
	    		"\nCREATED:\n$creation_date\n" .
	    		"\nLAST UPDATED:\n$last_updated\n" .
	    		"\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
			"\n";
			
	my %create_opts=(output_dir=>".", reference=>"", cwd=>"");
	GetOptions(\%create_opts, qw(output_dir|od=s reference|r=s cwd=s help|h log|l));
	if($create_opts{help}){
		print "$create_usage\n";
		exit;
	}elsif(! (-e $create_opts{reference})){
		print "Please enter a valid reference sequence.\n$create_usage\n";
		exit; 
	}elsif(! (-e $create_opts{cwd})){
		print "Please specify a cwd file.\n$create_usage\n";
		exit; 
	}else{
		my $logfile = "$create_opts{output_dir}/HLAProfiler.build.log";
		open (my $log, ">$logfile}");	
		open(my $commands, ">$create_opts{output_dir}/HLAProfiler.build.commands.txt");
		my $reference = mergeDuplicateAlleles($create_opts{reference}, $create_opts{output_dir},$create_opts{cwd}, $log, $commands, $logfile);
		createHLATaxonomy($reference,"$create_opts{output_dir}",$log, $commands, $logfile);
	}
	
}elsif($module eq "build_taxonomy"){
	my $build_usage = "\n$SCRIPT_NAME build_taxonomy\n" .
			"\nDESCRIPTION\n" .
			"A tool for building an HLA database using a reference and custom taxonomy\n" .
			"\nUSAGE: $SCRIPT_NAME build_taxonomy <options>\n" .
			"\nRequired Options:\n" .
			"-transcripts|t\t\tlocation of fasta file containing transcripts. Currently only GENCODE transcripts are supported.(required)\n" .
			"-transcript_gtf|g\tlocation of gtf file containing transcripts corresponding to the -transcripts option. Currently only GENCODE transcripts are supported.(required)\n" .
			"-exclusion_bed|e\tlocation of bed file containing the coordinated any regions to be excluded from the distractome. i.e. HLA region.(required)\n" .
			"-reference|r\t\tlocation of fasta file containing HLA reference. IPD-IMGT/HLA reference recommended.(required)\n" .
			"\nOutput Options:\n" .
			"-output_dir|o\t\tlocation of database directory(default:\".\")\n" .
			"-database_name|db\tname of the HLA database to be created(default:hla)\n" .
			"-kraken_path|kp\t\tbase directory of kraken installation. (default:base directory of path returned by `which kraken`)\n" . 
			"\nHLA database creation\n" .
			"-k_mer|k\t\tsize of the k-mer used to create database.(default:31)\n" .
			"-minimizer|m\t\tsize of the k-mer minimizer used to crate database.(default:13)\n" .
			"\nGeneral options:\n" .
			"-threads|c\t\tnumber of threads to uses for processing.(default:1)\n" .
			"-help|h\t\t\tprints this help prompt\n" .
	            		"\nAUTHORS:\n" . 
	    		"Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	    		"Chad Brown:chad.brown\@q2labsolutions.com\n" .
	    		"\nCREATED:\n$creation_date\n" .
	    		"\nLAST UPDATED:\n$last_updated\n" .
	    		"\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
			"\n"; 
	
	my %build_opts = (exclusion_bed=>"", transcripts=>"",transcript_gtf=>"", reference=>"", output_dir=>".",log=>"", database_name=>"hla", threads=>1, k_mer=>31, minimizer=>13, kraken_path=>$kraken_path);
	GetOptions(\%build_opts, qw(reference|r=s transcript_gtf|g=s exclusion_bed|e=s output_dir|o=s database_name|db=s transcripts|t=s threads|c=s help|h minimizer|m=s k_mer|k=s kraken_path|kp=s));
	if($build_opts{help}){
		print "$build_usage\n";
		exit;
	}elsif($build_opts{reference} eq ""){
		print "Missing required option -reference. Please enter the location of the reference fasta.\n$build_usage\n";
		exit;
	}elsif (! -e $build_opts{reference}){
		print "File $build_opts{reference} does not exist. Please enter a valid location for the reference fasta.\n$build_usage\n";
		exit;
	}elsif ($build_opts{transcripts} eq ""){
		print "Missing required option -transcripts. Please enter a valid location for the transcripts fasta.\n$build_usage\n";
		exit;
	
	}elsif (! -e $build_opts{transcripts}){
		print "File $build_opts{transcripts} does not exist. Please enter a valid location for the transcripts fasta.\n$build_usage\n";
		exit;
	}elsif ($build_opts{transcript_gtf} eq ""){
		print "Missing required option -transcript_gtf. Please enter a valid location for the transcripts gtf.\n$build_usage\n";
		exit;
	}elsif (! -e $build_opts{transcript_gtf}){
		print "File $build_opts{transcript_gtf} does not exist. Please enter a valid location for the transcripts gtf.\n$build_usage\n";
		exit;
	}elsif ($build_opts{exclusion_bed} eq ""){
		print "Missing required option -exclusion_bed. Please enter a valid location for the exclusion bed.\n$build_usage\n";
		exit;
	}elsif (! -e $build_opts{exclusion_bed}){
		print "File $build_opts{exclusion_bed} does not exist. Please enter a valid location for the exclusion bed.\n$build_usage\n";
		exit;
	}else{
		my $mkdir_cmds = "";
		if(! (-e "$build_opts{output_dir}/$build_opts{database_name}/data" && -d "$build_opts{output_dir}/$build_opts{database_name}/data")){
			if(! mkdir "$build_opts{output_dir}/$build_opts{database_name}/data"){
				print STDERR "Fatal Error: Error creating directory $build_opts{output_dir}/$build_opts{database_name}\n";
				exit; 	
			}else{
				$mkdir_cmds .= "mkdir $build_opts{output_dir}/$build_opts{database_name}/data\n"
			}
		}
		if(! (-e "$build_opts{output_dir}/$build_opts{database_name}/data/reference" && -d "$build_opts{output_dir}/$build_opts{database_name}/data/reference")){
			if(! mkdir "$build_opts{output_dir}/$build_opts{database_name}/data/reference"){
				print STDERR "Fatal Error: Error creating directory $build_opts{output_dir}/$build_opts{database_name}\n";
				exit; 	
			}else{
				$mkdir_cmds .= "mkdir $build_opts{output_dir}/$build_opts{database_name}/data/reference\n"
			}
		}
		if(! (-e "$build_opts{output_dir}/$build_opts{database_name}/data/logs" && -d "$build_opts{output_dir}/$build_opts{database_name}/data/logs")){
			if(! mkdir "$build_opts{output_dir}/$build_opts{database_name}/data/logs"){
				print STDERR "Fatal Error: Error creating directory $build_opts{output_dir}/$build_opts{database_name}/data/logs\n";
				exit; 	
			}else{
				$mkdir_cmds .= "mkdir $build_opts{output_dir}/$build_opts{database_name}/data/logs\n"
			}
		}
		my $logfile = "$build_opts{output_dir}/$build_opts{database_name}/data/logs/HLAProfiler.build.log";
		open (my $log, ">$logfile}");	
		open(my $commands, ">$build_opts{output_dir}/$build_opts{database_name}/data/logs/HLAProfiler.build.commands.txt");
		print $commands "$mkdir_cmds";
		buildHLADatabase($build_opts{reference},$build_opts{transcripts},$build_opts{transcript_gtf},$build_opts{exclusion_bed},$build_opts{output_dir},$build_opts{database_name},$build_opts{k_mer},$build_opts{minimizer},$build_opts{threads},$build_opts{kraken_path}, $log,$commands);
}
}elsif($module eq "create_profiles"){
	my $build_usage = "\n$SCRIPT_NAME create_profiles\n" .
			"\nDESCRIPTION\n" .
			"A tool for creating a k-mer profile to used with HLAProfile.pl predict\n" .
			"\nUSAGE:\n" .
			"perl $SCRIPT_NAME create_profile <options>\n" .
			"\nRequired Options:\n" .
			"-reference|r\tlocation of HLA reference fasta file.(required)\n" .
			"\nOutput Options:\n" .
			"-output_dir|o\t\tlocation of output directory(default:\".\")\n" .
			"-database_dir|dd\tlocation of database parent directory(default:\".\")\n" .
			"-database_name|db\tname of the HLA database to be created(default:hla)\n" .
			"-kraken_path|kp\t\tbase directory of kraken installation. (default:base directory of path returned by `which kraken`)\n" . 
			"\nK-mer profile creation\n" .
			"-num_reads|nr\t\tnumber of reads to simulated per reference allele for k-mer profile creations.(default:500000)\n" .
			"-read_length|rl\t\tlength of reads simulated for k-mer profile. Same as the length of the k-mers in the profile.(default:50)\n" .
			"-filter_reads|f\t\ttoggle whether or not to filter reads using in the HLA database when building the k-mer profile.It is STRONGLY recommended to use the default for this setting. Possibile values 0 or 1. (default:1)\n" .
			"-intermediate_files|if\ttoggles flag to keep intermediate files (default:off)\n" .
			"-max_insert|mi\t\tmaximum size of insert (default:1000)\n" .
			"-scale|sc\t\tscale of pareto distribution to determine insert size (default:80)\n" .
			"-shape|sh\t\tshape of pareto distribution to determine insert size (default:0.7)\n" .
			"-seed\t\t\tseed of random number generator for simulation (default:1234)\n" .
			"\nGeneral options:\n" .
			"-threads|c\t\tnumber of threads to uses for processing.(default:1)\n" .
			"-help|h\t\t\tprints this help prompt\n" .
	            		"\nAUTHORS:\n" . 
	    		"Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	    		"Chad Brown:chad.brown\@q2labsolutions.com\n" .
	    		"\nCREATED:\n$creation_date\n" .
	    		"\nLAST UPDATED:\n$last_updated\n" .
	    		"\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
			"\n";  
	my %build_opts = (reference=>"", output_dir=>".", database_dir=>".", database_name=>"hla", num_reads=>500000,read_length=>50, threads=>1, filter_reads=>1,k_mer=>31, minimizer=>13, kraken_path=>$kraken_path, max_insert=>1000, seed=>1234, scale=>80, shape=>0.7);
	GetOptions(\%build_opts, qw(intermediate_files|if reference|r=s output_dir|o=s database_name|db=s num_reads|nr=s read_length|l=s threads|c=s filter_reads|f=s help|h k_mer|k=s kraken_path|kp=s seed|sd=s scale|sc=s shape|sh=s max_insert|mi=s));
	if($build_opts{help}){
		print "$build_usage\n";
		exit;
	}elsif($build_opts{reference} eq ""){
		print "Missing required option -reference. Please enter the location of the reference fasta.\n$build_usage\n";
		exit;
	}elsif (! -e $build_opts{reference}){
		print "File $build_opts{reference} does not exist. Please enter a valid location for the reference fasta.\n$build_usage\n";
	}else{
		my $mkdir_cmds = "";
		if(! (-e "$build_opts{output_dir}/logs" && -d "$build_opts{output_dir}/logs")){
			if(! mkdir "$build_opts{output_dir}/logs"){
				print STDERR "Fatal Error: Error creating directory $build_opts{output_dir}/logs\n";
				exit; 	
			}else{
				$mkdir_cmds .= "mkdir $build_opts{output_dir}/logs\n"
			}
		}
		my $logfile = "$build_opts{output_dir}/logs/HLAProfiler.build.log";
		open (my $log, ">$logfile}");	
		open(my $commands, ">$build_opts{output_dir}/logs/HLAProfiler.build.commands.txt");
		print $commands "$mkdir_cmds";
		createKmerProfiles($build_opts{reference},$build_opts{database_dir},$build_opts{database_name},$build_opts{output_dir},$build_opts{num_reads},$build_opts{read_length},$build_opts{filter_reads}, $build_opts{threads}, $build_opts{kraken_path}, $build_opts{intermediate_files},$build_opts{max_insert}, $build_opts{scale}, $build_opts{shape}, $build_opts{seed}, $log, $commands);
	}
}elsif($module eq "filter"){
	my $filter_usage = "\n$SCRIPT_NAME filter\n" .
			"\nDESCRIPTION\n" .
			"A tool for filter paired fastq reads using a HLA database\n" .  
			"\nUSAGE:\n" .
			"perl $SCRIPT_NAME filter <options>\n" .
			"\nRequired Options:\n" .
			"-fastq1|fq1\t\tread1 fastq.(required)\n" .
			"-fastq2|fq2\t\tread2 fastq.(required)\n" .
			"\nOutput Options:\n" .
			"-output_dir|od\t\tlocation of output directory. (default:\".\")\n" .
			"-database_dir|dd\tlocation of database directory(default:\".\")\n" .
			"-database_name|db\tname of the HLA database to be created(default:hla)\n" .
			"-kraken_path|kp\t\tbase directory of kraken installation. (default:base directory of path returned by `which kraken`)\n". 
			"\nGeneral options:\n" .
			"-threads|c\t\tnumber of threads to uses for processing.(default:1)\n" .
			"-help|h\t\t\tprints this help prompt\n" .
	            		"\nAUTHORS:\n" . 
	    		"Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
	    		"Chad Brown:chad.brown\@q2labsolutions.com\n" .
	    		"\nCREATED:\n$creation_date\n" .
	    		"\nLAST UPDATED:\n$last_updated\n" .
	    		"\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
			"\n";  
	my %filter_opts = (database_dir=>"", output_dir=>"", database_name=>"hla", threads=>1, kraken_path=>$kraken_path, fastq1=>"", fastq2=>"");
	GetOptions(\%filter_opts, qw(database_dir|dd=s output_dir|od=s database_name|db=s threads|c=s fastq1|fq1=s fastq2|fq2=s kraken_path|kp=s help|h));
	if($filter_opts{help}){
		print "$filter_usage\n";
		exit;
	}elsif(! (-e $filter_opts{fastq1} && $filter_opts{fastq2})){
		print "Please enter valid paired-end fastq files.\n$filter_usage\n";
		exit; 
	}else{
		my $fastq_prefix = $filter_opts{fastq1};
		$fastq_prefix =~ s/.*\///;
		$fastq_prefix =~ s/_1\.f(ast)*q(\.gz)*//;
		filterReads($filter_opts{database_dir}, $filter_opts{database_name}, "$filter_opts{output_dir}/filtered",$fastq_prefix, $filter_opts{threads},$filter_opts{fastq1}, $filter_opts{fastq2},$filter_opts{kraken_path},"$filter_opts{output_dir}/HLAProfiler.filter.log");	
	}
}elsif($module eq "count_reads"){
	my $count_usage = "\n$SCRIPT_NAME count_reads\n" .
			"\nDESCRIPTION\n" .
			"A tool for counting reads in fastq files for a sample in a directory.\n" .  
			"\nUSAGE:\n" .
			"perl $SCRIPT_NAME count_reads <options>\n" .
			"\nRequired Options:\n" .
			"-reads_directory\tlocation of directory containing filtered read fastqs. Please make sure to filter files using HLAProfiler.pl filter before counting (required)\n" .
			"-sample_name|sn\t\tname of the sample. This must perfect match the prefix of each of the read count files. i.e. The sample name for file NA12878.200.B_1.uniq.cnt would be NA12878.200 (required)\n" .
			"-output_directory|od\tlocation of directory containing filtered read fastqs. Please make sure to filter files using HLAProfiler.pl filter before counting (default:-reads_directory)\n" .
			"\nGeneral options:\n" .
			"-threads|c\t\tnumber of threads to uses for processing.(default:1)\n" .
			"-help|h\t\t\tprints this help prompt\n" .
				"\nAUTHORS:\n" . 
			"Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
			"Chad Brown:chad.brown\@q2labsolutions.com\n" .
			"\nCREATED:\n$creation_date\n" .
			"\nLAST UPDATED:\n$last_updated\n" .
			"\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
			"\n";  
	my %count_opts = (threads=>1, reads_directory=>"", output_directory=>"",sample_name=>"");
	GetOptions(\%count_opts, qw(reads_directory|rd=s output_directory|od=s sample_name|sn=s threads|c=s help|h));

	if($count_opts{help}){
		print "$count_usage\n";
		exit;
	}elsif(! (-e $count_opts{reads_directory} && -d $count_opts{reads_directory})){
		print "Please enter valid directory containing filtered reads.\n$count_usage\n";
		exit; 
	}elsif($count_opts{output_directory} eq ""){
		$count_opts{output_directory} = $count_opts{reads_directory};
	}
	else{
		countUniqueReads($count_opts{reads_directory},$count_opts{output_directory},$count_opts{sample_name},$count_opts{threads},"$count_opts{reads_directory}/HLAProfiler.count.log");
	}

}elsif($module eq "predict_only"){
	my $predict_usage = "\n$SCRIPT_NAME predict_only\n" .
			"\nDESCRIPTION\n" .
			"A tool for predicting the HLA type of Paired-end NGS data from paired fastqs or read count files.\n" .
			"\nUSAGE:\n" .
			"perl $SCRIPT_NAME predict <options>\n" .
			"\nRequired Options:\n" .
			"-counts_directory|cd\t\tlocation of directory containing filtered and paired read counts files. To generate these files from fastq files please run HLAProfiler.pl filter followed by HLAProfiler.pl count_reads (required)\n" .
			"-reads_directory|cd\t\tlocation of directory containing filtered and paired read fastqs.(required)\n" .
			"-profile_directory|sdir\t\tpath to directory containing the profile files (required)\n" .
			"-sample_name|sn\t\t\tname of the sample. This must perfect match the prefix of each of the read count files. i.e. The sample name for file NA12878.200.B_1.uniq.cnt would be NA12878.200 (required)\n" .
			"-reference|r\t\t\tHLA reference fasta. There must also be an allele map file in the sample directory as the reference fa. (required)\n" .
			"\nGeneral Options:\n" .
			"-allele_refinement|ar\tSpecifies the level to which the predicted alleles are to be refined based on the observed reads (default:all)\n" .
			"   Possible values:\n" .
			"\trefine_only\t\tRefines the allelle call by looking predicting the true allele sequence using observed reads and looking for a better match in the reference\n" .
			"\tpredict_only\t\tReports if the observe reads support a novel allele sequence not found in the reference\n" .
			"\trefineAndPredict\tRefines the allele call (-refine_only) and report novel alleles (-novel_only)\n" . 
			"\tall\t\t\tRefines the allele call (-refine_only) and report novel alleles (-novel_only), creates a profile for the refined/novel allele sequence and calculates prediction metrics.\n" . 
			"\tnone\t\t\tTurns off refinement and novel allele prediction.\n" . 
			"-kraken_db|db\t\tbase directory of kraken database.\n" .
			"-kraken_path|kp\t\tbase directory of kraken installation. (default:base directory of path returned by `which kraken`)\n". 
			"-minimum_reads|min\tminimum number of reads from a gene before attempting to call HLA types.(default:100)\n" .
			"-output_dir|od\t\toutput directory(default:'.')\n" .
			"-threads|c\t\tnumber of threads (default:1)\n" .
			"-help|h\t\t\tprints this help prompt\n" .
				"\nAUTHORS:\n" . 
			"Martin Buchkovich:martin.buchkovich\@q2labsolutions.com\n" .
			"Chad Brown:chad.brown\@q2labsolutions.com\n" .
			"\nCREATED:\n$creation_date\n" .
			"\nLAST UPDATED:\n$last_updated\n" .
			"\nCopyright. Q2 Solutions|EA Genomics. 2016\n" .
			"\n";  
	my %predict_opts=(allele_refinement=>"none",profile_directory=>"",kraken_db=>".",kraken_path=>".",counts_directory=>"",reads_directory=>"", threads=>1, reference=>"", sample_name=>"", output_dir=>".", threads=>1, num_reads=>500000,read_length=>50, max_insert=>1000, seed=>1234, scale=>80, shape=>0.7, minimum_reads=>100);
	GetOptions(\%predict_opts, qw(allele_refinement|ar=s kraken_db|db kraken_path|kp profile_directory|sdir=s counts_directory|cd=s reads_directory|rd=s sample_name|sn=s reference|r=s output_dir|od=s threads|c=s help|h minimum_reads|min=s));
	if($predict_opts{help}){
		print "$predict_usage\n";
		exit;
	}elsif(! (-e $predict_opts{profile_directory} && -d $predict_opts{profile_directory})){
		print "Please specific existing profile directory.\n$predict_usage\n";
		exit; 
	}elsif(! (-e $predict_opts{counts_directory} && -d $predict_opts{counts_directory})){
		print "Please specify existing counts directory.\n$predict_usage\n";
		exit;
	}elsif(! (-e $predict_opts{reads_directory} && -d $predict_opts{reads_directory})){
		print "Please specify existing reads directory.\n$predict_usage\n";
		exit;
	}elsif(! (-e $predict_opts{reference})){
		print "Please specify a valid reference file using -reference.\n$predict_usage\n";
	}elsif($predict_opts{sample_name} eq ""){
		print "Please specify sample name.\n$predict_usage\n";
		exit;
	}else{
		my $allele_map = $predict_opts{reference};
		$allele_map =~ s/\.fa/\.allele_map.txt/;
		if(! -e $allele_map){
			print "Could not find the allele map $allele_map in the reference directory.\n";
				exit;
			}
			my $out;
			open ($out, ">$predict_opts{output_dir}/$predict_opts{sample_name}.HLATypes.txt") or $out=*STDOUT;
			my %simOpts=(num_reads=>$predict_opts{num_reads},read_length=>$predict_opts{read_length},max_insert=>$predict_opts{max_insert},scale=>$predict_opts{scale},shape=>$predict_opts{shape},seed=>$predict_opts{seed});
			my $answers_ref = predictHLAType($predict_opts{profile_directory},$predict_opts{reads_directory},$predict_opts{counts_directory},$predict_opts{sample_name},$predict_opts{reference},$allele_map,$predict_opts{threads},$predict_opts{allele_refinement},\%simOpts,$predict_opts{output_dir},$predict_opts{kraken_db},$predict_opts{kraken_path},"$predict_opts{output_dir}/$predict_opts{sample_name}",$predict_opts{minimum_reads});
			my %answers=%{$answers_ref};
			for my $answer (sort keys %answers){
				print $out "$answers{$answer}";
			}
			close($out);
		}	
	}elsif($module =~m/^-*h(elp)*$/){
		print "$usage\n";
	}else{
		print "$module is not a valid HLAProfiler module.\n\n$usage\n";
		exit;
}
sub run_build{
	my $ref = shift;
	my %opts = %{$ref};
	my $mkdir_cmds = "";
	#Set-up database directory structure
	if(! (-e "$opts{output_dir}/$opts{database_name}" && -d "$opts{output_dir}/$opts{database_name}")){
		if(! mkdir "$opts{output_dir}/$opts{database_name}"){
			print STDERR "Fatal Error: Error creating directory $opts{output_dir}/$opts{database_name}\n";
			exit; 	
		}else{
			$mkdir_cmds .= "mkdir $opts{output_dir}/$opts{database_name}\n"
		}
	}
	if(! (-e "$opts{output_dir}/$opts{database_name}/data" && -d "$opts{output_dir}/$opts{database_name}/data")){
		if(! mkdir "$opts{output_dir}/$opts{database_name}/data"){
			print STDERR "Fatal Error: Error creating directory $opts{output_dir}/$opts{database_name}/data\n";
			exit; 	
		}else{
			$mkdir_cmds .= "mkdir $opts{output_dir}/$opts{database_name}/data\n"
		}
	}
	if(! (-e "$opts{output_dir}/$opts{database_name}/data/reference" && -d "$opts{output_dir}/$opts{database_name}/data/reference")){
		if(! mkdir "$opts{output_dir}/$opts{database_name}/data/reference"){
			print STDERR "Fatal Error: Error creating directory $opts{output_dir}/$opts{database_name}/data/reference\n";
			exit; 	
		}else{
			$mkdir_cmds .= "mkdir $opts{output_dir}/$opts{database_name}/data/reference\n"
		}
	}
	if(! (-e "$opts{output_dir}/$opts{database_name}/data/logs" && -d "$opts{output_dir}/$opts{database_name}/data/logs")){
		if(! mkdir "$opts{output_dir}/$opts{database_name}/data/logs"){
			print STDERR "Fatal Error: Error creating directory $opts{output_dir}/$opts{database_name}/data/logs\n";
			exit; 	
		}else{
			$mkdir_cmds .= "mkdir $opts{output_dir}/$opts{database_name}/data/logs\n"
		}
	}

	my $logfile = "$opts{output_dir}/$opts{database_name}/data/logs/HLAProfiler.build.log";
	open($log, ">$logfile");
	open(my $commands, ">$opts{output_dir}/$opts{database_name}/data/logs/HLAProfiler.build.commands.txt");
	print $commands "$mkdir_cmds\n";

	my $merge_dir = "$opts{output_dir}/$opts{database_name}/data/reference";
	##Merge duplicate sequences
	my $reference = mergeDuplicateAlleles($opts{reference}, $merge_dir,$opts{cwd},$log,$commands,$logfile);
	
	##Create the HLA Taxonomy
	createHLATaxonomy($reference,"$opts{output_dir}/$opts{database_name}",$log,$commands, $logfile);
	
	##Build the HLA database
	buildHLADatabase("$opts{output_dir}/$opts{database_name}/data/reference/hla.ref.forKraken.fa",$opts{transcripts},$opts{transcript_gtf},$opts{exclusion_bed},$opts{output_dir},$opts{database_name},$opts{k_mer},$opts{minimizer},$opts{threads},$opts{kraken_path},$log,$commands);
	
	##Simulate reads and creates the k-mer profiles 
	createKmerProfiles($reference,$opts{output_dir},$opts{database_name},"$opts{output_dir}/$opts{database_name}/data/",$opts{num_reads},$opts{read_length},$opts{filter_reads}, $opts{threads}, $opts{kraken_path}, $opts{intermediate_files}, $opts{max_insert}, $opts{scale},$opts{shape}, $opts{seed},$log, $commands);
}

sub run_predict{	
	my $ref = shift;
	my $allele_map = shift;
	my %opts = %{$ref};
	my %simOpts=(num_reads=>$opts{num_reads},read_length=>$opts{read_length},max_insert=>$opts{max_insert},scale=>$opts{scale},shape=>$opts{shape},seed=>$opts{seed});
	
	load "$SCRIPTS_DIR/modules/RunKraken.pm";
	my $fastq_prefix = $opts{fastq1};
	#Strip suffix and parent directories from the file name 
	$fastq_prefix =~ s/.*\///;
	$fastq_prefix =~ s/_1\.f(ast)*q(\.gz)*//;
	
	my $log = *STDERR;
	if($opts{log}){
		open($log, ">$opts{log}");
	}else{
		$opts{log}=""; 
	}
	
	#Set-up output directory structure with the sample name being the root folder of the structure
	if(!(-e $opts{output_dir} && -d $opts{output_dir})){
		print $log "Making directory  $opts{output_dir}...";
		if(! mkdir "$opts{output_dir}"){
			print $log "Fatal Error: Error creating directory $opts{output_dir}/\n";
			exit; 	
		}else{
			print $log "mkdir $opts{output_dir}\n";
			print $log "DONE\n";
		}
	}
	if(! (-e "$opts{output_dir}/$fastq_prefix" && -d "$opts{output_dir}/$fastq_prefix")){
		print $log "Making directory  $opts{output_dir}/$fastq_prefix...";
		if(! mkdir "$opts{output_dir}/$fastq_prefix"){
			print $log "Fatal Error: Error creating directory $opts{output_dir}/$fastq_prefix\n";
			exit; 	
		}else{
			print $log "mkdir $opts{output_dir}/$fastq_prefix\n";
			print $log "DONE\n";
		}
	}
	if(!(-e "$opts{output_dir}/$fastq_prefix/filtered" && -d "$opts{output_dir}/$fastq_prefix/filtered")){
		print $log "Making directory  $opts{output_dir}/$fastq_prefix/filtered...";
		if(! mkdir "$opts{output_dir}/$fastq_prefix/filtered"){
			print $log "Fatal Error: Error creating directory $opts{output_dir}/$fastq_prefix/filtered\n";
			exit; 	
		}else{
			print $log "mkdir $opts{output_dir}/$fastq_prefix/filtered\n";
			print $log "DONE\n";
		}
	}
	if(!(-e "$opts{output_dir}/$fastq_prefix/counts" && -d "$opts{output_dir}/$fastq_prefix/filtered")){
		print $log "Making directory  $opts{output_dir}/$fastq_prefix/counts...";
		if(! mkdir "$opts{output_dir}/$fastq_prefix/counts"){
			print $log "Fatal Error: Error creating directory $opts{output_dir}/$fastq_prefix/counts\n";
			exit; 	
		}else{
			print $log "mkdir $opts{output_dir}/$fastq_prefix/counts\n";
			print $log "DONE\n";
		}
	}

	##Filter reads using kraken 	
	my $command = filterReads($opts{database_dir}, $opts{database_name}, "$opts{output_dir}/$fastq_prefix/filtered",$fastq_prefix, $opts{threads},$opts{fastq1}, $opts{fastq2},$opts{kraken_path}, $opts{log});	
	#print $log "$command\n";

	##Count the filtered reads 
	countUniqueReads("$opts{output_dir}/$fastq_prefix/filtered","$opts{output_dir}/$fastq_prefix/counts",$fastq_prefix,$opts{threads},$opts{log});
	
	##Predict the HLA type
	my ($answer_ref,$command_ref) = predictHLAType("$opts{database_dir}/$opts{database_name}/data/kmer_profiles","$opts{output_dir}/$fastq_prefix/filtered","$opts{output_dir}/$fastq_prefix/counts",$fastq_prefix,$opts{reference},$allele_map,$opts{threads}, $opts{allele_refinement},\%simOpts,"$opts{output_dir}/$fastq_prefix/","$opts{database_dir}/$opts{database_name}",$opts{kraken_path},"$opts{output_dir}/$fastq_prefix/$fastq_prefix",$opts{minimum_reads});
	my %answers = %{$answer_ref};
	my $out;
	
	##Prints output
	open ($out, ">$opts{output_dir}/$fastq_prefix/$fastq_prefix.HLATypes.txt") or $out = *STDOUT;
	print $out "Allele1_Accession\tAllele2_Accession\tAllele1\tAllele2\tProportion_reads\tProportion_signal\tCorrelation\tError\tPair_score\tFinal_score\tAllele1 Comments\tAllele2 Comments\n";
	for my $answer (sort keys %answers){
		print $out "$answers{$answer}";
	}
	close($out);	

	##Cleans intermediate file
	if(! $opts{intermediate_files}){
		print $log "Cleaning out intermediate files...";
		system("rm -rf $opts{output_dir}/$fastq_prefix/filtered");
		print $log "DONE";
	}
}

sub mergeDuplicateAlleles{
	load "$SCRIPTS_DIR/modules/MergeDuplicates.pm";
	my $reference_fasta = shift;
	my $merge_dir = shift;
	my $cwd_file = shift;
	my $log = shift;
	my $commands = shift;
	my $logfile = shift;
	
	my $reference_base = $reference_fasta;
	$reference_base=~s/.*\///;
	$reference_base=~s/\.fa(sta)*$//;
	my $output_prefix = "$merge_dir/$reference_base";

	print $commands "perl $SCRIPTS_DIR/modules/MergeDuplicates.pm $reference_fasta $output_prefix $cwd_file $logfile\n";
	print $log "Calling merge duplicates module...";
	
	##Calls the merge duplicate modules
	MergeDuplicates->mergeDuplicates($reference_fasta, $output_prefix,$cwd_file, $logfile);
	print $log "DONE.\n";
	return "$output_prefix.merged.fa";
} 

sub createHLATaxonomy{
	load "$SCRIPTS_DIR/modules/HLATaxonomy.pm";
	load "$SCRIPTS_DIR/modules/TaxonomyDivisions.pm";
	my $fa = shift;
	my $output_dir = shift;
	my $log = shift;
	my $commands = shift;
	my $logfile = shift;

	print $commands "perl $SCRIPTS_DIR/modules/HLATaxonomy.pm $fa $output_dir $logfile\n";
	print $log "Creating the HLA Taxonomy...";
	
	##Creates the taxonomy
	HLATaxonomy->createHLATaxonomy($fa, $output_dir, $logfile);
	print $log "DONE\n";

	print $commands "perl $SCRIPTS_DIR/modules/TaxonomyDivisions.pm $output_dir/taxonomy/names.dmp $output_dir/taxonomy/nodes.dmp $output_dir/taxonomy/divisions.txt $logfile\n";
	print $log "Creating the Taxonomy Divisions file...";
	
	##Creates the divisions file needed for filtering fastqs
	TaxonomyDivisions->createTaxonomyDivisions("$output_dir/taxonomy/names.dmp","$output_dir/taxonomy/nodes.dmp","$output_dir/taxonomy/divisions.txt", $logfile);
	print $log "DONE\n";
}

sub buildHLADatabase{
	my $reference = shift;
	my $transcript = shift;
	my $transcript_gtf = shift;
	my $exclusion_bed = shift;
	my $output_dir = shift;
	my $database_name = shift;
	my $kmer = shift;
	my $minimizer = shift;
	my $threads = shift;
	my $kraken_path = shift;
	my $distractome = "$output_dir/$database_name/data/reference/distractome.fa";
	my $log = shift;
	my $commands = shift;
	
	my $build_log = "$output_dir/$database_name/data/logs/HLAProfiler.build.database.log";
	print $log "Building database. Database specific log file can be found in $build_log\n";
	
	#Creates the distractome
	createDistractome($transcript,$transcript_gtf,$exclusion_bed, $distractome,$reference, $build_log, $log, $commands);

	#This will clean out the library if the database already exists so that duplicate sequences are not introduced
	if(-e "$output_dir/$database_name/library/added" && -d "$output_dir/$database_name/library/added"){
		print $commands "rm $output_dir/$database_name/library/added/*.fna\n";
		print $log "Removing existing fasta files in library...";
		unlink glob("$output_dir/$database_name/library/added/*.fna");
		print $log "DONE\n";
	}

	load "$SCRIPTS_DIR/modules/RunKraken.pm";
	
	##Add merged hla reference to library
	print $commands "perl $SCRIPTS_DIR/modules/RunKraken.pm -module library -fasta $reference -database_directory $output_dir -database_name $database_name -log $build_log -kraken_path $kraken_path\n";
	print $log "Adding reference fasta to library...";
	RunKraken->addLibraryToKraken($output_dir,$database_name,$reference,$kraken_path,$build_log);
	print $log "DONE\n";
	
	##Add distractome to library
	print $commands "perl $SCRIPTS_DIR/modules/RunKraken.pm -module library -fasta $distractome -database_directory $output_dir -database_name $database_name -log $build_log -kraken_path $kraken_path\n";
	print $log "Adding distractome to library...";
	RunKraken->addLibraryToKraken($output_dir,$database_name,$distractome,$kraken_path,$build_log);
	print $log "DONE\n";

	if(! -X "$output_dir/$database_name/taxonomy/gi_taxid_nucl.dmp"){
		print $commands "echo \" \" >$output_dir/$database_name/taxonomy/gi_taxid_nucl.dmp\n";
		print $log "Touching $output_dir/$database_name/taxonomy/gi_taxid_nucl.dmp...";
		open(GI, ">$output_dir/$database_name/taxonomy/gi_taxid_nucl.dmp");
		print GI " ";
		close(GI);
		print $log "DONE\n";
	}
	my $build_options = "--minimizer-len $minimizer --kmer-len $kmer --jellyfish-hash-size 500000";
	
	##This will clean out existing database files to prevent errors during building and to trigger a completely new build
	cleanDatabaseFiles($output_dir,$database_name,$log,$commands);
	
	##Build the database with Kraken
	print $commands "perl $SCRIPTS_DIR/modules/RunKraken.pm -module build -database_directory $output_dir -database_name $database_name -threads $threads -build_options \"$build_options\" -log $build_log -kraken_path $kraken_path\n";
	print $log "Building the HLA Database...";
	RunKraken->buildKraken($output_dir,$database_name,$build_options,$threads,$kraken_path,$build_log); 
	print $log "DONE\n";
}

sub cleanDatabaseFiles{
	my $output_dir = shift;
	my $database_name = shift;
	my $log = shift;
	my $commands = shift;
	
	if(-e "$output_dir/$database_name/databse.kdb"){
		print $commands "rm $output_dir/$database_name/databse.kdb\n";
		print $log "Removing existing database file $output_dir/$database_name/databse.kdb...";	
		unlink "$output_dir/$database_name/databse.kdb";
		print $log "DONE\n";
	}	
	if(-e "$output_dir/$database_name/databse.jdb"){
		print $commands "rm $output_dir/$database_name/databse.jdb\n";
		print $log "Removing existing database file $output_dir/$database_name/databse.jdb...";	
		unlink "$output_dir/$database_name/databse.jdb";
		print $log "DONE\n";
	}	
	if(-e "$output_dir/$database_name/databse.idx"){
		print $commands "rm $output_dir/$database_name/databse.idx\n";
		print $log "Removing existing database file $output_dir/$database_name/databse.idx...";	
		unlink "$output_dir/$database_name/databse.idx";
		print $log "DONE\n";
	}	
	if(-e "$output_dir/$database_name/lca.complete"){
		print $commands "rm $output_dir/$database_name/lca.complete\n";
		print $log "Removing existing database file $output_dir/$database_name/lca.complete...";	
		unlink "$output_dir/$database_name/lca.complete";
		print $log "DONE\n";
	}	
	if(-e "$output_dir/$database_name/gi2seqid.map"){
		print $commands "rm $output_dir/$database_name/gi2seqid.map\n";
		print $log "Removing existing database file $output_dir/$database_name/giseqid.map...";	
		unlink "$output_dir/$database_name/gi2seqid.map";
		print $log "DONE\n";
	}	
	if(-e "$output_dir/$database_name/seqid2taxid.map"){
		print $commands "rm $output_dir/$database_name/seqid2taxid.map\n";
		print $log "Removing existing database file $output_dir/$database_name/seqid2taxid.map...";	
		unlink "$output_dir/$database_name/seqid2taxid.map";
		print $log "DONE\n";
	}	
}

sub createDistractome{
	my $transcript=shift;
	my $transcript_gtf=shift;
	my $exclusion_bed = shift;
	my $distractome=shift;
	my $build_log = shift;
	my $log = shift;
	my $commands = shift;
	load "$SCRIPTS_DIR/modules/HLADistractome.pm";
	
	print $commands "perl $SCRIPTS_DIR/modules/HLADistractome.pm -transcript_fa $transcript -transcript_gtf $transcript_gtf -exclusion_bed $exclusion_bed -output_fa $distractome -log $build_log\n";
	print $log "Creating distractome...";

	##Creates the list of genes to exclude
	HLADistractome->findExcludeGenes($transcript_gtf,$exclusion_bed, $build_log);	
	
	##Creates the distractome
	HLADistractome->createDistractome($transcript, $distractome, $build_log);	
	print $log "DONE\n";
	
}

sub createKmerProfiles{
	my $reference = shift;
	my $db_dir = shift;
	my $database_name = shift;
	my $output_dir = shift;
	my $num_reads = shift;
	my $read_length = shift;
	my $filter_reads = shift;
	my $threads = shift;
	my $kraken_path = shift;
	my $keep_files = shift;
	my $max_insert = shift;
	my $scale = shift;
	my $shape = shift;
	my $seed = shift;
	my $log = shift;
	
	my $filtered_output = "$output_dir/simulated";
	if($filter_reads == 1){
		$filtered_output = "$output_dir/filtered";
	}	
	
	my $profile_log = "$output_dir/logs/HLAProfiler.build.profile.log";
	my $commands_file = "$output_dir/logs/HLAProfiler.build.profile.commands.txt";	
	my $slog;
	my $scommands;	
	open($slog, ">$profile_log") || (print  $log "Cannot open log file $profile_log. Exiting.\n" and exit);
	open($scommands, ">$commands_file") || (print $log "Cannot open commands file $commands_file. Exiting.\n" and exit);
	

	print $log "Creating K-mer profiles. Profile specific commands are captured in $output_dir/logs/HLAProfiler.build.profile.commands.txt and logged in $output_dir/logs/HLAProfiler.build.profile.log\n";

		
	my $output_prefix = "$output_dir/simulated/simulatedReads";
	load "$SCRIPTS_DIR/modules/SimulateReads.pm";
	
		
	print $log "Setting simulation parameters...";
	print $slog "Setting simulation parameters...";
	SimulateReads->setSimulationOptions($reference, $num_reads, $read_length,$max_insert, $scale, $shape, $seed, $threads, $SCRIPTS_DIR);
	print $log "DONE\n";	
	print $slog "DONE\n";	

	if(! (-e "$output_dir/simulated/" && -d "$output_dir/simulated/")){
		print $scommands "mkdir $output_dir/simulated/\n";
		print $log "Making directory $output_dir/simulated/...";
		print $slog "Making directory $output_dir/simulated/...";
		mkdir "$output_dir/simulated/";
		print $log "DONE\n";
		print $slog "DONE\n";
	} 
	
	##Load the HLA database into memory in order to speed up downstream filtering
	print $scommands "cat $db_dir/$database_name/database.* /dev/null\n";
	print $log "Loading database into memory...";
	print $slog "Loading database into memory...";
	`cat $db_dir/$database_name/database.* /dev/null`;
	print $log "DONE\n";
	print $slog "DONE\n";
	print $log "Simulating Reads...";


	
	if(! (-e "$output_dir/logs/profile" && -d "$output_dir/logs/profile")){
		print $scommands "mkdir $output_dir/logs/profile\n";
		print $log "Making directory $output_dir/logs/profile/...";
		print $slog "Making directory $output_dir/logs/profile/...";
		mkdir "$output_dir/logs/profile/";
		print $log "DONE\n";
		print $slog "DONE\n";
	} 
	
	if(! (-e "$output_dir/filtered" && -d "$output_dir/filtered")){
		print $scommands "mkdir $output_dir/filtered\n";
		print $log "Making directory $output_dir/filtered/...";
		print $slog "Making directory $output_dir/filtered/...";
		mkdir "$output_dir/filtered/";
		print $log "DONE\n";
		print $slog "DONE\n";
	} 

	
	my $dh;
	open($dh, $reference) || die "Cannot open reference file: $reference\n";
	my $line  = <$dh>;
	chomp $line;
	my $fn = $line;
	my $fastaString = "";
	my $fm = Parallel::ForkManager->new($threads);
	my $allele_commands = "";

	##Iterate through the lines of the reference fasta and simulate reads, filter reads, and count reads for each allele
	FASTA:while($line  = <$dh>){
		chomp $line;
		if(substr($line, 0, 1) eq ">") {
			my $prev_fn = $fn;
			my $prev_fastaString = $fastaString;
			$fn = $line;
			$fastaString = "";
				
			my $pid = $fm->start and next FASTA;
				$allele_commands .= simulateAndProcessAllele($prev_fn, $prev_fastaString, $filtered_output, $db_dir, $database_name, $output_dir, $kraken_path, $filter_reads, $keep_files,"-numReads $num_reads -readLength $read_length -max_insert $max_insert -scale $scale -shape $shape -seed $seed -threads $threads");
				print $slog "Simulating and Processing allele $prev_fn...DONE\n";
				print $scommands "$allele_commands";
			$fm->finish;
		} else {
			$fastaString=$fastaString.$line;
		}
	}
	$fm->wait_all_children;
	$allele_commands .= simulateAndProcessAllele($fn, $fastaString, $filtered_output, $db_dir, $database_name, $output_dir, $kraken_path, $filter_reads, $keep_files);		
	print $scommands $allele_commands;
	print $slog "Simulating and Processing allele $fn...DONE\n";
	print $log "DONE\n";

	##Determine the k-mer profile using simulated and filtered read counts
	print $log "Determining Profiles...";
	print $slog "Determining Profiles...";
	determineProfile("$output_dir/filtered",$output_dir,"simulatedReads",$threads,$slog,$scommands,1,0,$SCRIPTS_DIR);
	print $log "DONE\n";
	print $slog "DONE\n";	

	close($slog);
	close($scommands);
}

sub simulateAndProcessAllele{
	my $fn = shift;
	my $fastaString = shift;
	my $filtered_output = shift;
	my $db_dir = shift;
	my $database_name = shift;
	my $output_dir = shift;
	my $kraken_path = shift;
	my $filter_reads = shift;
	my $keep_files = shift;
	my $sim_options = shift;	

	my @parts = split /[\t ]/, $fn;
	my $acc = $parts[0];
	$acc =~ s/>HLA://;
	$acc =~ s/>IPD://;

	my $gene = $parts[1];
	$parts[1]=~m/([A-Z0-9]*)\*/;
	$gene = $1;
	
	my $allele_commands = "";
	if(!(-e "$output_dir/logs/profile/$gene" && -d "$output_dir/logs/profile/$gene")){
		mkdir "$output_dir/logs/profile/$gene";
		$allele_commands .= "mkdir $output_dir/logs/profile/$gene\n";
	}


	my $logfile = "$output_dir/logs/profile/$gene/$acc.profile.log";
	open(my $log, ">$logfile");	

	#Directory structure set-up
	if(!(-e "$output_dir/simulated/$gene" && -d "$output_dir/simulated/$gene")){
		print $log "Creating directory $output_dir/simulated/$gene/...";
		mkdir "$output_dir/simulated/$gene";
		print $log "DONE\n";
		$allele_commands .= "mkdir $output_dir/simulated/$gene\n";
	}
	if(!(-e "$output_dir/simulated/$gene/$acc" && -d "$output_dir/simulated/$gene/$acc")){
		print $log "Creating directory $output_dir/simulated/$gene/$acc/...";
		mkdir "$output_dir/simulated/$gene/$acc";
		print $log "DONE\n";
		$allele_commands .= "mkdir $output_dir/simulated/$gene/$acc\n";
	}
	if(!(-e "$output_dir/filtered/$gene" && -d "$output_dir/filtered/$gene")){
		print $log "Creating directory $output_dir/filtered/$gene/...";
		mkdir "$output_dir/filtered/$gene";
		print $log "DONE\n";
		$allele_commands .= "mkdir $output_dir/filtered/$gene\n";
	}
	if(!(-e "$output_dir/filtered/$gene/counts" && -d "$output_dir/filtered/$gene/counts")){
		print $log "Creating directory $output_dir/filtered/$gene/counts/...";
		mkdir "$output_dir/filtered/$gene/counts";
		print $log "DONE\n";
		$allele_commands .= "mkdir $output_dir/filtered/$gene/counts\n";
	}
	

	##Simulate the read
	my $output_prefix = "$output_dir/simulated/$gene/$acc/simulatedReads";
	$allele_commands .= "perl $SCRIPTS_DIR/modules/SimulateReads.pm -sequence $fastaString -name $fn -outPrefix $output_prefix $sim_options\n";
	print $log "Simulating reads for $acc...";
	my $allele = $acc;
	$allele = SimulateReads->makeSim($fn, $fastaString,$output_prefix,$logfile);
	print $log "DONE\n";		

	##Filter the read if that is the option
	my $filtered_prefix = "simulatedReads.$allele";		
	if($filter_reads == 1){
		if(!(-e "$output_dir/filtered/$gene/$acc" && -d "$output_dir/filtered/$gene/$acc")){
			print $log "Creating directory $output_dir/filtered/$gene/$acc/...";
			mkdir "$output_dir/filtered/$gene/$acc";
			print $log "DONE\n";
			$allele_commands .= "mkdir $output_dir/filtered/$gene/$acc\n";
		}
		print $log "Filtering $acc...\n";
		$allele_commands .= filterReads($db_dir, $database_name, "$output_dir/filtered/$gene/$acc/", $filtered_prefix, 1,"$output_prefix.${allele}_1.fastq","$output_prefix.${allele}_2.fastq",$kraken_path, $logfile);
		print $log "DONE\n";
	}
	
	##Count the reads
	print $log "Counting reads...";
	$allele_commands .= countUniqueReads("$filtered_output/$gene/$acc","$filtered_output/$gene/$acc",$filtered_prefix,1,$logfile);
	print $log "DONE\n";

	##Moves the reads to a new directory for read ability and better file IO performance
	print $log "Moving count files to $output_dir/filtered/$gene/counts...";
	move("$output_dir/filtered/$gene/$acc/${filtered_prefix}.${gene}_1.uniq.cnts","$output_dir/filtered/$gene/counts/${filtered_prefix}.${gene}_1.uniq.cnts");
	move("$output_dir/filtered/$gene/$acc/${filtered_prefix}.${gene}_2.uniq.cnts","$output_dir/filtered/$gene/counts/${filtered_prefix}.${gene}_2.uniq.cnts");
	$allele_commands .= "mv $output_dir/filtered/$gene/$acc/${filtered_prefix}_[12].uniq.cnt $output_dir/filtered/$gene/counts/\n";
	print $log "DONE\n";
		
	##Will clean up intermediate files if not flagged by -if to be saved
	if(! $keep_files){	
		if(-e "$output_dir/filtered/$gene/$acc/" && -d "$output_dir/filtered/$gene/$acc/"){
			$allele_commands .= "rm -rf $output_dir/filtered//$gene/$acc/\n";
			system("rm -rf $output_dir/filtered/$gene/$acc/");
		}
		if(-e "$output_dir/simulated/$gene/$acc/" && -d "$output_dir/simulated/$gene/$acc/"){
			$allele_commands .= "rm -rf $output_dir/simulated/$gene/$acc/\n";
			system("rm -rf $output_dir/simulated/$gene/$acc/");
		}
	}else{
		if(-e "$output_dir/filtered/$gene/$acc/" && -d "$output_dir/filtered/$gene/$acc/"){
			$allele_commands .= "for i in `ls $output_dir/filtered/$gene/$acc/$filtered_prefix*.fq`; do gzip \$i; done\n";
			gzipFiles("$output_dir/filtered/$gene/$acc", "$filtered_prefix\..*\.fq");
		}
		$allele_commands .= "gzip $output_prefix.${allele}_1.fastq\n";
		gzip("$output_prefix.${allele}_1.fastq");
		$allele_commands .= "gzip $output_prefix.${allele}_2.fastq\n";
		gzip("$output_prefix.${allele}_2.fastq");
	}
	return $allele_commands;
	close($log);
}

sub gzip{
	my $file = shift;
	`gzip $file`;
}

sub gzipFiles{
	my $directory = shift;
	my $search_string = shift;
	opendir(my $dh, $directory);
	while (readdir $dh){
		if($_=~m/$search_string/ && $_!~m/.gz$/){
			gzip("$directory/$_");
		}
	}
	close($dh);
	
}

sub determineProfile{
	my $input_dir = shift;
	my $output_dir = shift;
	my $prefix = shift;
	my $threads = shift;
	my $log = shift;
	my $commands = shift;
	load "$SCRIPTS_DIR/modules/DetermineProfile.pm";

	if(! (-e "$output_dir/kmer_profiles/" && -d "$output_dir/kmer_profiles/")){
		print $commands "mkdir $output_dir/kmer_profiles/\n";
		print $log "Making directory $output_dir/kmer_profiles/..."; 
	mkdir "$output_dir/kmer_profiles/";
	print $log "DONE\n";
}

my $fm = Parallel::ForkManager->new($threads);
opendir(my $dh, "$input_dir");
##Iterate through files in the input directory
DIR:while (readdir $dh){
	if(-d "$input_dir/$_" && $_ ne "." && $_ ne ".."){ ##Only concerned about directorys because these are the genes
		my $pid = $fm->start and next DIR;
		my $protein = $_;
		print $commands "perl $SCRIPTS_DIR/modules/DetermineProfile.pm fastq $protein $input_dir/$protein/counts $prefix $output_dir/kmer_profiles\n";
		DetermineProfile->createProfile($protein, "$input_dir/$protein/counts", $prefix, "$output_dir/kmer_profiles");
		print $log "Creating profile for $protein...DONE\n";
		$fm->finish;
	}
}
$fm->wait_all_children;
closedir $dh; 
}

sub filterReads{
	my $db_dir = shift;
	my $database_name = shift;
	my $output_dir = shift;
	my $prefix = shift;
	my $threads = shift;
	my $fastq1 = shift;
	my $fastq2 = shift;
	my $kraken_path = shift;
	my $log = shift;
	if(! (-e "$output_dir" && -d "$output_dir")){
		mkdir "$output_dir";
	} 
	
	load "$SCRIPTS_DIR/modules/RunKraken.pm";
	RunKraken->filterReads($db_dir, $database_name, "$output_dir/$prefix", $threads,$fastq1,$fastq2,$kraken_path,"", $log);		
	return "perl $SCRIPTS_DIR/modules/RunKraken.pm -module filter -database_directory $db_dir -database_name $database_name -output_prefix $output_dir/filtered/$prefix -threads $threads -fastq1 $fastq1 -fastq2 $fastq2 -kraken_path $kraken_path -log $log\n";
}


sub countUniqueReads{
	my $directory = shift;
	my $out_directory =shift;
	my $prefix = shift;
	my $threads = shift;
	my $log = shift;
	load "$SCRIPTS_DIR/modules/ReadCounter.pm";
	opendir(my $dh, $directory);
	my $fm = Parallel::ForkManager->new($threads);
	my $count_commands = "";
	##Iterate through files in the directory
	FILE:while (readdir $dh){
		if($_=~m/($prefix\..*_[12]).(fq|fastq|fq\.gz|fastq\.gz)$/){ ##Find fastq files
			my $file=$1;
			my $pid = $fm->start and next FILE;
			ReadCounter->count_reads("$directory/$_", "$out_directory/$file.uniq.cnts",$log);
			$count_commands .= "perl $SCRIPTS_DIR/modules/ReadCounter.pm $directory/$_ $out_directory/$file.uniq.cnts $log\n";
			$fm->finish;
		}	
	}
	$fm->wait_all_children;
	closedir $dh;	
	return ($count_commands);
}


sub predictHLAType{
	my $profile_directory = shift;
	my $reads_dir = shift;
	my $counts_dir = shift;
	my $prefix = shift;
	my $reference = shift;
	my $allele_map = shift;
	my $threads = shift;
	my $predict_mode = shift;
	my $opts_ref = shift;
	my $output_dir = shift;
	my $kraken_db = shift;
	my $kraken_path = shift;
	my $log_prefix = shift;
	my $minReads = shift;
	load "$SCRIPTS_DIR/modules/HLAPredict.pm";
	HLAPredict->setOptions(2,20,100,1,20,$threads,0,$output_dir, $SCRIPTS_DIR, $predict_mode, $opts_ref,$kraken_db,$kraken_path,$minReads);
	opendir(my $dh, $counts_dir);
	my %answers;
	my $fm = Parallel::ForkManager->new($threads);

	##This section of code is necessary for the child processes spawned by fork manager to return data directly back to the program. This will run everytime a child process finishes and store the gene predictions in the prediction hash
	$fm -> run_on_finish (
		sub{
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $result) = @_;
			my ($gene, $answer)=@$result;
			$answers{$gene}=$answer;

		}
	);
	FILE:while (readdir $dh){
		my $file = $_;
		if($file=~m/$prefix\.(.*)_1.uniq.cnts$/){
			my $gene = $1;
			if($gene!~m/Class/ && $gene!~m/root/ && -e "$profile_directory/$gene.profile.ph"){
				print STDERR "Starting gene $gene\n";
				my $pid = $fm->start and next FILE;
				my $answer = HLAPredict->predictHLAType("$counts_dir/$file","$reads_dir","$profile_directory/$gene.profile.ph",$reference,$allele_map, "$log_prefix.$gene.prediction.log");
			
				my @result = ();
				@result = ($gene, $answer);
				$fm->finish(0,\@result); ##This is the syntax necessary to return results to run_on_finish
			}
		}
	}
	closedir $dh;
	$fm->wait_all_children;
	return (\%answers);	
}


