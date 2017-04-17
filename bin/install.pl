#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';

my %opts = (install_dir=>".",mode=>"wget",bin_dir=>"");
GetOptions(\%opts, qw(install_dir|d=s mode|m=s bin_dir|b=s));

if($opts{mode}=~m/wget/i){
	installWget();
}else{
	print STDERR "Please specify an installation mode using -mode|m option. Valid options are git or wget.\n";
}

sub installWget{
	$opts{install_dir} = abs_path($opts{install_dir});
	print "Moving into installation directory (\"cd $opts{install_dir}\")...";
	chdir "$opts{install_dir}";
	print "DONE\n";
	print "Downloading Jellyfish v1.1.10 (\"wget https://github.com/gmarcais/Jellyfish/archive/v1.1.10.tar.gz jellyfish.tar.gz\")...";
	#my $jellyfish = `wget https://github.com/gmarcais/Jellyfish/archive/v1.1.10.tar.gz -O jellyfish.tar.gz`;
	my $jellyfish = `wget http://www.cbcb.umd.edu/software/jellyfish/jellyfish-1.1.11.tar.gz -O jellyfish.tar.gz`;
	#if($jellyfish ne ""){
	#	print "\n$jellyfish\n";
	#}
	print "DONE\n";
	print "Unpacking tarball jellyfish.tar.gz (\"tar -xzf jellyfish.tar.gz\")...";
	`tar -xzf jellyfish.tar.gz`;
	print "DONE\n";
	installJellyfish("jellyfish-1.1.11");
	
	print "Moving into installation directory (\"cd $opts{install_dir}\")...";
	chdir "$opts{install_dir}";
	print "DONE\n";
	print "Downloading kraken (\"wget https://github.com/ExpressionAnalysis/kraken/archive/v0.10.5-beta.1-ea.tar.gz\")..."; 	
	`wget https://github.com/ExpressionAnalysis/kraken/archive/v0.10.5-beta-ea.1.tar.gz -O kraken_ea.tar.gz`;
	print "DONE\n";
	print "Unpacking tarball kraken_ea.tar.gz (\"tar -xzf kraken_ea.tar.gz\")...";
	`tar -xzf kraken_ea.tar.gz`;
	print "DONE\n";
	installKraken("kraken-0.10.5-beta-ea.1");	
}

sub installJellyfish{
	my $directory = shift;
	$directory = abs_path($directory);
	print "Moving into $directory...";
	chdir $directory;
	print "DONE\n";
	if(! -e "configure"){
		print "./configure not found. Running autoreconf (autoreconf -i -I /usr/share/aclocal)...";
		`autoreconf -i -I /usr/share/aclocal`;
		print "DONE\n";
	}
	print "Installing jellyfish\n";
	print "\t./configure --prefix=$directory\n";
	`./configure --prefix=$directory\n`;
	print "\tmake\n";
	`make`;
	print "\tmake install\n";
	`make install`;
	if($opts{bin_dir} eq ""){
		print "A PATH directory has not been specified with option -bin_dir|b. Either add $directory/bin/ to your path, or copy the $directory/bin/jellyfish into a folder in your path\n";
	}else{
		print "Copying executable to $opts{bin_dir}(cp $directory/bin/jellyfish $opts{bin_dir}\n...";
		`cp $directory/bin/jellyfish $opts{bin_dir}`;
	}
	print "INSTALLATION COMPLETE\n";
}

sub installKraken{
	my $directory = shift;
	$directory = abs_path($directory);
	print "Moving into $directory...";
	chdir $directory;
	print "DONE\n";
	print "Installing kraken (\"./install_kraken.sh .\")...";
	`./install_kraken.sh .`;
	print "DONE\n";
	if($opts{bin_dir} eq ""){
		print "A PATH directory has not been specified with option -bin_dir|b. Either add $directory/bin/ to your PATH, or copy $directory/bin/kraken and $directory/bin/kraken-build into a folder in your path\n";
	}else{
		print "Copying kraken to $opts{bin_dir}(cp $directory/kraken $opts{bin_dir}\n...";
		`cp $directory/kraken $opts{bin_dir}`;
		print "Copying kraken-build to $opts{bin_dir}(cp $directory/bin/kraken-build $opts{bin_dir}\n...";
		`cp $directory/kraken-build $opts{bin_dir}`;
	}
	print "INSTALLATION_COMPLETE\n";
}
