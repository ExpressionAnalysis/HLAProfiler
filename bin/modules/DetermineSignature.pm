#!/usr/bin/env perl
package DetermineSignature;

use Storable;
#use Statistics::R;


my %kmers = ();
my $kmer_keys = ();
my $allele_counts = ();
my %names = ();

sub runCommandline{
	my $module = shift;
	if($module eq "DetermineSignature"){
		$module = shift;
	}
	if($module eq "fastq"){
		createSignature(@_);
	}elsif($module eq "rdata"){
		my $rdata = shift;
		my $protein = shift;
		my $output_directory = shift;
		#printSignatureToHashFile(rdataToSig($rdata),"$output_directory/$protein.signature.ph");
	}
}
sub rdataToSig{
	my $rdata = shift;	
	
	#my $R = Statistics::R->new();
	#$R->start();
	#$R->run("load(\"$rdata\")");
	#my $bias = $R->run(bias);
	#my @blines = split /\n/, $bias;
	#my $sig = $R->get(sig);
	#my @slines = split /\n/, $sig;
	#my $tmp =1;
}

sub createSignature{
	my $protein = shift;
	if($protein eq "DetermineSignature"){
		$protein = shift;
	}
	my $simulation_directory = shift;
	my $prefix = shift;
	my $output_directory = shift;
	my $sig = determineSignature($protein, $simulation_directory, $prefix);
	printSignatureToHashFile($sig, "$output_directory/$protein.signature.ph");
	#printSignatureToMatrixFile($sig, "$output_directory/$protein.signature.txt");
}

sub printSignatureToHashFile{
	my $sig_ref = shift;
	my $file = shift;
	store $sig_ref, $file;
}

sub printSignatureToMatrixFile{
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
}

sub determineSignature{
	my $protein = shift;
	if($protein eq "DetermineSignature"){
		$protein = shift;
	}
	my $simulation_directory = shift;
	my $prefix = shift;

	$simulation_directory =~ s/\/$//;

	my $kmer_cnt = 0;
	opendir(my $dh, "$simulation_directory/");
	while (readdir $dh){
		if($_=~m/$prefix\.(.*)\.${protein}_([12]).uniq.cnts/){
			my $allele = $1;
			my $read = $2;
			open(IN, "$simulation_directory/$_");
			while(<IN>){
				chomp;
				my @parts = split /\t/, $_;
				my $key = "";
				$parts[0] = revcomp($parts[0]) if ($read == 1);
				if(defined $kmer_keys{$parts[0]}){
					$key = $kmer_keys{$parts[0]};
				}else{
					$key = $kmer_cnt;
					$kmer_keys{$parts[0]}=$key;
					$names{$key}=$parts[0];
					$kmer_cnt++;
				}
				$kmers{$allele}{$key}+=$parts[1];
				#$kmers{$allele}{total}+=$parts[1];
				$allele_counts{$allele}+=$parts[1];
			}
			close(IN);
		}		
	}
	close($dh);
	my %signature = ();
	$signature{sig}=\%kmers;
	$signature{totals}=\%allele_counts;
	$signature{names}=\%names;
	$signature{kmers}=\%kmer_keys;
	return(\%signature);	
}

sub updateSignature{
	my $protein = shift;
	if($protein eq "DetermineSignature"){
		$protein = shift;
	}
	my $simulation_directory = shift;
	my $signature_file = shift;

}

sub revcomp{
	my $read = shift;
	my $rev ="";
	foreach my $base (split //, $read){
		$base=~tr/ATCG/TAGC/;
		$rev = $base . $rev;
	}
	return $rev;
}



__PACKAGE__->runCommandline(@ARGV) unless caller;
