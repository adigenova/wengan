
###############################################################################
# Author: Alex Di Genova 
# Laboratory: ERABLE/Mathomics
# Copyright (c)
# year: 2017
###############################################################################
use Data::Dumper;
use Getopt::Std;
use strict;

sub usage {
	print "$0 usage : -a <list_of_fastq_files.txt> -b <memory mb> -c <cores> -p <prefix>\n";
	print "Error in use\n";
	exit 1;
}

my %opts = ();
getopts( "a:b:c:p:", \%opts );
if ( !defined $opts{a}  or !defined $opts{c} or !defined $opts{p} ) {
	usage;
}


#minia variables 
my $minia=$ENV{WROOTDIR}."/bin/minia";
my @kmers=(41,81,121);
my @abun=(2,2,2);
#user variables
my $cores=$opts{c};
my $prefix=$opts{p};
my $memory=5000;
# we set the memory if necesary
if(defined $opts{b} and $opts{b} > $memory){
	$memory=$opts{b};	
}
my $files=$opts{a};
#--kmer-sizes 41,81,121  --abundance-mins 2,2,2
my $reads=load_reads($opts{a});
#we run minia iteratively
for(my $i=0; $i< 3; $i++){
	my $newfile="round.$i.K$kmers[$i].txt";
	open(OUT,">".$newfile) or die "cannot open file $newfile\n";
	print OUT $_."\n" foreach(@$reads);
	if($i > 0){
      	 	my $contigs = $prefix.".".$kmers[$i-1].".contigs.fa";
		print OUT $contigs."\n"; 
		print OUT $contigs."\n"; 
		print OUT $contigs."\n"; 
	}
	close(OUT);
	run_minia($kmers[$i],$newfile,$prefix.".".$kmers[$i],$abun[$i],$memory,$cores,$minia);
}

sub load_reads{
	my ($file)=@_;
	open(FILE,$file) or die "cannot open file $file\n";
	my $files=();
	while(my $line=<FILE>){
		chomp $line;
		push(@{$files},$line);
	}
	return $files;
}


sub run_minia{
	my ($k,$files,$prefix,$abum, $mem, $cores, $minia)=@_;
	
	my @params=("-in",$files,
		"-kmer-size",$k,
		"-abundance-min",$abum, 
		"-out",$prefix,
		"-minimizer-size",10,
		"-max-memory",$mem,
		"-nb-cores",$cores,
		"2>",$prefix.".err",
		">",$prefix.".log");
	
	 my $cmd=join(" ",$minia,@params);
	 #system($cmd) or die "problem running $cmd\n";
	print $cmd."\n";
	system($cmd);
	if ($? == -1) {
    	print "failed to execute: $!\n";
   	exit(1);	
	}
	elsif ($? & 127) {
    	printf "child died with signal %d, %s coredump\n",
   	 ($? & 127),  ($? & 128) ? 'with' : 'without';
	exit(1);
	}
	else {
    		printf "child exited with value %d\n", $? >> 8;
	}

}


