#!/usr/bin/perl

###############################################################################
# Author: Alex Di Genova
# Laboratory: ERABLE/INRIA
# Copyright (c)
# year: 2019
###############################################################################

=head1 NAME

B<Wengan> - An accurate and ultrafast genome assembler

=head1 SYNOPSIS

  # Assembling Oxford nanopore and illumina reads with WenganM
   wengan.pl -x ontraw -a M -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l ont.fastq.gz -p asm1 -t 20 -g 3000

  # Assembling PacBio reads and illumina reads with WenganA
   wengan.pl -x pacraw -a A -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l pac.fastq.gz -p asm2 -t 20 -g 3000

  # Assembling ultra-long nanopore reads and BGI reads with WenganM
   wengan.pl -x ontlon -a M -s lib2.fwd.fastq.gz,lib2.rev.fastq.gz -l ont.fastq.gz -p asm3 -t 20 -g 3000

  # Non-hybrid assembly of PacBio Circular Consensus Sequence data with WenganM
   wengan.pl -x pacccs -a M -l ccs.fastq.gz -p asm4 -t 20 -g 3000

  # Assembling ultra-long nanopore reads and Illumina reads with WenganD (requires a high memory machine 600Gb)
   wengan.pl -x ontlon -a D -s lib2.fwd.fastq.gz,lib2.rev.fastq.gz -l ont.fastq.gz -p asm5 -t 20 -g 3000

  # Assembling pacraw reads with pre-assembled short-read contigs from Minia3
   wengan.pl -x pacraw -a M -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l pac.fastq.gz -p asm6 -t 20 -g 3000 -c contigs.minia.fa

  # Assembling pacraw reads with pre-assembled short-read contigs from Abyss
   wengan.pl -x pacraw -a A -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l pac.fastq.gz -p asm7 -t 20 -g 3000 -c contigs.abyss.fa

  # Assembling pacraw reads with pre-assembled short-read contigs from DiscovarDenovo
   wengan.pl -x pacraw -a D -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l pac.fastq.gz -p asm8 -t 20 -g 3000 -c contigs.disco.fa

=head1 DESCRIPTION


B<Wengan> is a new genome assembler that unlike most of the current long-reads assemblers avoid entirely the all-vs-all read comparison.
The key idea behind B<Wengan> is that long-read alignments can be B<inferred by building paths> on a sequence graph. To achieve this, B<Wengan> build a new sequence graph called the Synthetic Scaffolding Graph. The SSG is build from a spectrum of synthetic mate-pair libraries extracted from raw long-reads. Then, longer alignments are build by peforming a transitive reduction of the edges.
Another distinct feature of B<Wengan> is that is the only assembler that perform B<self-validation> by following the read information. B<Wengan> identify miss-assemblies at differents steps of the assembly process. For more information about the algorithmic ideas behind B<Wengan> please read the preprint available on bioRxiv.


=head2 ABOUT THE NAME

B<Wengan> is a Mapudungun word. The B<Mapudungun> is the language of the B<Mapuche> people, the largest indigenous inhabitants of south-central Chile. B<Wengan> means I<"Making the path">.
Thus, when you assemble a genome with B<Wengan>, you are I<making the genome path>.

=head1 AUTHOR - Alex Di Genova

Email digenova@gmail.com

=head1 SHORT-READ ASSEMBLY

B<Wengan> uses a de bruijn graph assembler to build the assembly backbone from short-read data.
Currently, B<Wengan> can use B<Minia3>, B<Abyss2> or B<DiscoVarDenovo>.  The recomended short-read coverage
is B<50-60X> of 2 x 150bp  or 2 x 250bp short reads.

=head2 WenganM [M]

This B<Wengan> mode use the B<Minia3> short-read assembler, this is the fastest mode of B<Wengan> and can assemble a complete human genome
in less than 210 CPU hours (~50Gb of RAM).

=head2 WenganA [A]

This B<Wengan> mode use the B<Abyss2> short-read assembler, this is the lowest memory mode of B<Wengan> and can assemble a complete human genome
in less than 40Gb of RAM (~900 CPU hours). This assembly mode takes  ~2 days when using 20 CPUs on a single machine.


=head2 WenganD [D]

This B<Wengan> mode use the B<DiscovarDenovo> short-read assembler, this is the greedier memory mode of B<Wengan> and for assembling a complete human genome need about 600Gb of RAM (~900 CPU hours).
This assembly mode takes ~2 days when using 20 CPUs on a single machine.

=cut


### We start coding
use strict;
use Data::Dumper;
use Getopt::Std;
use FindBin;
#load the wengan library
use lib "$FindBin::Bin/perl";
#local perl clases to control the wengan execution
use Wengan::Reads; # class to handle the read data.
use Wengan::Scheduler::Local; # the scheduler is make and control the execution end-to-end

sub usage {
   die(qq/
  Usage example :
    # Assembling Oxford nanopore and illumina reads with WenganM
    wengan.pl -x ontraw -a M -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l ont.fastq.gz -p asm1 -t 20 -g 3000

  Wengan options:

   Mandatory options***:
      -x preset [ontlon,ontraw,pacraw,pacccs]
      -a Mode [M,A,D]
      -s short-reads [fwd1.fastq.gz,rev1.fastq.gz..]
      -l long-reads.fq.gz
      -g 3000 [genome size in Mb]
      -p prefix

   General Options :
      -h [detail information]
      -t cores [1]
      -c <pre-assembled short-read contigs>
      -i <insert size lists>
      -n <show pipeline comands>;

   Advanced Options (Change the presets):
      FastMin-SG options:
        -k k-mer size [15-28]
        -w minimizer window [5-15]
        -q minimum mapping quality [20-60]
        -m moving window [150]
      IntervalMiss options:
        -d Minimum base coverage [def:7]
      Liger options:
        -M Minimum contig length in backbone [def:2000]
        -R Repeat copy number factor [def:1.5]
        -L Length of long mate-edges [def:100000]
        -N Number of long-read needed to keep a potencial erroneus mate-edge [def:5]
        -P Minimum length of reduced paths to convert them to physical fragments [def:20kb]

 $0 -h for detailed usage information.
   \n/);
   exit 1;
}

#x:s:l:p:t:g:a:c:i:k:w:q:m:d:M:L:N:P:hn

my %opts = ();

getopts( "x:s:l:p:t:g:a:c:i:k:w:q:m:d:M:L:N:P:R:hn", \%opts );
#display help usage using the perldoc
if($opts{h}){
   system("perldoc $0");
   exit 0;
}

#mandatory variables
if (!defined $opts{x} or !defined $opts{l} or !defined $opts{g}) {
   usage;
}

#object that handle the reads in fastq format and compressed with gzip
my $reads= new Wengan::Reads() or die "cannot create read object\n";
#print Dumper($reads);

#we check the preset for long-reads
if($opts{x} eq "ontlon"){
	wengan_ontlon($reads,%opts);
}elsif($opts{x} eq "ontraw"){
	wengan_ontraw($reads,%opts);
}elsif($opts{x} eq "pacraw"){
	wengan_pacraw($reads,%opts);
}elsif($opts{x} eq "pacccs"){
	wengan_pacccs($reads,%opts);
}else{
	print "Unkown preset $opts{x}\n";
	usage();
	exit 1;
}

if($opts{x} eq "pacccs" and $opts{a} ne "M"){
  print "Currently only the WenganM pipeline is available for Circular Consensus Sequence (CCS) data\n";
  exit 1;
}

#we define number of threads for the pipeline.
if(!defined $opts{t}){
  $opts{t}=1;#default is one
}


my $pipeline=();
# we check the pipeline called
if($opts{a} eq "M"){
      $pipeline=Wengan::Scheduler::Local->new($reads,"WenganM",%opts);
}elsif($opts{a} eq "A"){
      $pipeline=Wengan::Scheduler::Local->new($reads,"WenganA",%opts);
}elsif($opts{a} eq "D"){
      $pipeline=Wengan::Scheduler::Local->new($reads,"WenganD",%opts);
}else{
	print "Unkown pipeline $opts{a}\n";
  print "Available pipelines: WenganM, WenganA, WenganD\n";
	usage();
	exit 1;
}

#print Dumper($pipeline)
#we check if the user want to see the pipeline commands
if($opts{n}){
  $pipeline->show_pipeline();
}else{
$pipeline->run();
}


=head1 LONG-READ PRESETS

The presets define several variables of the wengan pipeline execution and depends on the long-read technology used to sequence the genome.
The recommended long-read coverage is 30X.

=cut

=head2	ontlon

preset for raw ultra-long-reads from Oxford Nanopore, tipically having an  N50 > 50kb.

=cut

sub wengan_ontlon{
	my ($reads,%opts)=@_;
	$reads->add_short_reads($opts{s});
	$reads->add_long_reads($opts{l});
}

=head2	ontraw

preset for raw long-reads Nanopore reads tipically having an  N50 ~[15kb-40kb].

=cut


sub wengan_ontraw{
  my ($reads,%opts)=@_;
	$reads->add_short_reads($opts{s});
	$reads->add_long_reads($opts{l});
  #print Dumper($reads);
}


=head2	pacraw

preset for raw long-reads from Pacific Bioscience (PacBio) tipically having an  N50 ~[8kb-60kb].

=cut


sub wengan_pacraw{
  my ($reads,%opts)=@_;
	$reads->add_short_reads($opts{s});
	$reads->add_long_reads($opts{l});
  #print Dumper($reads);

}


=head2	pacccs

preset for Circular Consensus Sequences from Pacific Bioscience (PacBio) tipically having an  N50 ~[15kb].

=cut

sub wengan_pacccs{
  my ($reads,%opts)=@_;
	$reads->add_long_reads($opts{l});
  #print Dumper($reads);
}


=head1 WENGAN ADVANCED OPTIONS

The following options allows to override the presets of B<Wengan> components.
Don't change this variables if you are not sure.


=head2 FastMin-SG

An alignment-free algorithm for ultrafast scaffolding graph construction from short or long reads.


=head3 FastMin-SG options (Override the presets)

   Indexing:
    -k INT       k-mer size (no larger than 28) [15]
    -w INT       minizer window size [10]

  Mapping synthetic read-pairs:

    -i list      Insert sizes for synthetic libraries [i.e 500,1000,2000,3000,4000,5000, .. ,20000]
    -q INT       Minimum quality score (no larger than 60) [40]
    -m INT       Moving windown [150]

=head2 IntervalMiss

IntervalMiss detect miss-assembled contigs and correct them when necessary.

=head3 IntervalMiss options

  -d INT Minimum base coverage [<10]


=head2  Liger

Liger use the Synthetic Scaffoding Graph to compute overlap among long reads,
order and orient short contigs, validate scaffols sequences, fill the gaps and
polishing of the assembly.

=head3 Liger options

Short-contigs options:

      -M INT     Minimum contig length in scaffolding [--mcs] (default=`2000', min=`1000')
      -R FLOAT   Repeat copy number factor [--rcn]  (default=`1.5')

Long-reads overlap options:

      -L         Length of long mate-edges [--lme] (default=`100000')

Validation of lines options:

      -N INT     Number of long-read needed to keep a potencial erroneus mate-edge [--nlm] (default=`5', min=`1')
      -P INT     Minimum length of reduced paths to convert them to physical fragments [--mlp] (default=`20000', min=`5000')

=cut

sub dirname {
	my $prog = shift;
	return '.' unless ($prog =~ /\//);
	$prog =~ s/\/[^\s\/]+$//g;
	return $prog;
}
