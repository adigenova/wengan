
###############################################################################
# Author: Alex Di Genova
# Laboratory: ERABLE/INRIA
# Copyright (c)
# year: 2019
###############################################################################

=head1 NAME

B<Wengan> - An accurate and ultrafast genome assembler

=head1 SYNOPSIS

  # Assembling Oxford nanopore and illumina reads wiht WenganM
   wengan.pl -x ontraw -a M -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l ont.fastq.gz -p asm1 -t 20 -g 3000

  # Assembling PacBio reads and illumina reads with WenganA
   wengan.pl -x pacraw -a A -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l pac.fastq.gz -p asm2 -t 20 -g 3000

  # Assembling ultra-long nanopore reads and BGI reads with WenganM
   wengan.pl -x ontlon -a M -s lib2.fwd.fastq.gz,lib2.rev.fastq.gz -l ont.fastq.gz -p asm3 -t 20 -g 3000

  # Non-hybrid assembly of PacBio Circular Consensus Sequence data with WenganM
   wengan.pl -x pacccs -a M -l ccs.fastq.gz -p asm4 -t 20 -g 3000

  # Assembling ultra-long nanopore reads and Illumina reads with WenganD (requires a high memory machine 600Gb)
   wengan.pl -x ontlon -a D -s lib2.fwd.fastq.gz,lib2.rev.fastq.gz -l ont.fastq.gz -p asm5 -t 20 -g 3000

  # Assembling pacraw reads wiht pre-assembled short-read contigs from Minia3
   wengan.pl -x pacraw -a M -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l pac.fastq.gz -p asm6 -t 20 -g 3000 -c contigs.minia.fa

  # Assembling pacraw reads wiht pre-assembled short-read contigs from Abyss
   wengan.pl -x pacraw -a A -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l pac.fastq.gz -p asm7 -t 20 -g 3000 -c contigs.abyss.fa

  # Assembling pacraw reads wiht pre-assembled short-read contigs from DiscovarDenovo
   wengan.pl -x pacraw -a D -s lib1.fwd.fastq.gz,lib1.rev.fastq.gz -l pac.fastq.gz -p asm8 -t 20 -g 3000 -c contigs.disco.fa

=head1 DESCRIPTION


B<Wengan> is a new genome assembler that unlike most of the current long-reads assemblers avoid entirely the all-vs-all read comparison.
The key idea behind B<Wengan> is that long-read alignments can be B<inferred by building paths> on a sequence graph. To achieve this, B<Wengan> build a new sequence graph called the Synthetic Scaffolding Graph. The SSG is build from a spectrum of synthetic mate-pair libraries extracted from raw long-reads. Then, longer alignments are build by peforming a transitive reduction of the edges.
Another distinct feature of B<Wengan> is that is the only assembler that perform B<self-validation> by following the read information. B<Wengan> identify miss-assemblies at differents steps of the assembly process. For more information about the algorithmic ideas behind B<Wengan> please read the preprint available on bioRxiv.


=head2 About the name

B<Wengan> is a Mapudungun word. The Mapudungun is the language of the B<Mapuche> people, the largest indigenous inhabitants of south-central Chile. B<Wengan> means I<"Making the path">.
Thus, when you assemble a genome with wengan, you are literally I<making the genome path>.

=head1 AUTHOR - Alex Di Genova

Email digenova@gmail.com

=head1 SHORT-READ ASSEMBLY

B<Wengan> uses a de bruijn graph assembler to build the assembly backbone from short-read data.
Currently, B<Wengan> can use B<Minia3>, B<Abyss2> or B<DiscoVarDenovo>.  The recomended short-read coverage
is B<50-60X> of 150bp x 2 or 250bp x 2 short reads.

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
   print "$0 usage : -x preset -a M [M,A,D] -s short-reads[fwd1.fastq.gz,rev1.fastq.gz,fwd2.fastq.gz,rev2.fastq.gz..] -l long-reads.fq.gz -p prefix -t cpus -g 3 [genomesize in Mb] -c <assembled short-read contigs> -i <insert size lists> -n <show pipeline comands>\n";
   #system("perldoc $0");
   print "$0 -h for detailed usage information\n";
   exit 1;
}


my %opts = ();

getopts( "x:s:l:p:t:g:a:c:i:I:L:F:hn", \%opts );
#display help usage using the perldoc
if($opts{h}){
   system("perldoc $0");
   exit 0;
}

#mandatory varibles
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

# we check the pipeline called
my $pipeline=();

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
#print Dumper(%opts);
#we check if the user want to see the comand being executed
if($opts{n}){
  $pipeline->show_pipeline();
}else{
$pipeline->run();
}

#We need to check for fast-sg exec and kmc
#my $fhome=dirname($0);

=head1 LONGREADS PRESETS

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


# not sure if allow this

=head1 WENGAN ADVANCED OPTIONS

The following options allows to override the presets of B<Wengan> components.
Don't change this variables if you are not sure.


=head2 FastMin-SG

An alignment-free algorithm for ultrafast scaffolding graph construction from short or long reads.

=head2 IntervalMiss

IntervalMiss detect miss-assembled contigs and correct them when necessary.


=head2  Liger

Liger use the Synthetic Scaffoding Graph to compute overlap among long reads,
order and orient short contigs, validate scaffols sequences, fill the gaps and
polishing of the assembly.

=cut


#get dirname
sub dirname {
	my $prog = shift;
	return '.' unless ($prog =~ /\//);
	$prog =~ s/\/[^\s\/]+$//g;
	return $prog;
}
