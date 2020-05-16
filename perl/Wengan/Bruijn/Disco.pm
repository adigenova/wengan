package Wengan::Bruijn::Disco;

use strict;
use Wengan::Common::GlobalConfig qw(DiscoVarDenovo_BIN SEQTK_BIN);

sub new{
  my ($packagename,%opts) = @_;


  #minimum variables for Disco2
  my $self = {contigs=>undef,cores=>$opts{t},prefix=>$opts{p}, preset=>$opts{x}};
  #we ask if the contigs are passed
  if(defined $opts{c}){
    $self->{contigs}=$opts{c};
  }

  #we complain if the contigs are not given and the DiscoVarDenovo_BIN is not found
  if (!defined DiscoVarDenovo_BIN and !defined $opts{c}){
    die "DiscoVarDenovo binary not found\n";
  }
  if(!defined SEQTK_BIN ){
      die "Seqtk binary not found for contig post procesing\n";
  }

  bless ($self, $packagename);
  return ($self);
}



#generic function is called from the WenganM pipeline
sub create_jobs{
  my ($self,$reads)=@_;
  #minia contigs were create previously
  if(defined $self->{contigs}){
      my $job=();
      push(@{$job->{target}},$self->{prefix}.".contigs-disco.fa");
      my $cmd="ln -s $self->{contigs} $job->{target}[0]";
      push(@{$job->{cmds}},$cmd);
      push(@{$self->{jobs}},$job);
      $self->{contigs}=$job->{target}[0];
      #clean the contigs
      $self->_create_seqtk_jobs(@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
      push(@{$self->{main_target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
  }else{
    #default method for create jobs
    $self->_create_jobs_short($reads);
    #push(@{$self->{main_target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
  }

}
sub _def_parameters{
      my ($self,$rlen)=@_;
      #Disco configuration used on the wengan ms
      # we use the discovar varsion discovarexp-51885
      # available on  ftp://ftp.broadinstitute.org/pub/crd/nightly/discovardenovo/discovarexp/discovarexp-51885.tar.gz
      #export MALLOC_PER_THREAD=1
      #${DISCO} READS="SRR891258_{1,2}.fastq.gz,SRR891259_{1,2}.fastq.gz" NUM_THREADS=44 OUT_DIR=/tmp/60XNA12878
      #cp -ra /tmp/60XNA12878 .
      #

      my $param="NUM_THREADS=$self->{cores} OUT_DIR=/tmp/$self->{prefix}D";
      return $param;
}

sub _create_jobs_short{
    my ($self,$reads)=@_;
    #we get the parameters for Disco2
    #my $params = $self->_def_parameters($reads->{sreads}[0]->{len});
    #preparing the reads for Disco2
    my @rfiles=();
    foreach my $r (@{$reads->{sreads}}){
            push(@rfiles,"$r->{fwd},$r->{rev}");
    }
    my $readfiles=join(",",@rfiles);
    #now we create the asbyss cmd
    my @mopt=(
      "READS=\"$readfiles\"",
      "OUT_DIR=/tmp/$self->{prefix}D",
      "NUM_THREADS=$self->{cores}",
      "2>",$self->{prefix}.".Disco_denovo.err",
      ">",$self->{prefix}.".Disco_denovo.log");
      my $job=();
      #We run the disco program
      my $cmd="export MALLOC_PER_THREAD=1";
      push(@{$job->{cmds}},$cmd);
      push(@{$job->{cmds}},join(" ",DiscoVarDenovo_BIN,@mopt));
      $cmd="cp -a /tmp/$self->{prefix}D $self->{prefix}D";
      push(@{$job->{cmds}},$cmd);
      $cmd="ln -s $self->{prefix}D/a.final/a.lines.fasta $self->{prefix}.contigs-disco.fa";
      push(@{$job->{cmds}},$cmd);

      push(@{$job->{target}},$self->{prefix}.".contigs-disco.fa");
      push(@{$self->{jobs}},$job);
    # we add the after assembly steps
    $self->_create_seqtk_jobs(@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    push(@{$self->{main_target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    #return $self->{jobs};
}

sub _create_seqtk_jobs{
  my ($self,$contigs)=@_;
  my $job=();
  push(@{$job->{deps}},$contigs);
  my $target=$contigs;
  $target=~s/.contigs-disco.fa/.contigs.disco.fa/;
  push(@{$job->{target}},$target);
  #SEQTK options to clean the Disco result file
  #${SEQTK} cutN -n 1 $(patsubst %.splitN.fa,%.fa,$@)  |
  #${SEQTK} seq -L 200 - | ${SEQTK} iupac2bases - |
  #${SEQTK} rename - A | fold > $@
  my $cmd=SEQTK_BIN." cutN -n 1   $contigs | ";
     $cmd.=SEQTK_BIN." seq -L 200 -  | ";
     $cmd.=SEQTK_BIN." iupac2bases -  | ";
     $cmd.=SEQTK_BIN." rename - D | ";
     $cmd.=SEQTK_BIN." seq -l 60 - > $target";
  push(@{$job->{cmds}},$cmd);
  push(@{$self->{jobs}},$job);
}

sub has_dependency{
   return 0;
}

sub main_target{
     my $self=shift;
    return $self->{main_target};
}

1;
