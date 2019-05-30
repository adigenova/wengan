package Wengan::Bruijn::Abyss;

use strict;
use Wengan::Common::GlobalConfig qw(ABYSS2_BIN SEQTK_BIN);

sub new{
  my ($packagename,%opts) = @_;


  #minimum variables for Abyss2
  my $self = {contigs=>undef,cores=>$opts{t},prefix=>$opts{p}, preset=>$opts{x}};
  #we ask if the contigs are passed
  if(defined $opts{c}){
    $self->{contigs}=$opts{c};
  }

  #we complain if the contigs are not given and the ABYSS2_BIN is not found
  if (!defined ABYSS2_BIN and !defined $opts{c}){
    die "ABYSS2 binary not found\n";
  }
  if(!defined SEQTK_BIN ){
      die "Seqtk binary not found for abyss post procesing\n";
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
      push(@{$job->{target}},$self->{prefix}.".abyss2-contigs.fa");
      my $cmd="ln -s $self->{contigs} $job->{target}[0]";
      push(@{$job->{cmds}},$cmd);
      push(@{$self->{jobs}},$job);
      $self->{contigs}=$job->{target}[0];
      #clean the contigs
      $self->_create_seqtk_jobs(@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
      $self->{main_target}=@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0];

  }else{
    #default method for create jobs
    $self->_create_jobs_short($reads);
  }

}

sub _def_parameters{
      my ($self,$rlen)=@_;
      #abyss configuration used on the wengan ms
      #${TIME} -v -o time_abyss_ILL150bp60X.txt abyss-pe  name=HG00733-abyss-ILL60X250
      #np=20 k=128
      #B=40G H=4 kc=3 v=-v contigs 2>abyss60X.err >abyss60X.log
      #lib="pea peb" pea="SRR5534476_1.fastq.gz SRR5534476_2.fastq.gz" peb="SRR5534475_1.fastq.gz SRR5534475_2.fastq.gz"
      my $param="k=96  B=40G H=4 kc=3 v=-v ";
      if($rlen >=250){
        $param="k=128 B=40G H=4 kc=3 v=-v ";
      }
      return $param;
}

sub _create_jobs_short{
    my ($self,$reads)=@_;
    #we get the parameters for abyss2
    my $params = $self->_def_parameters($reads->{sreads}[0]->{len});
    #preparing the reads for abyss2
    my $l=1;
    my @libs=();
    my @lnames=();
    foreach my $r (@{$reads->{sreads}}){
            my $lib="pe$l=\"$r->{fwd} $r->{rev}\"";
            push(@libs,$lib);
            push(@lnames,"pe$l");
            $l++;
    }
    #now we create the asbyss cmd
    my @mopt=("name=".$self->{prefix}.".abyss2",
      $params,
      "np=$self->{cores}",
      "lib=\"".join(" ",@lnames)."\"",
      join(" ",@libs),
      "contigs",
      "2>",$self->{prefix}.".abyss2.err",
      ">",$self->{prefix}.".abyss2.log");
      my $job=();
      #the abyss2 contigs
      push(@{$job->{target}},$self->{prefix}.".abyss2-contigs.fa");
      push(@{$job->{cmds}},join(" ",ABYSS2_BIN,@mopt));
      push(@{$self->{jobs}},$job);
    # we add the after assembly steps
    $self->_create_seqtk_jobs(@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    $self->{main_target}=@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0];

    #return $self->{jobs};
}

sub _create_seqtk_jobs{
  my ($self,$contigs)=@_;
  my $job=();
  push(@{$job->{deps}},$contigs);
  my $target=$contigs;
  $target=~s/.abyss2-contigs.fa/.abyss2.contigs.fa/;
  push(@{$job->{target}},$target);
  #SEQTK options to clean the abyss result file
  #${SEQTK} cutN -n 1 $(patsubst %.splitN.fa,%.fa,$@)  |
  #${SEQTK} seq -L 200 - | ${SEQTK} iupac2bases - |
  #${SEQTK} rename - A | fold > $@
  my $cmd=SEQTK_BIN." cutN -n 1   $contigs | ";
     $cmd.=SEQTK_BIN." seq -L 200 -  | ";
     $cmd.=SEQTK_BIN." iupac2bases -  | ";
     $cmd.=SEQTK_BIN." rename - A | ";
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
