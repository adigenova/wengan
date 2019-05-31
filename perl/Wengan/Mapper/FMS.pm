
package Wengan::Mapper::FMS;

=head1 NAME

Wengan::Mapper::FMS

=head1 DESCRIPTION

This Object perform the task fo mapping short-read to a set of given contigs.
The short-reads are mapped unsing the FastMin-SG progam.

=head2 Available methods


=cut

use strict;
use Wengan::Common::GlobalConfig qw(FM_BIN);

sub new{
  my ($packagename,%opts) = @_;

  if (!defined FM_BIN ){
		die "FastMin-SG binary not found\n";
	}

  #minimum variables for FastMin-SG
  my $self = {contigs=>undef,dependency=>undef,cores=>$opts{t},prefix=>$opts{p}, preset=>$opts{x}};
  #we ask if the contigs are passed
  if(defined $opts{c}){
    $self->{contigs}=$opts{c};
  }
  bless ($self, $packagename);
  return ($self);
}

sub set_init_dependency{
      my ($self,$dep)=@_;
      push(@{$self->{dependency}},@{$dep});
      $self->{contigs}=@$dep[0];
}

sub has_dependency{
   return 1;
}

sub main_target{
     my $self=shift;
    return $self->{main_target};
}


#generic function is called from the WenganM pipeline
sub create_jobs{
  my ($self,$reads)=@_;

  #Ouput files from FastMin-SG
  #Modified FM for long-read
  ### sprintf(filename, "%s.I%d.fm.sam", prefix,inserts[i]);
  ## FM modified for short-Reads
  #sprintf(filename, "%s.fm.sam",prefix); #prefix should be unique
 my $job=();
  push(@{$self->{main_target}},$self->{contigs});
  if($self->{preset} eq "pacccs"){
    #pacbio CCS reads, non-hybrid assembly
      push(@{$job->{deps}},@{$self->{dependency}});
      my $c=1;
      foreach my $r (@{$reads->{lreads}}){
          if($c == 1){
            push(@{$job->{cmds}},join(" ","\@echo",$self->{prefix}.$c.".ccs",$r->{long}," > ",$self->{prefix}.".fms.txt"));
          }else{
            push(@{$job->{cmds}},join(" ","\@echo",$self->{prefix}.$c,$r->{long}," >> ",$self->{prefix}.".fms.txt"));
          }
          #we save the targets that the command generate
          push(@{$job->{target}},$self->{prefix}.$c.".ccs.I500.fm.sam");
          push(@{$self->{main_target}},$self->{prefix}.$c.".ccs.I500.fm.sam");
          $c++;
      }
      #now we create the comand that run FastMin-SG
  }else{
    #normal short reads
    push(@{$job->{deps}},@{$self->{dependency}});
    my $c=1;
    foreach my $r (@{$reads->{sreads}}){
        if($c == 1){
          push(@{$job->{cmds}},join(" ","\@echo",$self->{prefix}.$c,$r->{fwd},$r->{rev}," > ",$self->{prefix}.".fms.txt"));
        }else{
          push(@{$job->{cmds}},join(" ","\@echo",$self->{prefix}.$c,$r->{fwd},$r->{rev}," >> ",$self->{prefix}.".fms.txt"));
        }
        #we save the targets that the command generate
        push(@{$job->{target}},$self->{prefix}.$c.".fm.sam");
        push(@{$self->{main_target}},$self->{prefix}.$c.".fm.sam");
        $c++;
    }
  }

  #FastMin-SG <preset> [options] <reference> <query>
  my $param=$self->_def_parameters($reads->get_read_length());
  my @mopt=($param,
  "-t",$self->{cores},
  $self->{contigs},
  $self->{prefix}.".fms.txt",
  "2>".$self->{prefix}.".fms.err",
  ">".$self->{prefix}.".fms.log");
  push(@{$job->{cmds}},join(" ",FM_BIN,@mopt));
  push(@{$self->{jobs}},$job);
}

sub _def_parameters{
      my ($self,$rlen)=@_;
      #shortr -t $(cpu) $(srf) $(ctg)
      # 250bp reads
      my $param = "shortr -c 30 -k 21 -w 10 -q 20 -r 50000";
      if($rlen <=200){
        #150 bp reads
          $param = "shortr -c 50 -k 21 -w 10 -q 20 -r 50000";
      }
      #preset for pacccs
      if($self->{preset} eq "pacccs"){
            $param = "pacccs -k 21 -w 10 -q 30 -r 500 -I 500";
      }
      return $param;
}

1;
