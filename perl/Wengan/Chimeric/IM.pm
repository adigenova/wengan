package Wengan::Chimeric::IM;

=head1 NAME

Wengan::Chimeric::IM

=head1 DESCRIPTION

This Object perform the task of detecting assembly errors on a set of contigs.
The assembly is validated using the IntervalMiss progam.

=cut

use strict;
use Wengan::Common::GlobalConfig qw(IM_BIN);

sub new{
  my ($packagename,%opts) = @_;

  if (!defined IM_BIN ){
		die "IntervalMiss binary not found\n";
	}

  #minimum variables for IntervalMiss
  my $self = {contigs=>undef,dependency=>undef,cores=>$opts{t},prefix=>$opts{p}, preset=>$opts{x}, pipeline=>$opts{a},opts=>\%opts};
  #we ask if the contigs are passed
  if(defined $opts{c}){
    $self->{contigs}=$opts{c};
  }
  bless ($self, $packagename);
  return ($self);
}

#methods for controling the dependencies between tools
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


  my $job=();
  push(@{$job->{deps}},@{$self->{dependency}});
    my $c=1;
  foreach my $d (@{$self->{dependency}}){
      #we skypt the contigs sequences
      next if($d eq $self->{contigs});
      my $t=0;#short
      $t=1 if($d=~m/.im./);
      if($c == 1){
            if($self->{pipeline} eq "D" or $self->{pipeline} eq "A"){
                  push(@{$job->{cmds}},join(" ","\@echo \"$d\t$t\""," > ",$self->{prefix}.".fms.sams.txt"));
            }else{
                  push(@{$job->{cmds}},join(" ","\@echo",$d," > ",$self->{prefix}.".fms.sams.txt"));
            }
      }else{
        if($self->{pipeline} eq "D" or $self->{pipeline} eq "A"){
              push(@{$job->{cmds}},join(" ","\@echo \"$d\t$t\""," >> ",$self->{prefix}.".fms.sams.txt"));
        }else{
            push(@{$job->{cmds}},join(" ","\@echo",$d," >> ",$self->{prefix}.".fms.sams.txt"));
        }
      }
      $c++;
  }
  #outputs from intervalmiss
  # string p = prefix + ".MBC" + to_string(ai.mbc_arg[i]);
  #//we create the fasta of the contigs
  #ofstream fasta; //save normal ctgs and quimeric splitted
  #fasta.open(prefix+".msplit.fa");
  my $param=$self->_def_parameters(undef);
  my @mopt=($param,
  "-t",$self->{cores},
  "-s",$self->{prefix}.".fms.sams.txt",
  "-c ",$self->{contigs},
  "-p ",$self->{prefix},
  "2>".$self->{prefix}.".im.err",
  ">".$self->{prefix}.".im.log");
  push(@{$job->{cmds}},join(" ",IM_BIN,@mopt));
  #we save the targets that the command generate
  push(@{$job->{target}},$self->{prefix}.".MBC".$self->_get_depth().".msplit.fa");
  push(@{$self->{jobs}},$job);
  $self->_create_coveragefile_jobs($self->{prefix}.".MBC".$self->_get_depth().".msplit.fa");
  #main results from intervalMiss
  push(@{$self->{main_target}},$self->{prefix}.".MBC".$self->_get_depth().".msplit.fa");
  push(@{$self->{main_target}},$self->{prefix}.".MBC".$self->_get_depth().".msplit.cov.txt");

}

sub _create_coveragefile_jobs{
  my ($self,$contigs)=@_;
  my $job=();
  push(@{$job->{deps}},$contigs);
  my $target=$contigs;
  $target=~s/.fa/.cov.txt/;
  push(@{$job->{target}},$target);
  my $cmd='grep ">" '.$contigs.' | sed \'s/>//\' | awk \'{print $$1" "$$2}\' | sed \'s/>//g\' > '.$target;
  push(@{$job->{cmds}},$cmd);
  push(@{$self->{jobs}},$job);
}

sub _get_depth{
    my $self=shift;
  return $self->{d};
}

sub _set_depth{
  my $self=shift;
  my $d=shift;
  $self->{d}=$d;
}
#some parameters for WenganM pipeline
sub _def_parameters{
      my ($self,$rlen)=@_;

      #default depth
      my $d = 7;
      my $param = "-d $d";
      $self->_set_depth($d);
      #we adjust for PACCCS data
      if($self->{preset} eq "pacccs"){
            $d=3;
            $param = "-d $d";
            $self->_set_depth($d);
      }
      #Discovar contigs are more accurate than Minia3 and Abyss2.
      if($self->{pipeline} eq "D"){
          $d = 1;
          $param = "-d $d --clib 1 ";
          $self->_set_depth($d);
      }

      if($self->{pipeline} eq "A"){
          $param = "-d $d --clib 1 ";
          $self->_set_depth($d);
      }
      #we ask if the -d option was specified
      if(defined $self->{opts}->{d}){
         $d = $self->{opts}->{d};
        $param = "-d $d --clib 1";#affect only A and D
        $self->_set_depth($d);
      }
      #means WenganM pipeline, we have to provide the coverage minia3 file
      if($self->{pipeline} eq "M"){
          my $fst=0.1;
          if (defined $self->{opts}->{f}){
            $fst=$self->{opts}->{f};
          }
          $param.=" --fst $fst -b ".$self->{prefix}.".minia.contigs.cov.txt";
      }
      return $param;
}






1;
