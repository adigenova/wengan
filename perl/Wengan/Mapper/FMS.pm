
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
  my $self = {contigs=>undef,dependency=>undef,cores=>$opts{t},prefix=>$opts{p}, preset=>$opts{x},gs=>$opts{g},pipeline=>$opts{a}};
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
  }elsif($self->{preset} eq "ccsont" or $self->{preset} eq "ccspac"){
          my $c=1;
          push(@{$job->{deps}},@{$self->{dependency}});
          push(@{$job->{cmds}},join(" ","\@echo",$self->{prefix}.$c.".ccs",$self->{prefix}.".ccs.ec.fa"," > ",$self->{prefix}.".fms.txt"));
          push(@{$job->{target}},$self->{prefix}.$c.".ccs.I500.fm.sam");
          push(@{$job->{target}},$self->{prefix}.$c.".ccs.I1000.fm.sam");
          push(@{$job->{target}},$self->{prefix}.$c.".ccs.I2000.fm.sam");
          push(@{$self->{main_target}},$self->{prefix}.$c.".ccs.I500.fm.sam");
          push(@{$self->{main_target}},$self->{prefix}.$c.".ccs.I1000.fm.sam");
          push(@{$self->{main_target}},$self->{prefix}.$c.".ccs.I2000.fm.sam");
          #FastMin-SG <preset> [options] <reference> <query>
          my $param=$self->_def_parameters_long();
          my @mopt=($param,
          "-t",$self->{cores},
          $self->{contigs},
          $self->{prefix}.".fms.txt",
          "2>".$self->{prefix}.".fms.err",
          ">".$self->{prefix}.".fms.log");
          push(@{$job->{cmds}},join(" ",FM_BIN,@mopt));
          if($self->{preset} eq "ccsont" or $self->{preset} eq "ccspac"){
            my $c=1;
            #delete the long-read created from fasmin-sg
            push(@{$job->{cmds}},join(" ","-rm -f","longreads.".$self->{prefix}.$c.".ccs.fa"));
          }
          push(@{$self->{jobs}},$job);
  }else{
    #normal short reads
    #push(@{$job->{deps}},@{$self->{dependency}}[0]);
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


    if($self->{pipeline} ne "M"){
    #we add the jobs for mapping mate pairs from long-reads
    $job=();#job for mp lib
    my $c=1;
    push(@{$job->{deps}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    foreach my $r (@{$reads->{lreads}}){
      next if($r->{type} eq "ccs");
            if($c == 1){
              push(@{$job->{cmds}},join(" ","\@echo",$self->{prefix}.".im.".$c,$r->{long}," > ",$self->{prefix}.".fml.im.txt"));
            }else{
              push(@{$job->{cmds}},join(" ","\@echo",$self->{prefix}.".im.".$c,$r->{long}," >> ",$self->{prefix}.".fml.im.txt"));
            }
            #we save the targets that the command generate
            my @ins=(500,1000,2000);
            foreach my $i(@ins){
              push(@{$job->{target}},$self->{prefix}.".im.".$c.".I$i.fm.sam");
              push(@{$self->{main_target}},$self->{prefix}.".im.".$c.".I$i.fm.sam");
            }
            $c++;
    }

    $param=$self->_def_parameters_long();
    @mopt=();
    @mopt=($param,
    "-t",$self->{cores},
    $self->{contigs},
    $self->{prefix}.".fml.im.txt",
    "2>".$self->{prefix}.".fml.im.err",
    ">".$self->{prefix}.".fml.im.log");
    push(@{$job->{cmds}},join(" ",FM_BIN,@mopt));
    #delete the long-read created from fasmin-sg
    $c--;
    push(@{$job->{cmds}},join(" ","-rm -f","longreads.".$self->{prefix}.".im.".$c.".fa"));
    push(@{$self->{jobs}},$job);
   }
  }


}
#def params for long-reads to create 500,1kb and 2kb libraries
sub _def_parameters_long{
  my ($self)=@_;
  #larger genomes
  my $k=20;#default k-mer size for >0.5Gb genomes
  #mid size genomes
  $k=19 if($self->{gs} <500);
  #bacterial genomes
  $k=15 if($self->{gs} <10);

  my $w=5;
  my $q=40;
  my $m=150;

  #we define the preset for mapping the long-reads to the set of contigs
  my  $param = "ontlon -k $k -w $w -q $q -m $m -r 300 -I 500,1000,2000";

  #preset for pacccs
  if($self->{preset} eq "pacccs"){
        $param = "pacccs -k 21 -w 10 -q 30 -r 500 -I 500,1000,2000 ";
  }
  #preset for ccsont or ccspac
  if($self->{preset} eq "ccsont" or $self->{preset} eq "ccspac"){
        $param = "pacccs -k 21 -w 10 -q 30 -r 500 -I 500,1000,2000 ";

  }elsif($self->{preset} eq "ontlon"){
       #$param = "ontlon -k $k -w 5 -q 40 -r 300";
       $w=5;
       # we override the variables if they are given by the user
       $k = (defined $self->{opts}->{k}) ? $self->{opts}->{k}:$k;
       $q = (defined $self->{opts}->{q}) ? $self->{opts}->{q}:$q;
       $m = (defined $self->{opts}->{m}) ? $self->{opts}->{m}:$m;
       $w = (defined $self->{opts}->{w}) ? $self->{opts}->{w}:$w;
      $param = "ontlon -k $k -w $w -q $q -m $m -r 300 -I 500,1000,2000 ";

    }elsif($self->{preset} eq "pacraw"){
          #$param = "pacraw -k $k -w 5 -q 40 -r 300";
          $k = (defined $self->{opts}->{k}) ? $self->{opts}->{k}:$k;
          $q = (defined $self->{opts}->{q}) ? $self->{opts}->{q}:$q;
          $m = (defined $self->{opts}->{m}) ? $self->{opts}->{m}:$m;
          $w = (defined $self->{opts}->{w}) ? $self->{opts}->{w}:$w;
         $param = "pacraw -k $k -w $w -q $q -m $m -r 300 -I 500,1000,2000 ";

    }elsif($self->{preset} eq "ontraw"){
      #$param = "ontraw -k $k -w 5 -q 40 -r 300";
      $k = (defined $self->{opts}->{k}) ? $self->{opts}->{k}:$k;
      $q = (defined $self->{opts}->{q}) ? $self->{opts}->{q}:$q;
      $m = (defined $self->{opts}->{m}) ? $self->{opts}->{m}:$m;
      $w = (defined $self->{opts}->{w}) ? $self->{opts}->{w}:$w;
     $param = "ontraw -k $k -w $w -q $q -m $m -r 300 -I 500,1000,2000 ";
    }
  return $param;
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
      return $param;
}

1;
