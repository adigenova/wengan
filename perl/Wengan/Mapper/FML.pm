
package Wengan::Mapper::FML;

=head1 NAME

Wengan::Mapper::FML

=head1 DESCRIPTION

This Object perform the task fo mapping long-read to a set of given contigs for create a spectrum of synthetic mate-pair libraries.
The synthetic mate-pair libraries are create  using the FastMin-SG progam.

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
  my $self = {contigs=>undef,dependency=>undef,cores=>$opts{t},prefix=>$opts{p}, preset=>$opts{x}, gs=>$opts{g}, uin=>$opts{i},uinccs=>$opts{I},opts=>\%opts};
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
  #we create the default jobs for 1 run of Fasmin-sg
  if($self->{preset} eq "ccspac"){
    $self->{preset}="pacraw";
    _def_create_jobs($self,$reads);
    $self->{preset}="ccspac";
    _ccs_create_jobs($self,$reads);
  }elsif($self->{preset} eq "ccsont"){
    $self->{preset}="ontraw";
    _def_create_jobs($self,$reads);
    $self->{preset}="ccsont";
    _ccs_create_jobs($self,$reads);
  }else{
    _def_create_jobs($self,$reads);
  }
}

#create default jobs for tools
sub _ccs_create_jobs{
  my ($self,$reads)=@_;

  #Ouput files from FastMin-SG
  #Modified FM for long-read
  ### sprintf(filename, "%s.I%d.fm.sam", prefix,inserts[i]);
  ## FM modified for short-Reads
  #sprintf(filename, "%s.fm.sam",prefix); #prefix should be unique
 my $job=();
  #push(@{$self->{main_target}},$self->{contigs});

  push(@{$job->{deps}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
  my $c=1;
  #push(@{$self->{main_target}},"longreads.".$self->{prefix}.".ccs2.fa"); #long read default file, we expect a single file
  #push(@{$job->{target}},"longreads.".$self->{prefix}.".ccs2.fa");
  #push(@{$job->{cmds}},join(" ","\@echo",$self->{prefix}.".ccs.ec.fa"," > ",$self->{prefix}.".fmccs.txt"));
  push(@{$job->{cmds}},join(" ","\@echo",$self->{prefix}.".ccs2",$self->{prefix}.".ccs.ec.fa"," > ",$self->{prefix}.".fmlccs.txt"));
  foreach my $i(@{$self->_get_inserts_sizes($self->{uinccs})}){
            push(@{$job->{target}},$self->{prefix}.".ccs2.I$i.fm.sam");
            push(@{$self->{main_target}},$self->{prefix}.".ccs2.I$i.fm.sam");
  }


  #FastMin-SG <preset> [options] <reference> <query>
  my $param=$self->_def_parameters(undef);
  my @mopt=($param,
  "-t ",$self->{cores},
  "-p",$self->{prefix},
  "-I ",join(",",@{$self->_get_inserts_sizes($self->{uinccs})}),
  $self->{contigs},
  $self->{prefix}.".fmlccs.txt",
  "2>".$self->{prefix}.".fmlccs.err",
  ">".$self->{prefix}.".fmlccs.log");
  push(@{$job->{cmds}},join(" ",FM_BIN,@mopt));
  push(@{$self->{jobs}},$job);
}




#create default jobs for tools
sub _def_create_jobs{
  my ($self,$reads)=@_;

  #Ouput files from FastMin-SG
  #Modified FM for long-read
  ### sprintf(filename, "%s.I%d.fm.sam", prefix,inserts[i]);
  ## FM modified for short-Reads
  #sprintf(filename, "%s.fm.sam",prefix); #prefix should be unique
 my $job=();
  push(@{$self->{main_target}},$self->{contigs});

  push(@{$job->{deps}},@{$self->{dependency}});
  my $c=1;
  push(@{$self->{main_target}},"longreads.".$self->{prefix}.$c.".fa"); #long read default file, we expect a single file
  push(@{$job->{target}},"longreads.".$self->{prefix}.$c.".fa");
  foreach my $r (@{$reads->{lreads}}){
    next if($r->{type} eq "ccs");
          if($c == 1){
            push(@{$job->{cmds}},join(" ","\@echo",$self->{prefix}.$c,$r->{long}," > ",$self->{prefix}.".fml.txt"));
          }else{
            push(@{$job->{cmds}},join(" ","\@echo",$self->{prefix}.$c,$r->{long}," >> ",$self->{prefix}.".fml.txt"));
          }
          #we save the targets that the command generate
          foreach my $i(@{$self->_get_inserts_sizes($self->{uin})}){
            push(@{$job->{target}},$self->{prefix}.$c.".I$i.fm.sam");
            push(@{$self->{main_target}},$self->{prefix}.$c.".I$i.fm.sam");
          }
          $c++;
      }

  #FastMin-SG <preset> [options] <reference> <query>
  my $param=$self->_def_parameters(undef);
  my @mopt=($param,
  "-t ",$self->{cores},
  "-p",$self->{prefix},
  "-I ",join(",",@{$self->_get_inserts_sizes($self->{uin})}),
  $self->{contigs},
  $self->{prefix}.".fml.txt",
  "2>".$self->{prefix}.".fml.err",
  ">".$self->{prefix}.".fml.log");
  push(@{$job->{cmds}},join(" ",FM_BIN,@mopt));
  push(@{$self->{jobs}},$job);
}

sub _get_inserts_sizes{
  my ($self, $inserts)=@_;
  my $is=();
  my $s="500,1000,2000,3000,4000,5000,6000,7000,8000,10000,15000,20000,30000,40000,50000";
  #FastMin-SG presets
  if($self->{preset} eq "pacccs"){
     $s="500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000";
  }elsif($self->{preset} eq "ontlon"){
     $s="500,1000,2000,3000,4000,5000,6000,7000,8000,10000,15000,20000,30000,40000,50000,60000,70000,80000,90000,100000,120000,150000,180000,200000";
  }elsif($self->{preset} eq "pacraw"){
    $s="500,1000,2000,3000,4000,5000,6000,7000,8000,10000,15000,20000";
  }elsif($self->{preset} eq "ontraw"){
    $s="500,1000,2000,3000,4000,5000,6000,7000,8000,10000,15000,20000,30000,40000,50000";
  }elsif($self->{preset} eq "ccsont" or $self->{preset} eq "ccspac"){
    $s="500,1000,2000,3000,4000,5000,6000,7000,8000,10000";#def for 10kb lib
    #$s="500,1000,2000,3000,4000,5000,6000,7000,8000,10000,12000,15000,18000";#def for 20kb lib
  }

  #user given list of insert-sizes
  if(length($inserts) > 0){
    $s=$inserts;
  }
  #we save each insert size on the lib
  push(@{$is},$_) foreach(split(",",$s));
  return $is;
}


sub _def_parameters{
      my ($self,$rlen)=@_;

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
      my  $param = "ontlon -k $k -w $w -q $q -m $m -r 300";
      #FastMin-SG presets
      if($self->{preset} eq "pacccs" or $self->{preset} eq "ccspac" or $self->{preset} eq "ccsont"){
            #preset for pacccs
            $k=21;
            $w=10;
            # we override the variables if they are given by the user
            $k = (defined $self->{opts}->{k}) ? $self->{opts}->{k}:$k;
            $q = (defined $self->{opts}->{q}) ? $self->{opts}->{q}:$q;
            $m = (defined $self->{opts}->{m}) ? $self->{opts}->{m}:$m;
            $w = (defined $self->{opts}->{w}) ? $self->{opts}->{w}:$w;

           $param = "pacccs -k $k -w $w -q $q -m $m -r 500";
           if($self->{preset} eq "ccspac" or $self->{preset} eq "ccsont"){
              #we start counting from 100M reads for the HiFi Mapping
             $param = "pacccs -k $k -w $w -q $q -m $m -r 500 -n 100000000"
           }

      }elsif($self->{preset} eq "ontlon"){
         #$param = "ontlon -k $k -w 5 -q 40 -r 300";
         $w=5;
         # we override the variables if they are given by the user
         $k = (defined $self->{opts}->{k}) ? $self->{opts}->{k}:$k;
         $q = (defined $self->{opts}->{q}) ? $self->{opts}->{q}:$q;
         $m = (defined $self->{opts}->{m}) ? $self->{opts}->{m}:$m;
         $w = (defined $self->{opts}->{w}) ? $self->{opts}->{w}:$w;
        $param = "ontlon -k $k -w $w -q $q -m $m -r 300";

      }elsif($self->{preset} eq "pacraw"){
            #$param = "pacraw -k $k -w 5 -q 40 -r 300";
            $k = (defined $self->{opts}->{k}) ? $self->{opts}->{k}:$k;
            $q = (defined $self->{opts}->{q}) ? $self->{opts}->{q}:$q;
            $m = (defined $self->{opts}->{m}) ? $self->{opts}->{m}:$m;
            $w = (defined $self->{opts}->{w}) ? $self->{opts}->{w}:$w;
           $param = "pacraw -k $k -w $w -q $q -m $m -r 300";

      }elsif($self->{preset} eq "ontraw"){
        #$param = "ontraw -k $k -w 5 -q 40 -r 300";
        $k = (defined $self->{opts}->{k}) ? $self->{opts}->{k}:$k;
        $q = (defined $self->{opts}->{q}) ? $self->{opts}->{q}:$q;
        $m = (defined $self->{opts}->{m}) ? $self->{opts}->{m}:$m;
        $w = (defined $self->{opts}->{w}) ? $self->{opts}->{w}:$w;
       $param = "ontraw -k $k -w $w -q $q -m $m -r 300";
      }
      return $param;
}

1;
