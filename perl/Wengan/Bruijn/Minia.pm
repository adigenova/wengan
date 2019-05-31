package Wengan::Bruijn::Minia;

use strict;
use Wengan::Common::GlobalConfig qw(MINIA3_BIN);

sub new{
  my ($packagename,%opts) = @_;

  if (!defined MINIA3_BIN){
		die "Minia3 binary not found\n";
	}
  #minimum variables for Minia3
  my $self = {contigs=>undef,cores=>$opts{t},prefix=>$opts{p}, preset=>$opts{x}, main_target=>undef};
  #we ask if the contigs are passed
  if(defined $opts{c}){
    $self->{contigs}=$opts{c};
  }
  bless ($self, $packagename);
  return ($self);
}

#method to mark that contigs were already assembled previously
sub contigs_done{
      my ($self,$file)=@_;
      $self->{contigs}=$file;
}

sub _def_parameters{
      my $param=();
      #define the multiple kmer sizes
      @{$param->{MK}}=(41,81,121);
      @{$param->{abun}}=(2,2,2);
      return $param;
}

#generic function is called from the WenganM pipeline
sub create_jobs{
  my ($self,$reads)=@_;
  #minia contigs were create previously
  if(defined $self->{contigs}){
      my $job=();
      push(@{$job->{target}},$self->{prefix}.".minia.121.contigs.fa");
      my $cmd="ln -s $self->{contigs} $job->{target}[0]";
      push(@{$job->{cmds}},$cmd);
      push(@{$self->{jobs}},$job);
      $self->{contigs}=$job->{target}[0];
      #$self->{main_target}=@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0];
      push(@{$self->{main_target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
      $self->_create_coveragefile_jobs($self->{contigs});
      push(@{$self->{main_target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);

  }elsif($self->{preset} eq "pacccs"){
    #create the normal jobs
    $self->_create_jobs_long($reads);
  }else{
    #default method for create jobs
    $self->_create_jobs_short($reads);
  }

}

sub _create_coveragefile_jobs{
  my ($self,$contigs)=@_;
  my $job=();
  push(@{$job->{deps}},$contigs);
  my $target=$contigs;
  $target=~s/.fa/.cov.txt/;
  push(@{$job->{target}},$target);
  ##grep ">" $1 | sed 's/km:f://' | awk '{print $1" "$4}' | sed 's/>//g' > $1.cov
  my $cmd='grep ">" '.$contigs.' | sed \'s/km:f://\' | awk \'{print $$1" "$$4}\' | sed \'s/>//g\' > '.$target;
  push(@{$job->{cmds}},$cmd);
  push(@{$self->{jobs}},$job);
}

sub _create_jobs_long{
    my ($self,$reads)=@_;
    my $params = _def_parameters();
    for(my $i=0; $i<scalar(@{$params->{MK}}); $i++){
      my $job=();
      my $k=$params->{MK}[$i];
      #we iterate the reads to add the comands to file
      my $create=0;
      foreach my $r (@{$reads->{lreads}}){
        if($create == 0){
            push(@{$job->{cmds}},join(" ","\@echo",$r->{long}," > ",$self->{prefix}.".minia_reads.$k.txt"));
            #push(@{$job->{cmds}},join(" ","\@echo",$r->{rev}," >> ",$self->{prefix}.".minia_reads.$k.txt"));
          }else{
          push(@{$job->{cmds}},join(" ","\@echo",$r->{long}," >> ",$self->{prefix}.".minia_reads.$k.txt"));
          #push(@{$job->{cmds}},join(" ","\@echo",$r->{rev}," >> ",$self->{prefix}.".minia_reads.$k.txt"));
        }
        $create++;
      }
      #we have to add the previous contigs
      if($i > 0){
        my $contigs = $self->{prefix}.".minia.".$params->{MK}[$i-1].".contigs.fa";
        #we add three times the minia contigs of the previous kmer
        push(@{$job->{cmds}},join(" ","\@echo",$contigs," >> ",$self->{prefix}.".minia_reads.$k.txt"));
        push(@{$job->{cmds}},join(" ","\@echo",$contigs," >> ",$self->{prefix}.".minia_reads.$k.txt"));
        push(@{$job->{cmds}},join(" ","\@echo",$contigs," >> ",$self->{prefix}.".minia_reads.$k.txt"));
        push(@{$job->{target}},$self->{prefix}.".minia.".$params->{MK}[$i].".contigs.fa");
        push(@{$job->{deps}},$self->{prefix}.".minia.".$params->{MK}[$i-1].".contigs.fa");
      }else{
        push(@{$job->{target}},$self->{prefix}.".minia.".$params->{MK}[$i].".contigs.fa");
      }
      #now we create the minia cmd
      my @mopt=("-in",$self->{prefix}.".minia_reads.$k.txt",
        "-kmer-size",$params->{MK}[$i],
        "-abundance-min",$params->{abun}[$i],
        "-out",$self->{prefix}.".minia.".$params->{MK}[$i],
        "-minimizer-size",10,
        "-max-memory",5000,
        "-nb-cores",$self->{cores},
        "2>",$self->{prefix}.".minia.".$params->{MK}[$i].".err",
        ">",$self->{prefix}.".minia.".$params->{MK}[$i].".log");
        push(@{$job->{cmds}},join(" ",MINIA3_BIN,@mopt));
        push(@{$self->{jobs}},$job);
    }
    # we add the target coverage file
    #$self->{main_target}=@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0];
    push(@{$self->{main_target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    # we add the target coverage file
    $self->_create_coveragefile_jobs(@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    push(@{$self->{main_target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    #return $self->{jobs};
}




sub _create_jobs_short{
    my ($self,$reads)=@_;
    my $params = _def_parameters();
    for(my $i=0; $i<scalar(@{$params->{MK}}); $i++){
      my $job=();
      my $k=$params->{MK}[$i];
      #we iterate the reads to add the comands to file
      my $create=0;
      foreach my $r (@{$reads->{sreads}}){
        if($create == 0){
            push(@{$job->{cmds}},join(" ","\@echo",$r->{fwd}," > ",$self->{prefix}.".minia_reads.$k.txt"));
            push(@{$job->{cmds}},join(" ","\@echo",$r->{rev}," >> ",$self->{prefix}.".minia_reads.$k.txt"));
          }else{
          push(@{$job->{cmds}},join(" ","\@echo",$r->{fwd}," >> ",$self->{prefix}.".minia_reads.$k.txt"));
          push(@{$job->{cmds}},join(" ","\@echo",$r->{rev}," >> ",$self->{prefix}.".minia_reads.$k.txt"));
        }
        $create++;
      }
      #we have to add the previous contigs
      if($i > 0){
        my $contigs = $self->{prefix}.".minia.".$params->{MK}[$i-1].".contigs.fa";
        #we add three times the minia contigs of the previous kmer
        push(@{$job->{cmds}},join(" ","\@echo",$contigs," >> ",$self->{prefix}.".minia_reads.$k.txt"));
        push(@{$job->{cmds}},join(" ","\@echo",$contigs," >> ",$self->{prefix}.".minia_reads.$k.txt"));
        push(@{$job->{cmds}},join(" ","\@echo",$contigs," >> ",$self->{prefix}.".minia_reads.$k.txt"));
        push(@{$job->{target}},$self->{prefix}.".minia.".$params->{MK}[$i].".contigs.fa");
        push(@{$job->{deps}},$self->{prefix}.".minia.".$params->{MK}[$i-1].".contigs.fa");
      }else{
        push(@{$job->{target}},$self->{prefix}.".minia.".$params->{MK}[$i].".contigs.fa");
      }
      #now we create the minia cmd
      my @mopt=("-in",$self->{prefix}.".minia_reads.$k.txt",
        "-kmer-size",$params->{MK}[$i],
        "-abundance-min",$params->{abun}[$i],
        "-out",$self->{prefix}.".minia.".$params->{MK}[$i],
        "-minimizer-size",10,
        "-max-memory",5000,
        "-nb-cores",$self->{cores},
        "2>",$self->{prefix}.".minia.".$params->{MK}[$i].".err",
        ">",$self->{prefix}.".minia.".$params->{MK}[$i].".log");
        push(@{$job->{cmds}},join(" ",MINIA3_BIN,@mopt));
        push(@{$self->{jobs}},$job);
    }
    # we set the main target to the contigs files
    #$self->{main_target}=@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0];
    push(@{$self->{main_target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    # we add the target coverage file
    $self->_create_coveragefile_jobs(@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    #coverage file for minia
    push(@{$self->{main_target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    #$self->{main_target}=@{$self->{jobs}}[$#{$self->{jobs}}]->{target};
    #return $self->{jobs};
}

sub has_dependency{
   return 0;
}

sub main_target{
     my $self=shift;
    return $self->{main_target};
}

1; #EOM
