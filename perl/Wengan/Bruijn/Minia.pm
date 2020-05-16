package Wengan::Bruijn::Minia;

use strict;
use warnings;

use Data::Dumper;
use Wengan::Common::GlobalConfig qw(MINIA3_BIN SEQTK_BIN);

sub new{
  my ($packagename,%opts) = @_;

  if (!defined MINIA3_BIN){
		die "Minia3 binary not found\n";
	}

 if(!defined SEQTK_BIN ){
      die "Seqtk binary not found for Minia3 post procesing\n";
  }
  #pipeline for ccsont
  #minimum variables for Minia3
  my $self = {contigs=>undef,cores=>$opts{t},prefix=>$opts{p},ploidy=>$opts{D},preset=>$opts{x}, main_target=>undef};
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

sub _def_parameters_ccs{
      my $param=();
      #define the multiple kmer sizes
      @{$param->{MK}}=(41,81,121,161,201,251,301,351);
      @{$param->{abun}}=(2,2,2,2,2,2,2,2);
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
      $self->_create_seqtk_jobs(@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
      $self->{contigs}=@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0];
      push(@{$self->{main_target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
      $self->_create_coveragefile_jobs($self->{contigs});
      push(@{$self->{main_target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);


  }elsif($self->{preset} eq "pacccs"){
    #create the normal jobs
    #print Dumper($reads);
    $self->_create_jobs_long($reads);

  }elsif($self->{preset} eq "ccsont" or $self->{preset} eq "ccspac"){
    if($self->{ploidy} == 1){
        #$self->_create_jobs_ccslong_haploid($reads);#haploid
        $self->_create_jobs_ccslong_diploid($reads);#diploid
      }else{
        $self->_create_jobs_ccslong_diploid($reads);#diploid
      }
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

sub _create_seqtk_jobs{
  my ($self,$contigs)=@_;
  my $job=();
  push(@{$job->{deps}},$contigs);
  my $target=$contigs;
  $target=~s/.minia.121.contigs.fa/.minia.contigs.fa/;
  $target=~s/.minia.351.contigs.fa/.minia.contigs.fa/;
  $target=~s/.minia.451.contigs.fa/.minia.contigs.fa/;
  $target=~s/.minia.251.contigs.fa/.minia.contigs.fa/;
  push(@{$job->{target}},$target);
  my $cmd.=SEQTK_BIN." seq -L 200 $contigs  | ";
     $cmd.=SEQTK_BIN." seq -l 60 - > $target";
  push(@{$job->{cmds}},$cmd);
  push(@{$self->{jobs}},$job);
}


sub _create_jobs_ccslong_diploid{
    my ($self,$reads)=@_;
    my $params = _def_parameters_ccs();

    my @hifi_files=();
    foreach my $r (@{$reads->get_hifi()}){
          push(@hifi_files,$r->{long});
    }
    my $job=();
    #we merge and transform the hifi files to fasta
    push(@{$job->{cmds}},join(" ","zcat",@hifi_files," | ",SEQTK_BIN." seq  -l 60 -A -C -  >",$self->{prefix}.".ccs.ec.fa"));
    push(@{$job->{target}},$self->{prefix}.".ccs.ec.fa");
    push(@{$self->{jobs}},$job);
     @hifi_files=();
    push(@hifi_files,$self->{prefix}.".ccs.ec.fa");

    for(my $i=0; $i<scalar(@{$params->{MK}}); $i++){
      my $job=();
      my $k=$params->{MK}[$i];
      #we iterate the reads to add the comands to file
      my $create=0;

      foreach my $r (@hifi_files){
        if($create == 0){
            push(@{$job->{cmds}},join(" ","\@echo",$r," > ",$self->{prefix}.".minia_reads.$k.txt"));
            #push(@{$job->{cmds}},join(" ","\@echo",$r->{rev}," >> ",$self->{prefix}.".minia_reads.$k.txt"));
          }else{
          push(@{$job->{cmds}},join(" ","\@echo",$r," >> ",$self->{prefix}.".minia_reads.$k.txt"));
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
        push(@{$job->{deps}},$self->{prefix}.".ccs.ec.fa");
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
	#we add commands to clean minia tmp files
	push(@{$job->{cmds}},join(" ","-rm -f","$self->{prefix}.minia.$k.unitigs.fa.glue*",
				               "$self->{prefix}.minia.$k.h5",
					       "$self->{prefix}.minia.$k.unitigs.fa"));
                 
        push(@{$self->{jobs}},$job);
    }
    # we set the main target to the contigs files
    #$self->{main_target}=@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0];
    push(@{$self->{target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    # we remove short contigs
    $self->_create_seqtk_jobs(@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    push(@{$self->{main_target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    # we add the target coverage file
    $self->_create_coveragefile_jobs(@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    #coverage file for minia
    push(@{$self->{main_target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    #$self->{main_target}=@{$self->{jobs}}[$#{$self->{jobs}}]->{target};
    #return $self->{jobs};
}


#alternative pipeline to explore the use of HiFi reads
sub _create_jobs_ccslong_diploid2{
    my ($self,$reads)=@_;
    my $params = _def_parameters();

    my @hifi_files=();
    foreach my $r (@{$reads->get_hifi()}){
          push(@hifi_files,$r->{long});
    }
    my $job=();
    #we merge and transform the hifi files to fasta
    push(@{$job->{cmds}},join(" ","zcat",@hifi_files," | ",SEQTK_BIN." seq  -l 60 -A -C -  >",$self->{prefix}.".ccs.ec.fa"));
    push(@{$job->{target}},$self->{prefix}.".ccs.ec.fa");
    push(@{$self->{jobs}},$job);
     @hifi_files=();
    push(@hifi_files,$self->{prefix}.".ccs.ec.fa");

    for(my $i=0; $i<scalar(@{$params->{MK}}); $i++){
      my $job=();
      my $k=$params->{MK}[$i];
      #we iterate the reads to add the comands to file
      my $create=0;

      foreach my $r (@hifi_files){
        if($create == 0){
            push(@{$job->{cmds}},join(" ","\@echo",$r," > ",$self->{prefix}.".minia_reads.$k.txt"));
            #push(@{$job->{cmds}},join(" ","\@echo",$r->{rev}," >> ",$self->{prefix}.".minia_reads.$k.txt"));
          }else{
          push(@{$job->{cmds}},join(" ","\@echo",$r," >> ",$self->{prefix}.".minia_reads.$k.txt"));
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
        push(@{$job->{deps}},$self->{prefix}.".ccs.ec.fa");
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
	#we add commands to clean minia tmp files
	push(@{$job->{cmds}},join(" ","-rm -f","$self->{prefix}.minia.$k.unitigs.fa.glue*",
				               "$self->{prefix}.minia.$k.h5"));
					       #"$self->{prefix}.minia.$k.unitigs.fa"));
        push(@{$self->{jobs}},$job);
    }
    # we set the main target to the contigs files
    #$self->{main_target}=@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0];
    push(@{$self->{target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    # we remove short contigs
    $self->_create_seqtk_jobs(@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    push(@{$self->{main_target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    # we add the target coverage file
    $self->_create_coveragefile_jobs(@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    #coverage file for minia
    push(@{$self->{main_target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);

    #We map the long-reads to the assembled contigs at 0.5 and 1kb

    #ran liger with 2 libs of 0.5 and 1kb
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
	#we add commands to clean minia tmp files
	push(@{$job->{cmds}},join(" ","-rm -f","$self->{prefix}.minia.$k.unitigs.fa.glue*",
				               "$self->{prefix}.minia.$k.h5",
					       "$self->{prefix}.minia.$k.unitigs.fa"));
        push(@{$self->{jobs}},$job);
    }
    #print Dumper($reads);
    # we add the target coverage file
    #$self->{main_target}=@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0];
    push(@{$self->{target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    #print Dumper($reads);
    # we remove short contigs
    $self->_create_seqtk_jobs(@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    push(@{$self->{main_target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    # we add the target coverage file
    $self->_create_coveragefile_jobs(@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    push(@{$self->{main_target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    #print Dumper($self);
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
	#we add commands to clean minia tmp files
	push(@{$job->{cmds}},join(" ","-rm -f","$self->{prefix}.minia.$k.unitigs.fa.glue*",
				               "$self->{prefix}.minia.$k.h5",
					       "$self->{prefix}.minia.$k.unitigs.fa"));
        push(@{$self->{jobs}},$job);
    }
    # we set the main target to the contigs files
    #$self->{main_target}=@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0];
    push(@{$self->{target}},@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
    # we remove short contigs
    $self->_create_seqtk_jobs(@{$self->{jobs}}[$#{$self->{jobs}}]->{target}[0]);
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
