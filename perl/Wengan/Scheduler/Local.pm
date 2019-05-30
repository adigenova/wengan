package Wengan::Scheduler::Local;
=head1 NAME

Wengan::Scheduler::Local

=head1 DESCRIPTION

This Module execute locally a collection of jobs.

=head2 Available methods


=cut

use strict;
use warnings;
#use Wengan::Common::GlobalConfig qw(Wengan_PROJECT_DIR Wengan_BIN_DIR);
use Module::Load;
use Data::Dumper;

my %data=();#hash storing groups of reads

sub new {
	#create an array of tools object
   	my ($class,$reads,$pipeline,%opts)=@_;
    #the default pipeline is WenganM
  if(!defined $pipeline){
	     $pipeline="WenganM";
	}
	#load the pipeline passed as argument
	load "Wengan::Pipeline::$pipeline";
		my $name="Wengan::Pipeline::$pipeline";
		$data{Pipeline} = new $name(%opts);

	#load tools of the given pipeline
  my $i=0;
  my @tools=@{$data{Pipeline}->get_tools()};
  #create the jobs and its dependeny between the pipeline tools
  foreach my $t(@tools){
    #we ask if the tool has dependecy on  a previous tool
    if($t->has_dependency()){
        #we get the main targets of the previous tool, usually the contigs sequence is the first one
        $t->set_init_dependency($tools[$i-1]->main_target());
    }
    $t->create_jobs($reads);#return the jobs associated to the current tool
		#push(@{$data{jobs}},$jobs);
    $i++;
	}
  #bless the Scheduler object
	#print Dumper($data{jobs});
	my $self = {%data,%opts,};
   	bless $self,$class;

   	return $self;
};


sub run{
     my ($self)=@_;
     $self->_create_make();
     system("make -f ".$self->{p}.".mk all"); #or die "cannot run the makefile\n";
     if ($? == -1) {
       	print "failed to execute: $!\n";
      	exit(1);
   	}
   	elsif ($? & 127) {
       	printf "child died with signal %d, %s coredump\n",
      	 ($? & 127),  ($? & 128) ? 'with' : 'without';
   	exit(1);
   	}
   	#else {
    #   		printf "child exited with value %d\n", $? >> 8;
   	#}
}

sub show_pipeline{
  my ($self)=@_;
  $self->_create_make();
  system("make -f ".$self->{p}.".mk -n all"); #or die "cannot run the makefile\n";
  if ($? == -1) {
     print "failed to execute: $!\n";
     exit(1);
 }
 elsif ($? & 127) {
     printf "child died with signal %d, %s coredump\n",
      ($? & 127),  ($? & 128) ? 'with' : 'without';
 exit(1);
 }

}


#function that create a makefile from the tools
sub _create_make{
    my ($self)=@_;
    #we use the prefix to create a file
    open(FILE,">".$self->{p}.".mk") or die "cannot create makefile\n";
    print FILE ".DELETE_ON_ERROR:\n#Wengan automatic generated makefile\n";#Instruction to delete target in an error
    my $jobs=$self->{Pipeline}->get_pipeline_jobs();
    my $last_target="";
    my $targets=();
    foreach my $job(@{$jobs}){
          print FILE _build_make_rules($job)."\n";
          $last_target=$job->{target}[0];
          push(@{$targets},@{$job->{target}});
    }
    #special rules
    print FILE "all : $last_target\n";
    print FILE "clean :\n\t-rm -f ".join(" ",@{$targets})."\n";
}

sub _build_make_rules{
    my ($job)=@_;
    #print Dumper($job);
    #is a simple job
    my $mp="";
    #is a simple job
    if(scalar(@{$job->{target}}) == 1){
        $mp.=join(" : ",$job->{target}[0], defined $job->{deps} ? join(" ",@{$job->{deps}})."\n":"\n");
        #we print the job comands
        $mp.="\t".$_."\n" foreach(@{$job->{cmds}});
    }elsif(scalar(@{$job->{target}}) > 1){
    #we create a dependency among all the jobs this is not the best solution
    #cause deleting any dependency don't rebuild the comand.
    #todo: create a rule capable of rebuild the comand when a target is deleted
    #potential solutions are discused here : https://www.gnu.org/software/automake/manual/html_node/Multiple-Outputs.html
    my @targets=@{$job->{target}};
    #we establish a dependency among the tools
    for(my $i=1; $i<scalar(@targets); $i++){
          $mp.=join(" : ",$targets[$i],$targets[$i-1])."\n";
    }
    $mp.=join(" : ",$job->{target}[0], defined $job->{deps} ? join(" ",@{$job->{deps}})."\n":"\n");
    $mp.="\t".$_."\n" foreach(@{$job->{cmds}});

    }else{
        print STDERR "The job don't produce any target\n";
       exit 1;
    }

    return $mp;
}


1;
