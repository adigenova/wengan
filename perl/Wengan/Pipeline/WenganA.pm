package Wengan::Pipeline::WenganA;

=head1 NAME

Wengan::Pipeline::WenganA

=head1 DESCRIPTION

This Module load the WenganA pipeline implemented by the Wengan assembler.

=head2 Available methods


=cut

use strict;
use warnings;
use Module::Load;
#use Wengan::Common::GlobalConfig qw(Wengan_PROJECT_DIR);
#use Data::Dumper;

my @tools_name=("Wengan::Bruijn::Abyss","Wengan::Mapper::FMS","Wengan::Chimeric::IM","Wengan::Mapper::FML","Wengan::SSG::LG");

my %data=(
	tools=>[],
);

sub new {
	#create an array of tools object
   	my ($class,%opts)=@_;
	my $tools=();
	foreach my $t(@tools_name){
		load $t;
		my $name=$t;
		my $tmp= new $name(%opts);
    #if($t eq "TopHat"){
		 #$tmp->set_gtf(1); #set a parameter for this tool
		#}
		push(@{$tools},$tmp);
	}


	$data{tools}=$tools;#general tools
	my $self = {%data,};
   	bless $self,$class;
	return $self;
};

#return the read level tools
sub get_tools{
	my ($self)=@_;
	return $self->{tools};
}

#return all the jobs that have to be executed by the pipeline
sub get_pipeline_jobs{
    my ($self)=shift;
    my $jobs=();
  foreach my $t(@{$self->{tools}}){
      push(@{$jobs},@{$t->{jobs}});
  }
  return $jobs;
}

#return some specials tools
sub get_tools_once{
	my ($self)=@_;
	return $self->{tools_once};
}

1;
