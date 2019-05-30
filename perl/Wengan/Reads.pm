package Wengan::Reads;


=head1 NAME

Wengan::Reads

=head1 DESCRIPTION

This Object handles the reads given as input for assembly.
It check that all of them exists as well as guess some parameters and dumb user errors.


=head2 Available methods

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=over 5

=item new()

=item add_long_reads()

=item add_long()

=item add_short_reads()

=item add_pair()

=item is_short()

=item get_read_length()

=item _get_dir()

=item _guess_read_length()

=back

=cut

use strict;

sub new{
  my ($packagename) = @_;

  my $self = { sreads => [], lreads=> [] };

  bless ($self, $packagename);
  return ($self);
}

sub add_long_reads{
	my ($self,$files)=@_;
	chomp $files;
	my @files=split(",",$files);
	#we add the short-reads files
	for(my $i=0; $i<=$#files; $i++){
		$self->add_long($files[$i]);
	}
}

sub add_long{
	my ($self,$file)=@_;
	my $tmp=();
	$tmp->{long}=$file;
	$tmp->{type}="long";
	$tmp->{path}=_get_dir($file);
	$tmp->{len}=_guess_read_length($file);
	#we store the pair end information
	push(@{$self->{lreads}},$tmp);
}


sub add_short_reads{
	my ($self,$files)=@_;
	chomp $files;
	my @files=split(",",$files);
	#we add the short-reads files
	for(my $i=0; $i<$#files; $i+=2){
		$self->add_pair($files[$i],$files[$i+1]);
	}
}


sub add_pair{
	my ($self,$fwd,$rev)=@_;
	if($fwd eq $rev){
		print $fwd." ".$rev." are the same\n";
		exit 1;
	}

	my $tmp=();
	$tmp->{fwd}=$fwd;
	$tmp->{rev}=$rev;
	$tmp->{type}="short";
	$tmp->{path}=_get_dir($fwd);
	$tmp->{len}=_guess_read_length($fwd);
	#we store the pair end information
	push(@{$self->{sreads}},$tmp);
}

#we wait a fastq file compressed with gzip
sub _guess_read_length{
	my ($file)=@_;

	open(FILE,"gzip -dc ".$file." | ") or die "cannot open file $file\n";
	my $readc=0;
	my $l=0;
	#we assume fastq format
	while(my $line=<FILE>){
		$l+=length(<FILE>);
		#we skyp the next 2 line
		<FILE>;<FILE>;
		if($readc >= 1000){
			last;
		}
	        $readc++;
	}
	close(FILE);
	return int($l/$readc);
}

sub is_short{
    my $self=shift;
    if(defined $self->{sreads}){
        return 1;
    }else{
      return 0;
    }
}

sub get_read_length{
    my $self=shift;
  if($self->is_short()){
      return $self->{sreads}[0]->{len};
  }else{
    return  250;
  }
}
#get dirname
sub _get_dir{
	my $prog = shift;
	return '.' unless ($prog =~ /\//);
	$prog =~ s/\/[^\s\/]+$//g;
	return $prog;
}

1; #EOM
