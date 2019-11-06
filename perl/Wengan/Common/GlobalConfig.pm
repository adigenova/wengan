package Wengan::Common::GlobalConfig;

=head1 NAME

Wengan::Common::GlobalConfig

=head1 DESCRIPTION

Global configuration settings used for Wengan assembler


=head2 Binaries Configuration

=over 4

=item * seqtk

The seqtk  tools is use to FASTA/FASTQ files preprocessing
Here we define the binary locations.

=item *  Minia3

Minia3 is a ultrafast and low memory short-read assembler.

=item * Abyss2

Abyss2 is a low memory short-read assembler.

=item * DiscoVarDenovo

DiscoVarDenovo is a high memory short-read assembler that exploit the pair-end information.

=item * FastMin-SG

Fastmin-SG is an ultrafast alignment-free algorihtm for pseudo-alignment fo pair-end reads from short or long-reads.

=item * IntervalMiss

Identify and Error correct short-read contigs from  pair-end pseudo-alignments.

=item * Liger

Liger build and operate over the Synthetic Scaffolding Graph.

=item * Additional perl scripts

Set of perl scripts that performs diferents task for the correct work of Wengan.

=back

=cut

use constant SEQTK_BIN => $ENV{WROOTDIR}.'/bin/seqtk';
use constant MINIA3_BIN =>$ENV{WROOTDIR}.'/bin/minia'; #put undef for stop a execution
use constant ABYSS2_BIN => $ENV{WROOTDIR}.'/bin/abyss-pe';
use constant DiscoVarDenovo_BIN => $ENV{WROOTDIR}.'/bin/DiscovarExp';
use constant FM_BIN => $ENV{WROOTDIR}.'/bin/fastmin-sg';
use constant IM_BIN => $ENV{WROOTDIR}.'/bin/intervalmiss';
use constant LIGER_BIN => $ENV{WROOTDIR}.'/bin/liger';
use constant AUXPERL_DIR => $ENV{WROOTDIR}.'/aux_scripts/';

use strict;
use warnings;
use base qw(Exporter);
our $VERSION = 1.0;

our @EXPORT_OK=qw(
      MINIA3_BIN
      ABYSS2_BIN
      DiscoVarDenovo_BIN
      FM_BIN
      IM_BIN
      LIGER_BIN
      AUXPERL_DIR
      SEQTK_BIN
);
1;
