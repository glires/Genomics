#!/usr/bin/env perl

=head1 NAME

  ns0010 - check paired-end FASTQ data and reprint

=head1 SYNOPSIS

  ns0010 -h
  ns0010 -a sample_R1_0.fastq -b sample_R2_0.fastq \
         -c sample_R1_1.fastq -d sample_R2_1.fastq \
         -e sample_R3_1.fastq[-x]

=head1 DESCRIPTION

  This script reads a couple of paired-end FASTQ data files
  and check the format and information.
  It checks the following.
    - One FASTQ entry should consist of four lines.
    - One sequence should consists of T, C, A, G, or N.
    - Lengths of a read and its scores should be identical.
    - All reads should have the same lengtgh.
    - All entries should be paired.
  CASAVA-filtered reads which have Y in $8 are not printed.
  Unpaired reads are printed onto the third output file,
  which will be later deleted if the size is zero.
  For single-end read data, an identical file name should be
  specified with options -a and -b.
  So are options -c and -d.
  Do not forget to add option -s for single-end data.

=head1 OPTIONS

  -h  Print usage
  -a  Input file 1
  -b  Input file 2
  -c  Output file 1
  -d  Output file 2
  -e  Output file 3 (single end data)
  -s  Single-end read data
  -x  Avoid read length check

=head1 EXAMPLES

  $ perldoc ns0010
  $ ns0010 -h
  $ ns0010 -a sample_R1_0.fastq -b sample_R2_0.fastq \
           -c sample_R1_1.fastq -d sample_R2_1.fastq \
           -e sample_R3_1.fastq

=head1 AUTHOR

  Coded by Kohji OKAMURA, Ph.D.

=head1 HISTORY

  Jun 07, 2012  The first version from scratch
  Jun 08, 2012  POD finished
  Jan 31, 2013  Renamed from jb0001.pl to jf0010.pl
  Jan 31, 2013  Flowcell ID can contain hyphens
  Jul 03, 2013  Instrument ID can contain hyphens (HiSeq 2500)
  Jul 29, 2013  Renamed to ns0010.pl; delete $opts{e} if zero
  Oct 11, 2013  Renamed to ns0010; called from ns0001
  Nov 21, 2013  Overlook 1-bp difference in sequence lengths
  Dec 18, 2014  Add option -x to avid length check
  Mar 04, 2015  Add option -s for single-end data

=cut


use warnings;
use strict;
use Getopt::Std;

my %opts;
my ($line11, $line12, $line13, $line14);
my ($line21, $line22, $line23, $line24);
my ($lane1, $tile1, $x1, $y1, $r1, $filtered1);
my ($lane2, $tile2, $x2, $y2, $r2, $filtered2);
my $len = 0;	# this initialized value should be replaced with the length

getopts 'ha:b:c:d:e:sx', \%opts;

if (defined $opts{h})                                     { exec "perldoc $0" }
unless (defined $opts{a})                       { die "Option -a is required" }
unless (defined $opts{b})                       { die "Option -b is required" }
unless (defined $opts{c})                       { die "Option -c is required" }
unless (defined $opts{d})                       { die "Option -d is required" }
unless (defined $opts{e})                       { die "Option -e is required" }

open INPUT1, $opts{a}                                    or die "$!: $opts{a}";
open INPUT2, $opts{b}                                    or die "$!: $opts{b}";
open OUTPUT1, "> $opts{c}"                               or die "$!: $opts{c}";
open OUTPUT2, "> $opts{d}"                               or die "$!: $opts{d}";
open OUTPUT3, "> $opts{e}"                               or die "$!: $opts{e}";

while ($line11 = <INPUT1>)
{
  $line12 = <INPUT1>;
  $line13 = <INPUT1>;
  $line14 = <INPUT1>;
  $line21 = <INPUT2>;
  $line22 = <INPUT2>;
  $line23 = <INPUT2>;
  $line24 = <INPUT2>;

  if ($line11 =~ /^\@[\w\-]+:\d+:[\w\-]+:(\d+):(\d+):(\d+):(\d+) (\d):([YN]):/)
  { ($lane1, $tile1, $x1, $y1, $r1, $filtered1) = ($1, $2, $3, $4, $5, $6) }
  else                                        { die "Format error 1: $line11" }
  if ($line21 =~ /^\@[\w\-]+:\d+:[\w\-]+:(\d+):(\d+):(\d+):(\d+) (\d):([YN]):/)
  { ($lane2, $tile2, $x2, $y2, $r2, $filtered2) = ($1, $2, $3, $4, $5, $6) }
  else                                        { die "Format error 2: $line21" }

  if ($filtered1 eq 'Y' and $filtered2 eq 'Y') { next }

  if ($filtered1 eq 'Y')
  { output_single($line21, $line22, $line23, $line24); next }

  if ($filtered2 eq 'Y')
  { output_single($line11, $line12, $line13, $line14); next }

  if ($filtered1 eq 'N' and $filtered2 eq 'N')
  {
    unless ($lane1 == $lane2)      { die "Paired error 1: \n$line11\n$line21" }
    unless ($tile1 == $tile2)      { die "Paired error 2: \n$line11\n$line21" }
    unless ($x1 == $x2)            { die "Paired error 3: \n$line11\n$line21" }
    unless ($y1 == $y2)            { die "Paired error 4: \n$line11\n$line21" }
    if (defined $opts{s}) { ; }
    elsif ($r1 == $r2)             { die "Paired error 5: \n$line11\n$line21" }
    if (check_four_lines($line11, $line12, $line13, $line14))
    { print OUTPUT1 "$line11$line12$line13$line14" }
    if (check_four_lines($line21, $line22, $line23, $line24))
    { print OUTPUT2 "$line21$line22$line23$line24" }
  }
  else                            { die "Format error: 3: \n$line11\n$line21" }
}

$line21 = <INPUT2>;
if (defined $line21)  { print STDERR "Warning: $opts{b} has more lines.\n" }

close INPUT1;
close INPUT2;
close OUTPUT1;
close OUTPUT2;
close OUTPUT3;

if (-z $opts{e}) { unlink $opts{e} }
else { warn "A fastq file $opts{e} has been created." }

exit;


sub output_single
{
  my ($l1, $l2, $l3, $l4) = @_;
  if (check_four_lines($l1, $l2, $l3, $l4)) { print OUTPUT3 "$l1$l2$l3$l4" }
}


sub check_four_lines
{
  my ($l1, $l2, $l3, $l4) = @_;

  $l1 =~ m/^\@./                                 or die "Format error 51: $l1";
  $l2 =~ m/^[TCAGN]+\n$/                         or die "Format error 52: $l2";
  $l3 =~ m/^\+\n$/                               or die "Format error 52: $l3";
  my $length = length($l2);
  unless (defined $opts{x})
  {
    unless ($len == ($length - 1))
    {
      if ($len == 0) { $len = $length - 1 }	# for the first time
      elsif (1 == (abs ($len - ($length - 1))))	# overlook 1-bp difference
      { $len = $length - 1 }
      else       { $length--; die "Length error 1: \n$len\n$length\n$l2\n$l4" }
    }
    unless ($length == length($l4))        { die "Length error 2: \n$l2\n$l4" }
  }
  my $score;
  for (my $i = 0; $i < ($length - 1); $i++)	# decrement for the new line
  {
    $score = ord(substr $l4, $i, 1) - 33;
    unless (1 < $score and $score < 42)          { die "Format error 53: $l4" }
	# min score: 2, max score: 41
  }
  return 1;	# if it did not die (if it passed this filter)
}
