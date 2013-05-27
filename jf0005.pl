#!/usr/bin/env perl

=head1 NAME

  jf0005.pl - remove multiple-hit reads from a SAM file

=head1 SYNOPSIS

  $ perldoc jf0005.pl
  $ jf0005.pl -h
  $ jf0005.pl all_hits.sam > unique_hits_only.sam

=head1 DESCRIPTION

  This script reads a SAM file which was created by bwa sampe,
  and removes multiple-hit reads; only unique hits remain.
  The reduced output is also a SAM file that is printed onto
  the standard output.

=head1 OPTIONS

  -h  Print usage

=head1 AUTHOR

  Coded by Kohji OKAMURA, Ph.D.

=head1 HISTORY

  Aug 31, 2012  The first version from scratch
  May 27, 2013  POD is modified

=cut


use warnings;
use strict;
use English;
use Getopt::Std;

my %opts;
my ($line1, $line2);  # to read a pair of reads (R1 and R2)
my ($id1, $id2);
my ($mapq1, $mapq2);
my ($data1, $data2);

getopts 'h', \%opts;

if (defined $opts{h})                                     { exec "perldoc $0" }

while ($line1 = <>)
{
  if ($line1 =~ m/^\@/) { print $line1; next }	# header lines
  else
  {
    $line2 = <>;	# paired one
    $line1 =~ m/^(\S+)\t\S+\t\S+\t\S+\t(\d+)\t(.+)\n$/
                                               or die "Format error 1: $line1";
    ($id1, $mapq1, $data1) = ($1, $2, $3);
    $line2 =~ m/^(\S+)\t\S+\t\S+\t\S+\t(\d+)\t(.+)\n$/
                                               or die "Format error 2: $line2";
    ($id2, $mapq2, $data2) = ($1, $2, $3);
    unless ($id1 eq $id2)                   { die "Format error 3: $id1 $id2" }
    if ($mapq1 == 0 or $mapq2 == 0) { next }
    $data1 =~ m/\sXT:A:U\s/ or next;
    $data2 =~ m/\sXT:A:U\s/ or next;

    print $line1, $line2;
  }
}
