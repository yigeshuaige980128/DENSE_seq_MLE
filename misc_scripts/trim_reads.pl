#!/usr/bin/env perl
use strict;
use warnings;

my $pollen = $ARGV[0] // die "Usage: perl trim_reads.pl file1.fq file2.fq\n";
my $leaf   = $ARGV[1] // die "Usage: perl trim_reads.pl file1.fq file2.fq\n";

my %pollen_lengths;
my $total = 0;

print STDERR "IMPORTING POLLEN DISTRIBUTION\n";
open my $pol, '<', $pollen or die "Cannot open $pollen: $!\n";
while (my $h = <$pol>) {
    my $read  = <$pol>;
    my $plus  = <$pol>;
    my $qual  = <$pol>;
    last if !defined $qual;

    chomp $read;
    $pollen_lengths{length($read)}++;
    $total++;
}
close $pol;

print STDERR "SUBSAMPLING LEAF READS\n";
open my $lf, '<', $leaf or die "Cannot open $leaf: $!\n";
while (my $h = <$lf>) {
    my $read  = <$lf>;
    my $plus  = <$lf>;
    my $qual  = <$lf>;
    last if !defined $qual;

    chomp $read;
    chomp $qual;

    my $length = length($read);
    for my $l (reverse 20 .. $length) {
        if (exists $pollen_lengths{$l} && $pollen_lengths{$l} > 0) {
            print $h, substr($read, 0, $l), "\n", $plus, substr($qual, 0, $l), "\n";
            $pollen_lengths{$l}--;
            delete $pollen_lengths{$l} if $pollen_lengths{$l} == 0;
            last;
        }
    }

    last if keys(%pollen_lengths) == 0;
}
close $lf;