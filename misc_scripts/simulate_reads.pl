#!/usr/bin/env perl
use strict;
use warnings;

my $error = 1e-2;

my $out_dir = 'BOOT';
mkdir $out_dir unless -d $out_dir;

my (%somatic_depth, %germline_depth, %cm);

while (my $line = <STDIN>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    my @split = split /\t/, $line;
    next if @split < 7;

    $somatic_depth{$split[1]} = $split[3] + $split[4];
    $germline_depth{$split[1]} = $split[5] + $split[6];
    $cm{$split[1]} = 0.0 + $split[2];
}

my @positions = keys %cm;
for my $effect (0, 0.01) {
    for my $rep (1) {
        my $site = $positions[int(rand(@positions))];
        open my $out, '>', "$out_dir/${rep}_${effect}_${site}.txt" or die "Cannot write output: $!\n";

        for my $position (sort { $a <=> $b } keys %somatic_depth) {
            my $drawA = random_binomial($somatic_depth{$position}, 0.5);
            print {$out} join("\t", 'scaffold_X', $position, $cm{$position}, $drawA, $somatic_depth{$position} - $drawA);
            print {$out} "\t";

            my $rec = (1 - exp(-2 * abs($cm{$site} - $cm{$position}) / 100)) / 2;

            $drawA = random_binomial($germline_depth{$position}, 0.5 + $effect);

            my $swapA = random_binomial($drawA, $rec);
            my $swapa = random_binomial($germline_depth{$position} - $drawA, $rec);
            $drawA = $drawA - $swapA + $swapa;

            $swapA = random_binomial($drawA, $error);
            $swapa = random_binomial($germline_depth{$position} - $drawA, $error);
            $drawA = $drawA - $swapA + $swapa;

            print {$out} join("\t", $drawA, $germline_depth{$position} - $drawA), "\n";
        }

        close $out;
    }
}

sub random_binomial {
    my ($n, $p) = @_;
    return 0 if $n <= 0;
    $p = 0 if $p < 0;
    $p = 1 if $p > 1;

    my $count = 0;
    for (1 .. $n) {
        $count++ if rand() < $p;
    }
    return $count;
}