#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use POSIX qw(floor);

my %options = (
    data_file => 'sites_ind1.cm.txt',
    rate_file => '',
    min_size  => 1000,
    error     => 1e-2,
    tolerance => 1e-4,
    window    => 0,
    bootstrap => 0,
);

GetOptions(
    'd=s' => \$options{data_file},
    'r=s' => \$options{rate_file},
    'm=i' => \$options{min_size},
    'w=i' => \$options{window},
    'b=i' => \$options{bootstrap},
    'e=f' => \$options{error},
    't=f' => \$options{tolerance},
) or die "Usage: perl map_sd.pl -d <data_file> [-r rate_file] [-m int] [-w int] [-b int] [-e float] [-t float]\n";

my ($germline, $somatic) = import_data($options{data_file});
my $rates = {};
$rates = import_rates($options{rate_file}) if $options{rate_file} ne '';

for my $chrom (sort keys %{$germline}) {
    my $g_sites = $germline->{$chrom};
    my $s_sites = $somatic->{$chrom};

    next if @{$g_sites} == 0;

    if ($options{window} != 0 && $options{bootstrap} == 0) {
        run_window($chrom, $g_sites, $s_sites, \%options);
    }
    elsif ($options{bootstrap} > 0) {
        run_bootstrap($chrom, $g_sites, $s_sites, \%options);
    }
    else {
        run_single_fit($chrom, $g_sites, $s_sites, \%options);
    }
}

exit 0;

sub run_single_fit {
    my ($chrom, $g_sites, $s_sites, $opts) = @_;

    my $optimum_site = golden_search_position($opts, $s_sites, $g_sites);
    my $somatic_k = golden_search_site($opts, $s_sites, $optimum_site);
    my $germline_k = golden_search_site($opts, $g_sites, $optimum_site);

    my $lnl_s = compute_likelihood($opts->{error}, $g_sites, $optimum_site, $somatic_k);
    my $lnl_g = compute_likelihood($opts->{error}, $g_sites, $optimum_site, $germline_k);

    print join("\t", $chrom, fmt($optimum_site), fmt($somatic_k), fmt($germline_k), fmt($lnl_s), fmt($lnl_g)), "\n";
}

sub run_bootstrap {
    my ($chrom, $g_sites, $s_sites, $opts) = @_;

    my $n = scalar @{$g_sites};
    return if $n == 0;

    for (my $boot_idx = 0; $boot_idx < $opts->{bootstrap}; $boot_idx++) {
        my @bootstrap_germline;
        my @bootstrap_somatic;

        for (my $i = 0; $i < $n; $i++) {
            my $index = int(rand($n));
            push @bootstrap_germline, { %{$g_sites->[$index]} };
            push @bootstrap_somatic, { %{$s_sites->[$index]} };
        }

        @bootstrap_germline = sort { $a->{position} <=> $b->{position} } @bootstrap_germline;
        @bootstrap_somatic = sort { $a->{position} <=> $b->{position} } @bootstrap_somatic;

        my $boot_site = golden_search_position($opts, \@bootstrap_somatic, \@bootstrap_germline);

        # Keep behavior consistent with upstream implementation: k is optimized on original arrays.
        my $somatic_k = golden_search_site($opts, $s_sites, $boot_site);
        my $germline_k = golden_search_site($opts, $g_sites, $boot_site);

        my $lnl_s = compute_likelihood($opts->{error}, \@bootstrap_germline, $boot_site, $somatic_k);
        my $lnl_g = compute_likelihood($opts->{error}, \@bootstrap_germline, $boot_site, $germline_k);

        print join("\t", $chrom, $boot_idx, fmt($boot_site), fmt($somatic_k), fmt($germline_k), fmt($lnl_s), fmt($lnl_g)), "\n";
    }
}

sub run_window {
    my ($chrom, $g_sites, $s_sites, $opts) = @_;

    my $half = int($opts->{window} / 2);
    my $end = $g_sites->[-1]{position} - $half;

    for (my $site = $half; $site <= $end; $site += $opts->{window}) {
        my $somatic_k = golden_search_site($opts, $s_sites, $site);
        my $germline_k = golden_search_site($opts, $g_sites, $site);

        my $lnl_s = compute_likelihood($opts->{error}, $g_sites, $site, $somatic_k);
        my $lnl_g = compute_likelihood($opts->{error}, $g_sites, $site, $germline_k);

        print join("\t", $chrom, $site, fmt($somatic_k), fmt($germline_k), fmt($lnl_s), fmt($lnl_g), fmt($lnl_s - $lnl_g)), "\n";
    }
}

sub import_data {
    my ($file) = @_;

    open my $fh, '<', $file or die "Cannot open input file '$file': $!\n";

    my (%germline, %somatic);

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/^\x{FEFF}//;
        next if $line =~ /^\s*$/;

        my ($chrom, $position, $cm, $sA, $sa, $gA, $ga) = split /\t/, $line;
        next if !defined $ga;

        push @{$somatic{$chrom}}, {
            position => int($position),
            cm       => 0.0 + $cm,
            A        => int($sA),
            a        => int($sa),
        };

        push @{$germline{$chrom}}, {
            position => int($position),
            cm       => 0.0 + $cm,
            A        => int($gA),
            a        => int($ga),
        };
    }

    close $fh;

    for my $chrom (keys %somatic) {
        my @s = sort { $a->{position} <=> $b->{position} } @{$somatic{$chrom}};
        my @g = sort { $a->{position} <=> $b->{position} } @{$germline{$chrom}};
        $somatic{$chrom} = \@s;
        $germline{$chrom} = \@g;
    }

    return (\%germline, \%somatic);
}

sub import_rates {
    my ($file) = @_;
    my %rates;

    open my $fh, '<', $file or die "Cannot open rate file '$file': $!\n";
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^\s*$/;
        my ($chrom, $rate) = split /\s+/, $line;
        next if !defined $rate;
        $rates{$chrom} = 0.0 + $rate;
    }
    close $fh;

    return \%rates;
}

sub compute_likelihood {
    my ($error, $sites, $site_position, $k) = @_;

    my $site_i = int($site_position);
    my $cm_at_site = find_cm_at_site($sites, $site_i);
    my $lnl = 0.0;

    for my $s (@{$sites}) {
        my $rec = (1 - exp(-2 * (abs($s->{cm} - $cm_at_site) / 100))) / 2;

        my $lA = 0.0;
        $lA += (1 - $error) * (1 - $rec) * $k;
        $lA += (1 - $error) * $rec * (1 - $k);
        $lA += $error * (1 - $rec) * (1 - $k);
        $lA += $error * $rec * $k;

        my $la = 0.0;
        $la += (1 - $error) * (1 - $rec) * (1 - $k);
        $la += (1 - $error) * $rec * $k;
        $la += $error * (1 - $rec) * $k;
        $la += $error * $rec * (1 - $k);

        $lA = clamp_prob($lA);
        $la = clamp_prob($la);

        $lnl += log($lA) * $s->{A};
        $lnl += log($la) * $s->{a};
    }

    return $lnl;
}

sub find_cm_at_site {
    my ($sites, $site_position) = @_;

    my $cm = $sites->[0]{cm};
    my $min = 10**15;

    for my $s (@{$sites}) {
        my $d = abs($s->{position} - $site_position);
        if ($d < $min) {
            $min = $d;
            $cm = $s->{cm};
        }
    }

    return $cm;
}

sub golden_search_site {
    my ($opts, $sites, $site_position) = @_;

    my $phi = (sqrt(5) - 1) / 2;

    my $low_bracket = 0.01;
    my $high_bracket = 0.99;
    my $param_low = $high_bracket + $phi * ($low_bracket - $high_bracket);
    my $param_high = $low_bracket + $phi * ($high_bracket - $low_bracket);

    my ($lnl_low, $lnl_high);
    my $lnl_diff = 1e9;
    my $iter = 0;
    my $max_iter = 2000;

    while ($opts->{tolerance} < $lnl_diff && $iter < $max_iter) {
        if (!defined $lnl_low) {
            $lnl_low = compute_likelihood($opts->{error}, $sites, $site_position, $param_low);
        }
        if (!defined $lnl_high) {
            $lnl_high = compute_likelihood($opts->{error}, $sites, $site_position, $param_high);
        }

        $lnl_diff = abs($lnl_low - $lnl_high);

        if ($lnl_high >= $lnl_low) {
            $low_bracket = $param_low;
            $param_low = $param_high;
            $param_high = $low_bracket + ($high_bracket - $low_bracket) * $phi;
            $lnl_low = $lnl_high;
            undef $lnl_high;
        }
        else {
            $high_bracket = $param_high;
            $param_high = $param_low;
            $param_low = $high_bracket + ($low_bracket - $high_bracket) * $phi;
            $lnl_high = $lnl_low;
            undef $lnl_low;
        }

        $iter++;
    }

    return ($param_low + $param_high) / 2;
}

sub golden_search_position {
    my ($opts, $somatic, $germline) = @_;

    my $phi = (sqrt(5) - 1) / 2;

    my $low_bracket = $somatic->[0]{position};
    my $high_bracket = $somatic->[-1]{position};
    my $param_low = $high_bracket + $phi * ($low_bracket - $high_bracket);
    my $param_high = $low_bracket + $phi * ($high_bracket - $low_bracket);

    my ($lnl_low, $lnl_high);
    my $lnl_diff = 1e9;
    my $iter = 0;
    my $max_iter = 2000;

    while ($opts->{tolerance} < $lnl_diff && $iter < $max_iter) {
        if (!defined $lnl_low) {
            my $test_site = int($param_low);
            my $somatic_k = golden_search_site($opts, $somatic, $test_site);
            my $germline_k = golden_search_site($opts, $germline, $test_site);
            $lnl_low = compute_likelihood($opts->{error}, $germline, $test_site, $germline_k)
                     - compute_likelihood($opts->{error}, $germline, $test_site, $somatic_k);
        }

        if (!defined $lnl_high) {
            my $test_site = int($param_high);
            my $somatic_k = golden_search_site($opts, $somatic, $test_site);
            my $germline_k = golden_search_site($opts, $germline, $test_site);
            $lnl_high = compute_likelihood($opts->{error}, $germline, $test_site, $germline_k)
                      - compute_likelihood($opts->{error}, $germline, $test_site, $somatic_k);
        }

        $lnl_diff = abs($lnl_low - $lnl_high);

        if ($lnl_high >= $lnl_low) {
            $low_bracket = $param_low;
            $param_low = $param_high;
            $param_high = $low_bracket + ($high_bracket - $low_bracket) * $phi;
            $lnl_low = $lnl_high;
            undef $lnl_high;
        }
        else {
            $high_bracket = $param_high;
            $param_high = $param_low;
            $param_low = $high_bracket + ($low_bracket - $high_bracket) * $phi;
            $lnl_high = $lnl_low;
            undef $lnl_low;
        }

        $iter++;
    }

    return ($param_low + $param_high) / 2;
}

sub clamp_prob {
    my ($x) = @_;
    return 1e-300 if $x <= 0;
    return 1 - 1e-15 if $x >= 1;
    return $x;
}

sub fmt {
    my ($x) = @_;
    return sprintf('%.10g', $x);
}
