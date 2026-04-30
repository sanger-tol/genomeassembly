#!/usr/bin/env perl
#
#    Copyright (C) 2018-2019 Genome Research Ltd.
#
#    Author: Shane McCarthy <sm15@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

use strict;
use warnings;
use Carp;
use Getopt::Long;

my $opts = parse_params();
stats($opts,'SCAFFOLD');

# 'seqtk comp' output format:
# 0        1       2       3       4       5       6       7       8   9       10  11  12
# name     length  A       C       G       T       IUPAC2  IUPAC3  N   CpG     ts  tv  CpG-ts
# Contig0  7097294 2186969 1367324 1358672 2183704 0       0       625 253444  0   0   0

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    die
        "About: Generate stats for fasta assembly file\n",
        "Usage: asmstats [OPTIONS] <fasta|fastq>\n",
        "Options:\n",
        "   -g, --genome-size INT[kMG]  estimated genome size for NGX stats\n",
        "       --seqtk <path>          path to seqtk executable [seqtk]\n",
        "   -h, --help                  this help message.\n";
}


sub parse_params
{
    my $opts = { fa => '-', seqtk => 'seqtk', iter => 0 };
    while (defined(my $arg=shift(@ARGV)))
    {
        if (                 $arg eq '--seqtk' ) { $$opts{seqtk} = shift(@ARGV); next; }
        if ( $arg eq '-g' || $arg eq '--genome-size' ) { $$opts{genome_size} = shift(@ARGV); next; }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        if (scalar @ARGV == 0) { $$opts{fa} = $arg; next; }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    error() if ($$opts{fa} eq '-' && $$opts{contigs} && $$opts{scaffolds});
    unless ($$opts{fa} eq '-')
    {
        error(qq[$$opts{fa} does not exist]) unless (-s $$opts{fa});
    }
    if (exists $$opts{genome_size} && $$opts{genome_size} =~ m/(\d+)([kmgKMG])$/)
    {
        if (uc($2) eq 'K') { $$opts{genome_size} = $1 * 1e3; }
        if (uc($2) eq 'M') { $$opts{genome_size} = $1 * 1e6; }
        if (uc($2) eq 'G') { $$opts{genome_size} = $1 * 1e9; }
    }
    return $opts;
}

sub stats
{
    my ($opts,$type) = @_;

    my @sums;
    my %n50;

    # transform input depending on type of sequence input and look at composition with 'seqtk comp'
    my $input;
    if ($type eq 'GAP')
    {
        # convert to a fasta file of just the N sequence
        $input = qq[$$opts{seqtk} cutN -n1 -g $$opts{fa} | ] . q[awk '{print ">"NR; i = 0; while (i++ < $3-$2) { printf "N"; } print "\n"}' | ] . qq[$$opts{seqtk} comp - | ];
    }
    elsif ($type eq 'CONTIG')
    {
        # split sequence at any N to get contigs
        $input = qq[$$opts{seqtk} cutN -n1 $$opts{fa} | $$opts{seqtk} comp - | ];
    }
    else
    {
        $input = qq[$$opts{seqtk} comp $$opts{fa} | ];
    }

    # read in 'seqtk comp' output
    open(my $fh, "$input") or die("Could not open '$input' for reading");
    while (my $line = <$fh>)
    {
        chomp $line;
        my (@a) = split "\t", $line;
        for (my $i=1; $i<scalar @a; $i++)
        {
            $sums[$i-1]+=$a[$i];
        }
        $n50{$a[0]}{len} = $a[1]; # scaffold length
        $n50{$a[0]}{seq} = $a[1]-$a[8] if ($a[8]); # amount of non-N sequence in scaffold
    }
    close($fh) or die("Could not close '$input'");

    my $n = scalar keys %n50;
    my $sum = $sums[0];
    my $ns = $sums[7];
    my $sum50 = 0;
    my $seqsum50 = 0;
    my $vp = 0.5;
    my $wp = 0.5;

    # if there are no Ns, we'll only produce contig stats
    $type = 'CONTIG' unless ($ns);

    if ($n > 0)
    {
        my $result = $sum / $n;
        my @keys = sort { $n50{$b}{len} <=> $n50{$a}{len} } keys %n50;
        unless ($$opts{iter})
        {
            # use Data::Dumper;
            # print Dumper(\@sums);
            # on the first pass through, print fasta input path, composition, Ns (if any) and IUPAC related counts (if any)
            print "$$opts{fa}\n";
            printf(qq[A = $sums[1] (%.1f%%), C = $sums[2] (%.1f%%), G = $sums[3] (%.1f%%), T = $sums[4] (%.1f%%)], 100*$sums[1]/$sum, , 100*$sums[2]/$sum, 100*$sums[3]/$sum, 100*$sums[4]/$sum);
            printf(qq[, N = $ns (%.1f%%)], 100*$ns/$sum) if ($ns);
            printf(qq[, CpG = $sums[8] (%.1f%%)], 100*$sums[8]/$sum) if ($sums[8]);
            print qq[, IUPAC3 = $sums[6], IUPAC2 = $sums[5], ts = $sums[9], tv = $sums[10], CpG-ts = $sums[11]] if ($sums[5] || $sums[6]);
            print "\n";
        }
        print "$type\t" unless ($$opts{iter} == 0 && $type eq 'CONTIG');
        print "sum = $sum, n = $n, ave = $result, largest = $n50{$keys[0]}{len}, smallest = $n50{$keys[-1]}{len}\n";
        return if ($type eq 'GAP'); # no NX stats for GAP

        my %res;
        my $m = 0;
        foreach my $key (@keys)
        {
            $sum50 += $n50{$key}{len};
            $seqsum50 += $n50{$key}{seq} if ($n50{$key}{seq});
            $m++;
            if ( $sum50/$vp > $sum )
            {
                # print "$type\t" unless ($$opts{iter} == 0 && $type eq 'CONTIG');
                # printf "N%.0f = $n50{$key}{len}, L%.0f = $m\n", $vp*100, $vp*100;
                $res{100*$vp}{N} = $n50{$key}{len};
                $res{100*$vp}{L} = $m;
                # printf "$type\tN%.0f = $n50{$key}{len}, L%.0f = $m, T%.0f = $seqsum50 (%.1f%%\)\n", $vp*100, $vp*100, $vp*100, 100*$seqsum50/$sum;
                $vp += .1;
            }
            if ( exists $$opts{genome_size} && $sum50/$wp > $$opts{genome_size} )
            {
                # print "$type\t" unless ($$opts{iter} == 0 && $type eq 'CONTIG');
                # printf "NG%.0f = $n50{$key}{len}, LG%.0f = $m\n", $wp*100, $wp*100, $wp*100;
                $res{100*$wp}{NG} = $n50{$key}{len};
                $res{100*$wp}{LG} = $m;
                # printf ", NG%.0f = $n50{$key}{len}, LG%.0f = $m, TG%.0f = $seqsum50 (%.1f%%\)\n", $wp*100, $wp*100, $wp*100, 100*$seqsum50/$$opts{genome_size};
                $wp += .1;
            }
        }
        if ($n > 1)
        {
            foreach my $val (50,60,70,80,90,100)
            {
                print "$type\t" unless ($$opts{iter} == 0 && $type eq 'CONTIG');
                printf "N%.0f = %d, L%.0f = %d", $val, $res{$val}{N}, $val, $res{$val}{L};
                if ($res{$val}{NG})
                {
                    printf ", NG%.0f = %d, LG%.0f = %d", $val, $res{$val}{NG}, $val, $res{$val}{LG};
                }
                print "\n";
            }
        }
        $$opts{iter}++;
        if ($ns && $type eq 'SCAFFOLD')
        {
            stats($opts,'CONTIG');
            stats($opts,'GAP');
        }
    }
}
