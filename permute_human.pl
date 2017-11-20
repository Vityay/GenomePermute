#!/usr/bin/perl -w
use strict;

## 20171120 Victor Guryev victor.guryev@gmail.com

my $usage = "Usage: perl $0 <feature_file.bed> <annotaion_file.bed>\n";

my $iterations = 1000;
my $flank = 0;
my $feat_file = shift or die $usage;
my $anno_file = shift or die $usage;

my $analysis = $feat_file.'_'.$anno_file.'_'.$iterations.'_permutations_'.$flank.'bp_flank';

# Reading features file
warn "Reading feature file $feat_file ...\n";
my %sces = ();
my %sces_anno = ();
my ( $i, $j, $max, $line ) = ( 0, 0, 0, 0 );
open F, $feat_file;
while ( <F> ) {
    $line++;
    next if m/^track/;
    next if m/^chromoso?me/i;
    next if m/^\#/;
    chomp;
    my ( $chr, $start, $end, @rest ) = split /\s+/;
    $chr =~ s/^chr//;
    next unless $chr =~ m/^[\dXY]{1,2}$/;
    $start -= $flank;
    $start = 1 if $start < 1;
    $end += $flank;

    if ( !exists($sces{$chr}{$start}) or $sces{$chr}{$start} < $end ) {
        $sces{$chr}{$start} = $end;
        $sces_anno{$chr}{$start} = join( "\t", @rest );
    }
    $i+= $end-$start+1;
    $j++;
    $max = $end-$start+1 if $end-$start +1 > $max;
}
close F;
die "Seen NO features in Features file, terminating this analysis\n" if $j == 0;
warn "Feature set: $feat_file $j features with total size $i bp and biggest feature: $max bp\n";

# Reading annotation file
warn "Reading annotation file $anno_file ...\n";
my %anno = ();
my %anno_anno = ();
( $i, $j, $max, $line ) = ( 0, 0, 0, 0 );
open F, $anno_file;
while ( <F> ) {
    $line++;
    next if m/^track/;
    next if m/^chromoso?me/i;
    next if m/^\#/;
    chomp;
    my ( $chr, $start, $end, @rest ) = split /\s+/;
    $chr=~ s/^chr//;
    next unless $chr =~ m/^[\dXY]{1,2}$/;
    $i+= $end-$start+1;
    $j++;
    $max = $end-$start+1 if $end-$start +1 > $max;
    if ( !exists($anno{$chr}{$start}) or $anno{$chr}{$start}<$end ) {
        $anno{$chr}{$start}=$end;
        $anno_anno{$chr}{$start}=join( "\t", @rest );
    }
}
close F;
die "Seen NO features in Annotation file, terminating this analysis\n" if $j == 0;
warn "Annotation set: $anno_file $j features with total size $i bp, biggest feature $max bp\n";

# Generate random shifts for permutation analysis
warn "Generating random shift values ...\n";
my @shifts = ( 0 );
foreach ( 1 .. $iterations ) {
    my $shift = int(rand( 50_000_000 ));
    redo if $shift < 1_000_000;
    push @shifts, $shift;
}

warn "Loading chromosome lengths\n";
my %chrlen = ();
$i = 0;
open( F, 'GRCh37_chrlen.txt' ) or die 'Cannot find file GRCh37_chrlen.txt';
while ( <F> ) {
    chomp;
    my ( $chr, $len ) = split /\t/;
    $chrlen{$chr} = $len;
    $i+=$len;
}
close F;
warn scalar keys %chrlen, " chromosomes, total size: $i\n";

my %gaps = ();
warn "Reaging gap locations...\n";
open( F, 'GRCh37_gaps.bed' ) or die 'Cannot find file GRCh37_gaps.bed';
while ( <F> ) {
    chomp;
    my ( $chr, $start, $end ) = split /\t/;
    $gaps{$chr}{$start} = $end;
}
close F;

warn "Removing gaps from datasets...\n";
foreach my $chr ( keys %gaps ) {
    foreach my $s ( sort {$b<=>$a} keys %{$gaps{$chr}} ) {
        my $gap_size = $gaps{$chr}{$s} - $s + 1;
        my $e = $gaps{$chr}{$s};

        # Chr
        $chrlen{$chr} -= $gap_size;

        # Annotation Features
        foreach my $start ( sort {$a<=>$b} keys %{$anno{$chr}} ) {
            next if $anno{$chr}{$start} < $s;
            my $end = $anno{$chr}{$start};
            my $extended = $anno_anno{$chr}{$start};
            my $new_start = $start > $e ? ( $start - $gap_size ) : $start ;
            my $new_end = $end > $e ? ( $end - $gap_size ) : $end;
            delete $anno{$chr}{$start};
            $anno{$chr}{$new_start}=$new_end;
            delete $anno_anno{$chr}{$start};
            $anno_anno{$chr}{$new_start}=$extended;
        }

        # SCEs
        foreach my $start ( sort {$a<=>$b} keys %{$sces{$chr}} ) {
            next if $sces{$chr}{$start} < $s;
            my $end = $sces{$chr}{$start};
            my $extended = $sces_anno{$chr}{$start};
            my $new_start = $start > $e ? ( $start - $gap_size ) : $start ;
            my $new_end = $end > $e ? ( $end - $gap_size ) : $end ;
            delete $sces{$chr}{$start};
            $sces{$chr}{$new_start}=$new_end;
            delete $sces_anno{$chr}{$start};
            $sces_anno{$chr}{$new_start}=$extended;
        }
    }
}


warn "Sorting annotations\n";
my %anno_srt = ();
foreach my $chr ( keys %anno ) {
    @{$anno_srt{$chr}} = sort {$a<=>$b} keys %{$anno{$chr}};
}

open F, '>', $analysis.'.txt';
print F 'Shift', "\t", 'Overlaps', "\n";
$i = 0;
my @rand_overlaps = ();
my $real_overlap = q{};
my ( $direction, $p_value ) = ( '?', '?' );
warn "Processing\n";
foreach my $shift ( @shifts ) {
    warn 'Iteration: ',$i++,'/',$iterations,', shift:',$shift,"\n";
    my $hits = 0;
    foreach my $chr ( keys %sces ) {
        @{$anno_srt{$chr}} = () unless exists($anno_srt{$chr});
        my @anno_starts = @{$anno_srt{$chr}};
        my %sces_starts = ();
        foreach my $bstart ( keys %{$sces{$chr}} ) {
            my ( $s, $e ) = ( $bstart+$shift, $sces{$chr}{$bstart}+$shift);
            die unless $chrlen{$chr}; 
            while ( $s > $chrlen{$chr} ) {
                $s -= $chrlen{$chr};
                $e -= $chrlen{$chr};
            }
            $sces_starts{$s}=$e;
        }
        foreach my $bstart ( sort {$a<=>$b} keys %sces_starts ) {
            last unless @anno_starts;
            while ( @anno_starts and $anno_starts[0] < $bstart - 10_000 ) {
                shift @anno_starts;
            }
            my $gin = 0;
            foreach my $rstart ( @anno_starts ) {
                last if $rstart > $sces_starts{$bstart};
                next if $bstart > $anno{$chr}{$rstart} or $rstart > $sces_starts{$bstart};
                $gin=1;
                last if $shift;
            }
            $hits++ if $gin;
        }
    }
    if ( $shift ) {
        push @rand_overlaps, $hits;
    }
    else {
        $real_overlap = $hits;
    }
    my ( $p_over, $p_under ) = (1, 1);
    if ( @rand_overlaps ) {
        $p_over  = scalar( grep {$_ < $real_overlap } @rand_overlaps );
        $p_under = scalar( grep {$_ > $real_overlap } @rand_overlaps );
        $p_over /= @rand_overlaps;
        $p_under /= @rand_overlaps;
    }
    ( $direction, $p_value ) = $p_over < $p_under ? ( 'Depleted', $p_over ) : ( 'Enriched', $p_under );
    $p_value = $p_value == 0 ? 'p<'.(1/@rand_overlaps) : 'p='.$p_value;
    print F join( "\t", $shift, $hits ), "\n";
    warn ' ',$hits, " overlaps, Current estimate: $direction, $p_value\n";
}
close F;
print "=== FINAL CONCLUSION: $direction, $p_value ===\n";
