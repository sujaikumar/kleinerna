#!/usr/bin/env perl

# Script for taking a file with two columns - one with small RNA, one with target RNA, and finding 
# * seeds: 6mer, 6mer_offset, 7mer_m8, 7mer_A1, 8mer
# * whether it is low complexity (uses dustmasker from NCBI blast+ suite)
# * what the targetfinder pairing/score is for each pair. Rather than call targetfinder.pl separately, the code for targetfinder is reused in this script

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);
use List::Util qw/ sum /;
use File::Temp qw/ tempfile tempdir /;

my ($seeds, $dust, $targetfinder, $separator) = (1,1,1,",");

GetOptions (
  "seeds!"        => \$seeds,
  "dust!"         => \$dust,
  "targetfinder!" => \$targetfinder,
  "separator=s"   => \$separator,
);

## CODE FROM TARGETFINDER
## https://github.com/carringtonlab/TargetFinder/tree/848b2dddf065f1f84664c8b6fe11be5b2ed6b8e1

# miRNA-target alignment scoring matrix
my %bp;
$bp{"AU"} = 0;
$bp{"UA"} = 0;
$bp{"GC"} = 0;
$bp{"CG"} = 0;
$bp{"GU"} = 0.5;
$bp{"UG"} = 0.5;
$bp{"AC"} = 1;
$bp{"CA"} = 1;
$bp{"AG"} = 1;
$bp{"GA"} = 1;
$bp{"UC"} = 1;
$bp{"CU"} = 1;
$bp{"A-"} = 1;
$bp{"U-"} = 1;
$bp{"G-"} = 1;
$bp{"C-"} = 1;
$bp{"-A"} = 1;
$bp{"-U"} = 1;
$bp{"-G"} = 1;
$bp{"-C"} = 1;
$bp{"AA"} = 1;
$bp{"UU"} = 1;
$bp{"CC"} = 1;
$bp{"GG"} = 1;

while (<>) {
  next unless /^([atgcnu]+)[ \t,;]([atgcnu]+)$/i;
  my ($smallrna, $targetrna) = (uc $1, uc $2);
  
  my ($s_fh, $s_fn) = tempfile("stmpXXXXXXX", DIR => "."); print $s_fh ">small\n$smallrna";   close $s_fh;
  my ($t_fh, $t_fn) = tempfile("ttmpXXXXXXX", DIR => "."); print $t_fh ">target\n$targetrna"; close $t_fh;
  
  if ($targetfinder) {
    my @tf_results = &targetfinder ( $s_fn , $t_fn );
    if (@tf_results) { print $separator . join($separator, $tf_results[0]{target_seq}, $tf_results[0]{homology_string}, $tf_results[0]{miR_seq}, $tf_results[0]{score} ) }
    else { print $separator x 4 };
  }
  
  if ($dust) {
    my $dusted = `dustmasker -in $t_fn -outfmt fasta | seqtk seq | tail -n 1`;
    my $is_dust = 1;
    while ( $dusted=~/([ATGCN]+)/g) {
      if  ( length($1) >= 25 ) {
        $is_dust = 0;
        last;
      }
    }
    print $separator . $is_dust;
  }

  if ($seeds) {
    print ",". &is_6mer_match     ( $smallrna, $targetrna );
    print ",". &is_7mer_m8_match  ( $smallrna, $targetrna );
    print ",". &is_7mer_A1_match  ( $smallrna, $targetrna );
    print ",". &is_8mer_match     ( $smallrna, $targetrna );
    print ",". &is_6mer_off_match ( $smallrna, $targetrna );
  }

  print "\n";
  unlink $s_fn, $t_fn;

}


### END MAIN ###

#############################################################################


sub targetfinder {
  my ($s_fn, $t_fn) = @_;

## CODE FROM TARGETFINDER
  my @fasta        = fasta($s_fn, $t_fn);
  my @fasta_parsed = fasta_parser(@fasta);
  my @targets      = bp_score(@fasta_parsed);
  return @targets;
}

## CODE FROM TARGETFINDER
########################################
# Function: fasta
#      Execute the Smith-Waterman
#      alignment program
########################################
sub fasta {
        my $input = shift;
        my $db = shift;
        my @output;
        open FASTA, "ssearch36 -n -H -Q -f -16 -r +15/-10 -g -10 -w 100 -W 25 -E 100000 -i -U $input $db 1 2> /dev/null |";
        while (<FASTA>) {
                push (@output, $_);
        }
        close FASTA;
        return @output;
}
########################################
# Function: fasta_parser
#     Parse results from FASTA/SW output
########################################
sub fasta_parser {
        my @input = @_;
        my $hit_accession;
        my $id;
        my $spacer;
        my $miRNA;
        my $count = 0;
        my $length = 0;
        my $target;
        my $step = 0;
        my @output;
        foreach my $line (@input) {
                chomp $line;
                if ($step == 0) {
                        if ($line =~ /^\>{2}(.+)\(\d+ nt\)/) {
                                $hit_accession = $1;
                                if ($hit_accession =~ /(^.{6})/) {
                                        $id = $1;
                                }
                                while ($hit_accession =~ /\s$/) {
                                        $hit_accession =~ s/\s$//g;
                                }
                                $step = 1;

                        }
                } elsif ($step == 1) {
                        if ($line =~ /(small\-\s+)\b(\D+)\b/) {   #if ($line =~ /($name\-\s+)\b(\D+)\b/) {
                                $spacer = $1;
                                $miRNA = $2;
                                $count = ($spacer =~ tr/[A-Z][a-z][0-9][\-][\ ][_]//);
                                $miRNA =~ s/\s//g;
                                $length = length $miRNA;
                                $miRNA =~ tr/AUGC/UACG/;

                        } elsif ($line =~ /$id\s+/) {
                                if ($line =~ /.{$count}(\D{$length})/) {
                                        $target = $1;
                                        $target =~ s/\s//g;

                                }
                                push @output, "$hit_accession\t$miRNA\t$target";
                                $step = 0;
                        }
                }

        }
        return @output;
}
########################################
# Function: bp_score
#     Score alignments for base-pairing
########################################
sub bp_score {

        my @input = @_;
        my $counter = 0;
        my @output;
        foreach my $line (@input) {
                chomp $line;
                my $mismatch = 0;
                my $gu = 0;
                my $total_mispair = 0;
                my ($hit_accession, $miR_seq, $target_seq) = split /\t/, $line;
                my $miR_length = length $miR_seq;
                my $target_length = length $target_seq;
                next if ($target_seq !~ /^[AGCU-]+$/);
                next if ($miR_length != $target_length);
                next if ($miR_seq =~ /\-.{0,}\-/);
                next if ($target_seq =~ /\-.{0,}\-/);
                next if ($miR_seq =~ /\-/ && $target_seq =~ /\-/);
                my @miRNA = split //, $miR_seq;
                my @target = split //, $target_seq;
                my $cycle = 0;
                my $score = 0;
                my $old_score = 0;
                my $homology_string;
                for (1..$miR_length) {
                        $cycle++;
                        my $miR_base = pop @miRNA;
                        my $target_base = pop @target;
                        if ($cycle == 1) {
                                my $position = $bp{"$miR_base$target_base"};
                                if ($position == 1) {
                                        $mismatch++;
                                        $homology_string .= ' ';
                                } elsif ($position == 0.5) {
                                        $gu++;
                                        $homology_string .= '.';
                                } else {
                                        $homology_string .= ':';
                                }
                                $score = $position;
                        } elsif ($cycle > 13) {
                                my $position = $bp{"$miR_base$target_base"};
                                if ($position == 1) {
                                        $mismatch++;
                                        $homology_string .= ' ';
                                } elsif ($position == 0.5) {
                                        $gu++;
                                        $homology_string .= '.';
                                } else {
                                        $homology_string .= ':';
                                }
                                my $new_score = $position;
                                $old_score = $score;
                                $score = ($old_score + $new_score);
                        } else {
                                my $position = ($bp{"$miR_base$target_base"}*2);
                                if ($position == 2) {
                                        $mismatch++;
                                        $homology_string .= ' ';
                                } elsif ($position == 1) {
                                        $gu++;
                                        $homology_string .= '.';
                                } else {
                                        $homology_string .= ':';
                                }
                                my $new_score = $position;
                                $old_score = $score;
                                $score = ($old_score + $new_score);
                        }


                }
                $total_mispair = ($mismatch+$gu);
                next if ($total_mispair > 100 || $mismatch > 100 || $gu > 100); # next if ($total_mispair > 7 || $mismatch > 4 || $gu > 4);
                if ($score <= 100) { #if ($score <= $cutoff) {
                        my %hash;
                        $counter++;
                        $hash{'hit_accession'} = $hit_accession;
                        $hash{'miR_seq'} = $miR_seq;
                        $hash{'query_name'} = "query";
                        $hash{'target_seq'} = $target_seq;
                        $hash{'score'} = $score;
                        $hash{'homology_string'} = reverse $homology_string;
                        $hash{'target_start'} = 0;
                        $hash{'target_end'} = 0;
                        $hash{'target_strand'} = 0;
                        push @output, \%hash;
                }
        }
        return @output;
}

#############################################################################

sub revcomp {
  my $nuc = shift @_;
  $nuc =  uc($nuc);
  $nuc =~ tr /ATGCN/TACGN/;
  $nuc = scalar(reverse($nuc));
  return $nuc;
}

sub is_6mer_match {
  my ($srna, $trna) = @_;
  return index( $trna, &revcomp(substr($srna,1,6)) )
}

sub is_7mer_m8_match {
  my ($srna, $trna) = @_;
  return index( $trna, &revcomp(substr($srna,1,7)) )
}

sub is_7mer_A1_match {
  my ($srna, $trna) = @_;
  my $sevenmer_A1_start = index( $trna, &revcomp(substr($srna,1,6)) );
  if ($sevenmer_A1_start != -1 and substr ($trna, $sevenmer_A1_start + 6, 1) eq "A") {
    return $sevenmer_A1_start
  } else {
    return -1
  }
}

sub is_8mer_match {
  my ($srna, $trna) = @_;
  my $eightmer_A1_start = index( $trna, &revcomp(substr($srna,1,7)) );
  if ($eightmer_A1_start != -1 and substr ($trna, $eightmer_A1_start + 7, 1) eq "A") {
    return $eightmer_A1_start
  } else {
    return -1
  }
}

sub is_6mer_off_match {
  my ($srna, $trna) = @_;
  return index( $trna, &revcomp(substr($srna,2,6)) )
}
