#!/usr/bin/env perl

# FASTA to Stockholm conversion.
# Reads file names from stdin, prints to stdout.
# Last modified: tis sep 19, 2023  11:11
# Sign: Johan Nylander
# Version: 1.0

use warnings;
use strict;

my $usage = "Usage: $0 <FASTA alignment file(s)>\n";
if (! @ARGV) {
    print $usage;
    exit;
}

my @args;
while (@ARGV) {
    my $arg = shift;
    if ($arg =~ /^-/) {
        if ($arg eq "-h" || $arg eq "--help") {
            print STDERR $usage;
            exit;
        }
        else {
            die $usage;
        }
    }
    else {
        push @args, $arg;
    }
}
push @args, "-" unless @args;

foreach my $fasta (@args) {

    my %seq_hash;
    my @names;
    my $name;
    my $length;
    my $lname;

    open my $FASTA, "<", $fasta or die "Couldn't open '$fasta': $!";
    while (<$FASTA>) {
        if (/^\s*>\s*(\S+)/) {
            $name = $1;
            die "Duplicate name: $name" if defined $seq_hash{$name};
            push @names, $name;
        }
        else {
            if (/\S/ && !defined $name) {
                warn "Ignoring: $_";
            }
            else {
                s/\s//g;
                $seq_hash{$name} .= $_;
            }
        }
    }
    close $FASTA;

    foreach my $name (@names) {
        my $l = length $seq_hash{$name};
        if (defined $length) {
            die "Sequences not all same length ($lname is $length, $name is $l)"
              unless $length == $l;
        }
        else {
            $length = length $seq_hash{$name};
            $lname  = $name;
        }
    }

    print STDOUT "# STOCKHOLM 1.0\n";
    foreach my $name (@names) {
        print STDOUT $name, " ", $seq_hash{$name}, "\n";
    }
    print STDOUT "//\n";
}
exit;

