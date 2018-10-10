#!/usr/bin/perl 
#===============================================================================

=pod

=head1

         FILE: parse_nhmmer.pl

        USAGE: ./parse_nhmmer.pl -i nhmmer.out -f genome.fas

  DESCRIPTION: Parse HMMER output from command:

               $ nhmmer --tblout nhmmer.out --cpu 10 markers.hmm genome.fas > log.nhmmer.log 2> /dev/null &

               Input:
               HMMER-output file (from command above), Genome-fasta-file

               Expected fasta header format for this version (will use "1000"):
               1000.sate.default.pep2cds.removed.shortname.filtered
               
               Output:
               Will write separate fasta files to output directory.

      OPTIONS: -i <HMMER-output>
               -f <Genome-fasta-file>
               -d <Output-directory>

 REQUIREMENTS: ---

         BUGS: ---

        NOTES: ---

       AUTHOR: Johan Nylander (JN), Johan.Nylander@bils.se

      COMPANY: BILS/NRM

      VERSION: 1.0

      CREATED: 09/17/2015 11:06:04 PM

     REVISION: 09/21/2015 01:15:21 PM

=cut

#===============================================================================

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

exec( 'perldoc', $0 ) unless (@ARGV);

my %HoH         = (); # key:scaffold header, val:hash:key:query name, vals:
my %query_hash  = (); # key: name, val: count
my $infile      = q{};
my $fasta       = q{};
my $outfile     = q{};
my $directory   = 'Parse-HMMer-output';
my $VERBOSE     = 1;
my $r           = GetOptions(
                    "directory=s"   => \$directory,
                    "fasta=s"       => \$fasta,
                    "help"          => sub { exec( "perldoc", $0 ); exit(0); },
                    "infile=s"      => \$infile,
                    "outfile=s"     => \$outfile,
                    "verbose!"      => \$VERBOSE,
                    );

## Read HMM-output file
open my $INFILE, "<", $infile or die "$!\n";
print STDERR "Parsing HMMER output file $infile\n" if ($VERBOSE);
while (<$INFILE>) {
    chomp;
    next if /^#/;
    my (
        $target_name,     $target_accession, $query_name,
        $query_accession, $hmmfrom,          $hmm_to,
        $alifrom,         $ali_to,           $envfrom,
        $env_to,          $sq_len,           $strand,
        $E_value,         $score,            $bias,
        $description_of_target
    ) = split /\s+/, $_;

    ## Find first occurence of each query_name. This would be the best hit (assuming no ties)
    if ( $query_hash{$query_name}++ ) {
        next;
    }
    else {
        my $details = "/target name=$target_name /accession=$target_accession /query name=$query_name /accession=$query_accession /hmmfrom=$hmmfrom /hmm to=$hmm_to /alifrom=$alifrom /ali to=$ali_to /envfrom=$envfrom /env to=$env_to /sq len=$sq_len /strand=$strand /E-value=$E_value /score=$score /bias=$bias /description of target=$description_of_target";
        my $scaffold = $target_name;
        $HoH{$scaffold}{$query_name}{'alifrom'} = $alifrom;
        $HoH{$scaffold}{$query_name}{'ali to'}  = $ali_to;
        $HoH{$scaffold}{$query_name}{'strand'}  = $strand;
        $HoH{$scaffold}{$query_name}{'details'} = $details;
    }
}
close($INFILE);

## Create output directory for sequences
if ( -e $directory ) {
    die "Warning. Directory $directory already exists.\n";
}
else {
    unless ( mkdir $directory ) {
        die "Unable to create $directory\n";
    }
}
print STDOUT "Created output directory $directory\n" if ($VERBOSE);

## Read genome file
open my $FASTA, "<", $fasta or die "Cannot open fasta file: $fasta $! \n";
print STDERR "Reading fasta file $fasta\n" if ($VERBOSE);
my $def = $/;
$/ = ">";
while (<$FASTA>) {
    next if ( $_ eq '' );
    my ( $id, @sequencelines ) = split /\n/;
    $id =~ s/>//;
    if ( exists( $HoH{$id} ) ) {
        my $sequence = '';
        foreach my $line (@sequencelines) {
            $sequence .= $line;
        }
        foreach my $key ( keys %{ $HoH{$id} } ) {
            my ( $geneid, @rest ) = split /\./, $key; # highly format dependent!
            open my $OUTFILE, ">", "$directory/$geneid.fas"
              or die "Could not open fasta file for writing: $! \n";
            print STDERR "    writing $geneid.fas\n" if ($VERBOSE);
            print $OUTFILE ">$geneid $HoH{$id}{$key}{'details'}\n";
            my ( $from, $to ) = q{};
            if ( $HoH{$id}{$key}{'strand'} eq '+' ) {
                $from = $HoH{$id}{$key}{'alifrom'} - 1;
                $to   = $HoH{$id}{$key}{'ali to'};
            }
            elsif ( $HoH{$id}{$key}{'strand'} eq '-' ) {
                $from = $HoH{$id}{$key}{'ali to'} - 1;
                $to   = $HoH{$id}{$key}{'alifrom'};
            }
            my $len = $to - $from;
            my $seq = substr( $sequence, $from, $len );
            if ( $HoH{$id}{$key}{'strand'} eq '-' ) {
                $seq = revcomp($seq);
            }
            $seq =~ s/(.{60})/$1\n/g;
            print $OUTFILE $seq, "\n";
            close($OUTFILE);
        }
    }
}
close($FASTA);
$/ = $def;
print STDERR "\nEnd of script\n" if ($VERBOSE);

sub revcomp {
    ## Warning: does not handle IUPAC
    my $seq = shift;
    $seq = uc($seq);
    $seq = scalar reverse $seq;
    $seq =~ tr/GATC/CTAG/;
    return $seq;
}
__END__
