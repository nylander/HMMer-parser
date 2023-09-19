#!/usr/bin/env perl
#===============================================================================

=pod

=head1

         FILE: hmmer-parser.pl

        USAGE: ./hmmer-parser.pl [-v] [-s=<E-value|Score|Identity>] [-n <max>] -i hmmer-output-file-from-hmmsearch -o oufile.fasta

  DESCRIPTION: Parses output from HMMEer searches (hmmsearch or nhmmer).
               Prints hits in fasta format with descriptions added to the fasta header (tab separated).

               Example output fasta header:

               >TRINITY_DN52935_c2_g1_i1	Query:COI	Identity:83.35	E-value:1.2e-277	Score:926.9	Result:PRESENT

               A label "PRESENT" will be there if identity requirements are fulfilled:

                    if ( (percent identity in HSP >= PERCENTAGE) and (percent coverage of HSP to query >= COVERAGE) )

               If not "PRESENT", then the tag can be labeled as "ABSENT" or "TRUNCATED" depending on the values of
               percent identity in HSP and the percent coverage of HSP to query.

               The values of PERCENTAGE and COVERAGE can be set by options -p and -c and will only affect
               the tag "Result" in the output.

      OPTIONS:
               -v           Be verbose (or --noverbose).
               -s=<string>  sort on either "E-value", "Score", or "Identity". "Score" is default.
               -m, -n=<nr>  Maximum number of hits to show. Default is "1".
               -p=<integer> Minimum percentage for residual identity in alignment. Default is "80".
               -c=<integer> Minimum coverage ((length of query in alignment pair / original length of query) * 100).
                            Default is "80".
               -i <infile>  Infile.
               -o <oufile>  Outfile.


 REQUIREMENTS: BioPerl, Bio::SearchIO::hmmer

       AUTHOR: Johan Nylander

      COMPANY: NRM

      VERSION: 1.0.1

      CREATED: 09/17/2015

     REVISION: tis 19 sep 2023 15:37:16

=cut

#===============================================================================

use strict;
use warnings;
use Getopt::Long;
use Bio::SearchIO;

exec( 'perldoc', $0 ) unless (@ARGV);

my $infile             = q{};
my $outfile            = q{};
my $VERBOSE            = 1;
my $sort               = "Score";
my $max                = 1;
my $wraplength         = 80;
my $coverage           = q{};
my $coverage_default   = 80;
my $percentage         = q{};
my $percentage_default = 80;
my %res_hash           = ();
my $PRINT_FH;
my $using_outfile;

GetOptions(
    "infile=s"     => \$infile,
    "verbose!"     => \$VERBOSE,
    "outfile=s"    => \$outfile,
    "sort=s"       => \$sort,
    "max|n=i"      => \$max,
    "coverage=i"   => \$coverage,
    "percentage=i" => \$percentage,
    "help"         => sub { exec("perldoc", $0); exit(0); },
);

if (! $coverage) {
    $coverage = $coverage_default;
}
if (! $percentage) {
    $percentage = $percentage_default;
}

if (! $infile) {
    die "Error: need an infile.\nSee $0 -h for usage\n.";
}
if ($outfile) {
    open ($PRINT_FH, '>', $outfile) or die "$0 : Failed to open output file $outfile : $!\n\n";
    $using_outfile = 1;
}
else {
    $PRINT_FH = *STDOUT; # Using the typeglob notation in order to use STDOUT as a variable
}

my $in = Bio::SearchIO->new( -file => $infile, -format => 'hmmer' );

print STDERR "Reading file: $infile\n" if ($VERBOSE);

while ( my $result = $in -> next_result ) { # A Bio::Search::Result::HMMERResult object

    next unless ($result->query_name);
    print STDERR "Query name: " . $result->query_name . "\n" if ($VERBOSE);

    while( my $hit = $result->next_hit ) {

        while( my $hsp = $hit->next_hsp ) { # A Bio::Search::HSP::HSPI object
            my $percent_id = sprintf("%.2f", $hsp->percent_identity);

            my $percent_query_coverage = sprintf("%.2f", ((($hsp->length('query')/($result->query_length)))*100));

            my $hit_string = $hsp->hit_string;
            $hit_string =~ s/(.{$wraplength})/$1\n/g; # Sequence in HSP, wrapped

            my $Result = 'ABSENT';

            if ( ($percent_id >= $percentage) and ($percent_query_coverage >= $coverage) ) {
                $Result = "PRESENT";
            }
            elsif ( ($percent_id >= $percentage) and ($percent_query_coverage < $coverage) ) {
                $Result = "TRUNCATED";
            }

            if ($sort) {
                $res_hash{$hit->name()} = {
                    "Query" => $result->query_name,
                    "Identity" => $percent_id,
                    "E-value" => $hsp->evalue(),
                    "Score" => $hsp->score(),
                    "Result" => $Result,
                    "hit_string" => $hit_string,
                };
            }
        }
    }
}

if (! %res_hash) {
    die "No hits parsed\n";
}

if ($sort) {

    my @sorted_hitnames = ();

    if ($sort =~ /^E/i) {
        @sorted_hitnames = sort { $res_hash{$a}{"E-value"} <=> $res_hash{$b}{"E-value"} }  keys %res_hash;
    }
    elsif ($sort =~ /^I/i) {
        @sorted_hitnames = sort { $res_hash{$b}{"Identity"} <=> $res_hash{$a}{"Identity"} }  keys %res_hash;
    }
    elsif ($sort =~ /^S/i) {
        @sorted_hitnames = sort { $res_hash{$b}{"Score"} <=> $res_hash{$a}{"Score"} }  keys %res_hash;
    }
    else {
        die "Unknown sort argument for sort: $sort. Use Score, E-value, or Identity \n";
    }

    if ($max > scalar(@sorted_hitnames)) {
        $max = scalar(@sorted_hitnames);
        print STDERR "Max value adjusted. Only $max hits found.\n" if ($VERBOSE);
    }

    if ($max <= scalar(@sorted_hitnames)) {
        my $i = 0;
        foreach my $hitname (@sorted_hitnames) {
            last if ($i == $max);
            print $PRINT_FH ">", $hitname, "\t",
                  "Query:", $res_hash{$hitname}{"Query"}, "\t",
                  "Identity:", $res_hash{$hitname}{"Identity"}, "\t",
                  "E-value:", $res_hash{$hitname}{"E-value"}, "\t",
                  "Score:", $res_hash{$hitname}{"Score"}, "\t",
                  "Result:",$res_hash{$hitname}{"Result"}, "\n",
                  $res_hash{$hitname}{"hit_string"}, "\n";
            $i++;
        }
    }
}

if ($using_outfile) {

    close($PRINT_FH);

    if (-e $outfile) {
        print STDERR "Wrote file: $outfile\n" if ($VERBOSE);
    }
    else {
        die "Could not create outfile: $outfile $! \n";
    }
}

