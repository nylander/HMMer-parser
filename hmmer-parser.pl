#!/usr/bin/env perl
#===============================================================================

=pod

=head1

         FILE: hmmer-parser.pl

        USAGE: ./hmmer-parser.pl [-v] [-s=<E-value|Score|Identity>] [-n <max>] -i hmmer-output-file-from-hmmsearch -o oufile.fasta

  DESCRIPTION: Parses output from hmmer search. Prints hits in fasta format with descriptions
               added to the fasta header (tab separated).

               Example:

               >TRINITY_DN52935_c2_g1_i1	Query:COI	Identity:83.35	E-value:1.2e-277	Score:926.9	Result:PRESENT

               A label "PRESENT" will be there if identity requirements are fulfilled:

                    if ( ($percent_id >= 80) and ($percent_query_coverage >= 80) ) 

               Tested on output from nhmmer v.3.1b2.

      OPTIONS:
               -v           Be verbose (or --noverbose).
               -s=<string>  sort on either "E-value", "Score", or "Identity". "Score" is default.
               -m, -n=<nr>  Maximum number of hits to show. Default is one.
               -i <infile>  Infile.
               -o <oufile>  Outfile.


 REQUIREMENTS: BioPerl

         BUGS: ---

        NOTES: ---

       AUTHOR: Johan Nylander (JN), Johan.Nylander@nbis.se

      COMPANY: NRM

      VERSION: 1.0

      CREATED: 09/17/2015 10:39:21 PM

     REVISION: 03/24/2017 09:23:52 AM

=cut

#===============================================================================

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::SearchIO;

exec( 'perldoc', $0 ) unless (@ARGV);

my $infile     = q{};
my $outfile    = q{};
my $VERBOSE    = 1;
my $sort       = "Score";
my $max        = 1;
my $wraplength = 80;
my %res_hash   = ();
my $PRINT_FH;
my $using_outfile;


my $r = GetOptions(
    "infile=s"  => \$infile,
    "verbose!"  => \$VERBOSE,
    "outfile=s" => \$outfile,
    "sort=s"    => \$sort,
    "max|n=i"   => \$max,
    "help"      => sub { exec("perldoc", $0); exit(0); },
    );

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

while ( my $result = $in -> next_result ) { # This is a Bio::Search::Result::HMMERResult object

    next unless ($result->query_name);
    print STDERR "Query name: " . $result->query_name . "\n" if ($VERBOSE);

    while( my $hit = $result->next_hit ) {

        while( my $hsp = $hit->next_hsp ) { # A Bio::Search::HSP::HSPI object

            my $percent_id = sprintf("%.2f", $hsp->percent_identity);
            my $percent_query_coverage = sprintf("%.2f", ((($hsp->length('query')/($result->query_length)))*100));
            my $Result = 'ABSENT';
            my $hit_string = $hsp->hit_string;
            $hit_string =~ s/(.{$wraplength})/$1\n/g;

            if ( ($percent_id >= 80) and ($percent_query_coverage >= 80) ) {
                $Result = "PRESENT";
            }
            elsif ( ($percent_id >= 80) and ($percent_query_coverage < 80) ) {
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

__END__
print "\tNumber conserved residues: " . $hsp->num_conserved . "\n";
print "\tQuery length: " . $result->query_length . "\n";
print "\tAlignment length: " . $hsp->length('query') . "\n";
print "\tQuery coverage: $percent_q_coverage\%\n";

"Fraction identical:" , $res_hash{$hitname}{"Fraction_identical"}, "\t",

else {
    if ( ($percent_id >= 80) and ($percent_q_coverage >= 80) ) {
        ## Probably a good one
        print STDOUT ">", $hit->name(), "\t",
                     "Query:", $result->query_name, "\t",
                     "Identity:", $percent_id, "\t",
                     "Fraction identical:" , $frac_identical, "\t",
                     "E-value:", $hsp->evalue(), "\t",
                     "Score:", $hsp->score(), "\t",
                     "Result:PRESENT", "\n",
                     $hsp->hit_string, "\n";
        warn "\n HERE (hit return to continue)\n" and getc();
    }
    elsif ( ($percent_id >= 80) and ($percent_q_coverage < 80) ) {
        print STDERR "\t**", $hit->name(), "Result: TRUNCATED\n\n";
    }
    else {
        print STDERR "\t**", $hit->name(), "Result: ABSENT\n\n";
    }
}

"Fraction_identical" => $frac_identical,

my $frac_identical = sprintf("%.2f", $hsp->frac_identical('total'));

