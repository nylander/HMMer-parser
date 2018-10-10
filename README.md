# hmmer-parser.pl

## USAGE:

    ./hmmer-parser.pl [-v] [-s=<E-value|Score|Identity>] [-n <max>] -i hmmer-output-file-from-hmmsearch -o oufile.fasta

## DESCRIPTION:

Parses output from [hmmer](http://hmmer.org/) search.
Prints hits in fasta format with descriptions added to the fasta header (tab separated).

Example:

    >TRINITY_DN52935_c2_g1_i1        Query:COI       Identity:83.35  E-value:1.2e-277        Score:926.9     Result:PRESENT

A label "PRESENT" will be there if identity requirements are fulfilled:

    if ( ($percent_id >= 80) and ($percent_query_coverage >= 80) ) 


## OPTIONS:

    -v           Be verbose (or --noverbose).
    -s=<string>  Sort on either "E-value", "Score", or "Identity". "Score" is default.
    -m, -n=<nr>  Maximum number of hits to show. Default is one.
    -i <infile>  Infile.
    -o <oufile>  Outfile.


## REQUIREMENTS:

BioPerl

## BUGS:


## NOTES:

Tested on output from nhmmer v.3.1b2.

## AUTHOR:

Johan Nylander (JN), Johan.Nylander@nbis.se

## COMPANY:

NBIS/NRM

## VERSION:

1.0

## CREATED:

09/17/2015 10:39:21 PM

## REVISION:

03/24/2017 09:23:52 AM
