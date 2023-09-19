# hmmer-parser.pl

## USAGE

    ./hmmer-parser.pl [options] -i output-file-from-hmmsearch-or-nhmmer

## DESCRIPTION

Parses output from [hmmer](http://hmmer.org/) searches (hmmsearch or nhmmer).

Prints hits in fasta format with descriptions added to the fasta header (tab separated).

Will print to stdout or to file.

Prints hits in fasta format with descriptions added to the fasta header (tab separated).

Example output fasta header:

    >TRINITY_DN52935_c2_g1_i1        Query:COI       Identity:83.35  E-value:1.2e-277        Score:926.9     Result:PRESENT

A label "PRESENT" will be there if identity requirements are fulfilled:

    if ( (percent identity in HSP >= PERCENTAGE) and (percent coverage of HSP to query >= COVERAGE) )

If not "PRESENT", then the tag can be labeled as "ABSENT" or "TRUNCATED" depending on the
values of percent identity in HSP and the percent coverage of HSP to query.

The values of PERCENTAGE and COVERAGE can be set by options -p and -c and will only affect
the tag "Result" in the output fasta headers.

## OPTIONS

    -i <infile>  Infile. Mandatory.
    -m, -n=<nr>  Maximum number of hits to show.
                 Default is "1".
    -s=<string>  Sort output sequences on either "E-value", "Score", or "Identity".
                 "Score" is default.
    -p=<integer> Minimum percentage for residual identity in alignment.
                 Default is "80".
    -c=<integer> Minimum coverage ((length of query in alignment pair/original length of query)*100).
                 Default is "80".
    -o <oufile>  Outfile.
    -v           Be verbose (or --noverbose).

## REQUIREMENTS

BioPerl, Bio::SearchIO::hmmer.

Example installation on Ubuntu 22.04:

    $ sudo apt install \
          libbio-perl-perl \
          libbio-perl-run-perl \
          libbio-searchio-hmmer-perl

## WORKED EXAMPLES

### hmmsearch (HMMER 3.3.2); one coi sequence against 24 coi sequences

    $ scripts/fas2sto.pl data/ref-coi.fas > data/ref-coi.sto
    $ hmmbuild --cpu 4 data/ref-coi.hmm data/ref-coi.sto
    $ hmmsearch --cpu 4 data/ref-coi.hmm \
        data/coi.fas > data/coi-vs-ref-coi.hmmsearch.out

Parse

    $ ./hmmer-parser.pl \
        -i data/coi-vs-ref-coi.hmmsearch.out \
        -o data/coi-vs-ref-coi.hmmsearch.hmm-parser.fas

### nhmmer (HMMER 3.3.2); coi HMM-profile (calculated from a coi multiple sequence alignment) against one genome (nt)

    $ wget -O - \
       "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/335/GCF_000002335.3_Tcas5.2/GCF_000002335.3_Tcas5.2_genomic.fna.gz" | \
       gunzip -c > data/GCF_000002335.3.fna
    $ scripts/fas2sto.pl data/ref-coi.fas > data/ref-coi.sto
    $ hmmbuild --cpu 4 data/ref-coi.hmm data/ref-coi.sto
    $ nhmmer --cpu 4 \
        -o data/GCF_000002335.ref-coi-vs-GCF_000002335.3.nhmmer.out \
        data/ref-coi.hmm \
        data/GCF_000002335.3.fna

Parse

    $ ./hmmer-parser.pl \
        -i data/GCF_000002335.ref-coi-vs-GCF_000002335.3.nhmmer.out \
        -o data/GCF_000002335.ref-coi-vs-GCF_000002335.3.nhmmer.hmm-parser.fas

## NOTES

Tested on output from hmmsearch and nhmmer from HMMer v.3.1b2 and v.3.3.2.
Beware of change in output format between HMMer versions.

## AUTHOR

Johan Nylander

## COMPANY

NRM

## LICENSE

MIT. See [LICENSE file](LICENSE)

## DOWNLOAD

<https://github.com/nylander/HMMer-parser>
