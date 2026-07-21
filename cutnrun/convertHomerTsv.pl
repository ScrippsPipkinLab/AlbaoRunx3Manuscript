#!/usr/bin/env perl
use strict;
use warnings;
use utf8;
use open ':std', ':encoding(UTF-8)';

# homer_to_csv.pl
# Converts HOMER motif TSV -> cleaned CSV with columns:
# TF, FamilyName, Assay, Accession, Database, Consensus, PValue, QValue,
# TargetNumber, TargetPercent, BackgroundNumber, BackgroundPercent
#
# Behavior improvements:
# - FamilyName is extracted ONLY if a parentheses pair appears BEFORE the first '/'
#   (handles cases where assay contains parentheses like GSE accessions).
# - Assay may contain additional '/' characters; everything between the first and
#   last '/' is treated as the Assay.
# - Family strings like "A,B" or "A:B" are normalized to "A:B" (sorted, unique).
#
# Usage:
#   perl homer_to_csv.pl input.tsv output.csv
# or:
#   cat input.tsv | perl homer_to_csv.pl > output.csv

my ($infile, $outfile) = @ARGV;

# open input (or STDIN)
my $IN;
if (defined $infile) {
    open $IN, '<:encoding(UTF-8)', $infile or die "Can't open input '$infile': $!";
} else {
    binmode STDIN, ":encoding(UTF-8)";
    $IN = *STDIN;
}

# open output (or STDOUT)
my $OUT;
if (defined $outfile) {
    open $OUT, '>:encoding(UTF-8)', $outfile or die "Can't open output '$outfile': $!";
} else {
    binmode STDOUT, ":encoding(UTF-8)";
    $OUT = *STDOUT;
}

# -----------------------
# Helpers
# -----------------------
sub trim {
    my ($s) = @_;
    return '' unless defined $s;
    $s =~ s/^\s+|\s+$//g;
    return $s;
}

# CSV-escape fields (quotes doubled, quote whole field if contains comma/quote/newline)
sub csv_escape {
    my ($v) = @_;
    $v = '' unless defined $v;
    $v = trim($v);
    $v =~ s/"/""/g;
    return ($v =~ /[,\"\n]/) ? qq{"$v"} : $v;
}

# Normalize family string: split on comma/colon/semicolon, uniq (case-insensitive), sort alpha, join with colon
sub normalize_family {
    my ($fam) = @_;
    return '' unless defined $fam;
    $fam = trim($fam);
    return '' if $fam eq '';

    my @parts = split /[,:;]/, $fam;
    my %seen;
    my @clean;
    for my $p (@parts) {
        $p = trim($p);
        next if $p eq '';
        $p =~ s/\s+/ /g;
        my $key = lc $p;
        next if $seen{$key}++;
        push @clean, $p;
    }
    @clean = sort { lc($a) cmp lc($b) } @clean;
    return join(':', @clean);
}

# -----------------------
# Read header and map columns (flexible)
# -----------------------
my $header_line = <$IN>;
defined $header_line or die "Input file empty or read error";
chomp $header_line;
$header_line =~ s/^\x{FEFF}//;    # BOM safety

my @hdr = split(/\t/, $header_line, -1);
my %col_idx;
for my $i (0..$#hdr) {
    my $h = $hdr[$i];
    $h = '' unless defined $h;
    $h =~ s/^\s+|\s+$//g;
    my $lh = lc $h;
    if ($lh =~ /motif name/) {
        $col_idx{motif} = $i;
    } elsif ($lh =~ /consensus/) {
        $col_idx{consensus} = $i;
    } elsif ($lh =~ /p-?value/ && !exists $col_idx{pvalue}) {
        $col_idx{pvalue} = $i;
    } elsif ($lh =~ /q-?value|benjamini/) {
        $col_idx{qvalue} = $i;
    } elsif ($lh =~ /#.*target.*sequence|# of target/i) {
        $col_idx{n_target} = $i;
    } elsif ($lh =~ /%.*target/i) {
        $col_idx{pct_target} = $i;
    } elsif ($lh =~ /#.*background.*sequence|# of background/i) {
        $col_idx{n_background} = $i;
    } elsif ($lh =~ /%.*background/i) {
        $col_idx{pct_background} = $i;
    } elsif ($lh =~ /log p-?value/) {
        $col_idx{logp} = $i unless exists $col_idx{logp};
    }
}

# Fallback to expected sample order if matching failed
unless (defined $col_idx{motif} && defined $col_idx{consensus}) {
    @col_idx{qw(motif consensus pvalue logp qvalue n_target pct_target n_background pct_background)} = (0,1,2,3,4,5,6,7,8);
}

# -----------------------
# Output header (user-specified names)
# -----------------------
my @out_hdr = (
    "TF",
    "FamilyName",
    "Assay",
    "Accession",
    "Database",
    "Consensus",
    "PValue",
    "QValue",
    "TargetNumber",
    "TargetPercent",
    "BackgroundNumber",
    "BackgroundPercent",
);
print $OUT join(',', map { csv_escape($_) } @out_hdr) . "\n";

# -----------------------
# Main loop: parse each line
# Parsing strategy (robust to edge cases):
# 1) Find position of first '/' (first_slash).
# 2) If there is a '(... )' pair entirely before first_slash, treat that as TF(Family).
#    Otherwise TF is the substring before first_slash and FamilyName is blank.
# 3) Database is taken from substring after the LAST '/'.
# 4) Assay is everything between first_slash+1 and last_slash-1 (may contain additional '/').
# 5) If the assay ends with '(ACCESSION)', extract ACCESSION as Accession.
# -----------------------
while (my $line = <$IN>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    my @f = split(/\t/, $line, -1);

    my $motif_raw = defined $col_idx{motif} ? ($f[$col_idx{motif}] // '') : '';
    $motif_raw = trim($motif_raw);
    $motif_raw =~ s/^\x{FEFF}//;

    my ($tf, $family, $assay, $accession, $database) = ('','','','','');

    # find first and last slash indexes
    my $first_slash = index($motif_raw, '/');
    my $last_slash  = rindex($motif_raw, '/');

    if ($first_slash >= 0) {
        # candidate TF-part is everything before first slash
        my $tf_part = substr($motif_raw, 0, $first_slash);
        $tf_part = trim($tf_part);

        # check if TF-part contains a parentheses pair like TF(Family)
        if ($tf_part =~ /^(.+?)\s*\(([^)]+)\)\s*$/) {
            # parentheses are within the TF portion -> valid TF(Family)
            $tf     = trim($1);
            $family = trim($2);
        } else {
            # no parentheses in TF portion -> treat TF as the tf_part and family blank
            $tf = $tf_part;
            $family = '';
        }

        # database (after last slash)
        if ($last_slash > $first_slash) {
            $database = trim(substr($motif_raw, $last_slash + 1));
            $assay    = trim(substr($motif_raw, $first_slash + 1, $last_slash - $first_slash - 1));
        } else {
            # only one slash present: assay is everything after first slash
            $assay = trim(substr($motif_raw, $first_slash + 1));
            $database = '';
        }

        # if assay ends with (ACCESSION), extract accession
        if ($assay =~ /^(.*)\s*\(([^()]+)\)\s*$/) {
            $assay     = trim($1);
            $accession = trim($2);
        } else {
            # sometimes accession sits attached to the database part instead
            if ($database =~ /^(.*)\s*\(([^()]+)\)\s*$/) {
                $database = trim($1);
                $accession = trim($2) unless $accession;
            }
        }
    } else {
        # no slash at all; try to parse TF(Family) or treat whole as TF
        if ($motif_raw =~ /^(.+?)\s*\(([^)]+)\)\s*$/) {
            $tf     = trim($1);
            $family = trim($2);
        } else {
            $tf = $motif_raw;
        }
        $assay = '';
        $database = '';
    }

    # Normalize family string if present (A,B -> A:B)
    $family = normalize_family($family) if $family ne '';

    # extract other columns (safe)
    my $consensus  = (defined $col_idx{consensus} && defined $f[ $col_idx{consensus} ]) ? trim($f[ $col_idx{consensus} ]) : '';
    my $pval       = (defined $col_idx{pvalue}   && defined $f[ $col_idx{pvalue}   ]) ? trim($f[ $col_idx{pvalue}   ]) : '';
    my $qval       = (defined $col_idx{qvalue}   && defined $f[ $col_idx{qvalue}   ]) ? trim($f[ $col_idx{qvalue}   ]) : '';

    my $n_target   = (defined $col_idx{n_target}    && defined $f[ $col_idx{n_target}    ]) ? trim($f[ $col_idx{n_target}    ]) : '';
    my $pct_target = (defined $col_idx{pct_target}  && defined $f[ $col_idx{pct_target}  ]) ? trim($f[ $col_idx{pct_target}  ]) : '';
    my $n_bg       = (defined $col_idx{n_background}&& defined $f[ $col_idx{n_background}] ) ? trim($f[ $col_idx{n_background}] ) : '';
    my $pct_bg     = (defined $col_idx{pct_background}&& defined $f[ $col_idx{pct_background}]) ? trim($f[ $col_idx{pct_background}] ) : '';

    # clean numeric/percent fields: remove commas and trailing percent sign
    for my $v_ref (\$n_target, \$n_bg) {
        ${$v_ref} =~ s/,//g;
    }
    for my $v_ref (\$pct_target, \$pct_bg) {
        ${$v_ref} =~ s/%$//;
        ${$v_ref} =~ s/,//g;
    }

    # print CSV row in requested order
    my @out = (
        $tf,
        $family,
        $assay,
        $accession,
        $database,
        $consensus,
        $pval,
        $qval,
        $n_target,
        $pct_target,
        $n_bg,
        $pct_bg,
    );
    print $OUT join(',', map { csv_escape($_) } @out) . "\n";
}

close $IN if defined $infile;
close $OUT if defined $outfile;
__END__
