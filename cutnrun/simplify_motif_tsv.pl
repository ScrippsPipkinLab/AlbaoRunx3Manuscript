#!/usr/bin/env perl
use strict;
use warnings;

# ---------------------------------------------------------------------------
# simplify_motif_tsv.pl
#
# Reformats a HOMER-style motif TSV file:
#   - Drops the header line
#   - Truncates "Motif Name" to just the part before the first "/"
#   - Reorders columns to:
#       PositionID  Sequence  Motif Name  Offset  Strand  MotifScore
#
# Usage:
#   perl simplify_motif_tsv.pl input.tsv > output.tsv
#   perl simplify_motif_tsv.pl file1.tsv file2.tsv ... > output.tsv
#   perl simplify_motif_tsv.pl dir/*.instances > output.tsv
# ---------------------------------------------------------------------------

# Expect one or more file paths (e.g. from a shell glob)
die "Usage: $0 <input.tsv> [more.tsv ...]\n" unless @ARGV;

# Process every file given, stripping the header line from EACH one
foreach my $infile (@ARGV) {
    open(my $fh, '<', $infile)
        or die "Could not open '$infile': $!\n";

    # Read and discard this file's header line
    my $header = <$fh>;

    while (my $line = <$fh>) {
        chomp $line;
        next unless length $line;  # skip blank lines

        # Original column order:
        # PositionID, Offset, Sequence, Motif Name, Strand, MotifScore
        my ($position_id, $offset, $sequence, $motif_name, $strand, $motif_score) =
            split(/\t/, $line);

        # Isolate the motif name before the first "/"
        my ($short_motif_name) = split(/\//, $motif_name, 2);

        # Print in the new column order:
        # PositionID, Sequence, Motif Name, Offset, Strand, MotifScore
        print join("\t",
            $position_id,
            join(",",
                $short_motif_name,
                $sequence,
                $offset,
                $strand,
                $motif_score
            )
        ), "\n";
    }

    close($fh);
}