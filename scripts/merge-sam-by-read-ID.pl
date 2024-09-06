#!/usr/bin/perl

# Andrea Marquard
# Mar 16, 2015
#
#
# Usage: <script.pl> <sam-file-1> <sam-file-2>
#
#
# Aim: take two SAM files, and note the id of the alignments for each read from both files.
#      If a read only aligned to something in one file, it is not reported.
#
# Output: prints the read id, seq id of what the read aligned to in file 1, and in file 2, the read sequence and the read quality string.
# Prints to STDOUT.
#
# Get filenames from command line.
# Required the ID in the first field.
# Requires lines to be pre-sorted by ID!!!  (use --reorder in bowtie - requires that the input to bowtie was also ordered.)
# There should be no SAM headers (begin with @) (use --no-hd in bowtie)
# SAM files should only contain successful alignments (use --no-unal in bowtie)

############### COMMAND LINE ARGUMENTS ###################

my ($f1, $f2) = @ARGV;


############### MAIN PROGRAM ###################

open(IN1, "<", $f1) or die "Cannot open file $f1\n";
open(IN2, "<", $f2) or die "Cannot open file $f2\n";

my $line1;
my $line2;
my $end = 0;

while(defined($line1 = <IN1>) & defined($line2 = <IN2>)) {
	my @l1 = split("\t", $line1);
	my @l2 = split("\t", $line2);
	
	# As long as the read IDs in the two files are not the same,
	# keep reading new lines from the file with the "lower" read ID (alphabetically smaller)
	while($l1[0] ne $l2[0] and $end == 0) {
	  if ($l1[0] lt $l2[0] && defined($line1 = <IN1>)) {
	 	@l1 = split("\t", $line1);
	  } elsif (defined($line2 = <IN2>)) {
	  	@l2 = split("\t", $line2);
	  } else {
        $end = 1;
	  }
	}
	# Now we have a match between read IDs
    unless($end) {
    	print "$l1[0]\t$l1[2]\t$l2[2]\t$l1[9]\t$l1[10]\t$l1[4]\t$l2[4]\n";  # also print the mapping quality
    }
}

