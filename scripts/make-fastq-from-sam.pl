#!/usr/bin/perl

# Andrea Marquard
# Mar 16, 2015
#
#
# Usage: <script.pl> <sam-file-1> <sam-file-2> ... <sam-file-n> <threshold>
#
#
# Aim: take many SAM files (eg. output from bowtie) and keep only those reads that aligned in (as a minimum) a specified number of files.
# Output: prints the output in FASTQ format. Prints to STDOUT.
#
# Files are assumed to have results for exactly the same reads, in the same order (Use -reorder in bowtie to ensure this).
# Files should therefore obviously have both aligned and unaligned reads (ie. DON'T use --no-unal in bowtie).
# SAM files should NOT have headers (begin with @) (use --no-hd in bowtie)
#
# Get filenames from command line, followed by the minimum number of files to be found in for a read to be kept.
#


############### COMMAND LINE ARGUMENTS ###################

# Get the FASTX format to print:
my $fastx_format = pop @ARGV;
chomp $fastx_format;

# Get the minimum number of successful alignments
my $threshold = pop @ARGV;

# The rest of ARGV is the file names
my @Files = @ARGV;



############### MAIN PROGRAM ###################

my $nfiles = scalar @Files;
my @Handles = ();

# Open all files
foreach my $file (@Files) {
	open(my $fh, "<", $file) or die "Cannot open file $file\n";
	push(@Handles, $fh);
}

# Read the next four lines from all filehandles
while(defined(my $Lines = ReadManyFiles(@Handles))) {
	# Reset match counter to max number of files
	my $matches = $nfiles;

    # Go through files, decrease match counter if this read failed to align in a file
	foreach my $line (@{$Lines}) {
		my @array = split("\t", $line);
		$matches-- if $array[2] eq "*";
	}
	
	# Only print the FASTQ entry if the match count is at least the threshold
    if($matches >= $threshold) {
    	my @array = split("\t", ${$Lines}[0]);
    	# The read id, sequence and quality is the same no matter which file you look in.

    	if ($fastx_format eq 'fastq') {
	    	print "\@$array[0]\n";
	    	print "$array[9]\n";
	    	print "+\n";
	    	print "$array[10]\n";
	    }
	    if ($fastx_format eq 'fasta') {
	    	print ">$array[0]\n";
	    	print "$array[9]\n";
	    }

    }
}

# Close file handles
foreach my $handle (@Handles) {
	close $handle;
}




############### SUB-ROUTINES ###################

sub ReadManyFiles {
	my @Hands = @_;
	my @Lines = ();
	foreach my $handle (@Hands) {
		if(defined(my $line = <$handle>)) {
			push(@Lines, $line);
		} else {
			return(undef);
		}
	}
	return(\@Lines);
}