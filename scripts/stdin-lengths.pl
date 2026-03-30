#!/usr/bin/perl

# Andrea Marquard
# May 5, 2015
#

use strict;

my %LEN = ();

while (defined(my $seq = <STDIN>) ) {
	chomp $seq;
	my $l = length($seq);
	if (exists($LEN{$l})) {
		$LEN{$l}++;
	} else {
		$LEN{$l} = 1;
	}
}

print "length\treads\n";
foreach my $l (sort { $a <=> $b } keys %LEN) {
	print "$l\t$LEN{$l}\n";
}



