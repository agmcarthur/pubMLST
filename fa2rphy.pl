#!/usr/bin/perl -w

unless (@ARGV) {
	print STDERR "\tSingle argument is FASTA file to convert to Relaxed PHYLIP format\n";
	exit;
}

unless (-e $ARGV[0]) {
	print STDERR "\tCannot find file $ARGV[0]\n";
	exit;
}

$count = 0;
open (INPUT,"< $ARGV[0]");
while (defined($line=<INPUT>)) {
	chomp($line);
	if ($line =~ /^>/) {
		$count++;
		if (defined($sequence)) {
			$data{$name}=$sequence;
			$seqlength = length($sequence);		
			undef $sequence;
		}
		($name) = $line =~ /^>(.*)/;
	} else {
		$sequence .= $line;
	}
}
close (INPUT);
$data{$name}=$sequence;
$seqlength = length($sequence);		

print "$count $seqlength\n";
foreach $name (keys %data) {
	print "$name $data{$name}\n";
}

print STDERR "$count taxa, $seqlength positions\n";
