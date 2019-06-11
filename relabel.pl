#!/usr/bin/perl -w

unless (@ARGV) {
	print STDERR "\tFirst argument is the file to relabel (case-insenstive) (output going to STDOUT)\n";
	print STDERR "\tSecond argument is a tsv file with existing label in first column, new label in second column\n";
	exit;
}

unless (-e $ARGV[0]) {
	print STDERR "\tCannot find file $ARGV[0]\n";
	exit;
}

unless (-e $ARGV[1]) {
	print STDERR "\tCannot find file $ARGV[1]\n";
	exit;
}

open (INPUT,"< $ARGV[1]");
while (defined($line=<INPUT>)) {
	chomp($line);
	@temp = split(/\t/,$line);
	$relabel{$temp[0]}=$temp[1];
	$relabelcount{$temp[1]}++;
}
close (INPUT);

$fail = 0;
foreach $entry (sort keys %relabelcount) {
	if ($relabelcount{$entry} > 1) {
		print STDERR "ERROR, Non-Unique Label: $entry\n";
		$fail = 1;
	}
}
if ($fail == 1) {
	exit;
}

open (INPUT,"< $ARGV[0]");
while (defined($line=<INPUT>)) {
	chomp($line);
	foreach $entry (keys %relabel) {
		$line =~ s/$entry/$relabel{$entry}/gi;
	}
	print "$line\n";
}	
close (INPUT);

