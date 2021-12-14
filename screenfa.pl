#!/usr/bin/perl -w

# report FASTA using other than AGCTUN nucleotide codes (not liked by parsnp/MUSCLE)

chomp(@file1=<*.fna>);
chomp(@file2=<*.fa>);
chomp(@file3=<*.fasta>);

foreach $entry (@file1) {
	push(@files,$entry);
}
foreach $entry (@file2) {
	push(@files,$entry);
}
foreach $entry (@file3) {
	push(@files,$entry);
}

print "Characters other than AGCTUN:\n";
foreach $entry (@files) {
	$error = 0;
	open (INPUT,"< $entry");
	while (defined($line=<INPUT>)) {
		chomp($line);
		if ($line =~ /^>/) {
			next;
		}
		$line2 = uc($line);
		if ($line2 =~ /[BDEFHIJKLMOPQRSVWXYZ]/) {
			$error = 1;
		}
	}
	close (INPUT);
	if ($error == 1) {
		print "$entry\n";
	}
}
