#!/usr/bin/perl -w

# PubMLST2 analysis, now with generation of a PHYLIP file

unless (@ARGV) {
	print STDERR "\tPubMLST2 will analyze all .fasta files in the current directory.\n";
	print STDERR "\tFirst argument is name of analysis, used in output file names.\n";
	print STDERR "\tSecond argument is the path to the directory containing PubMLST allele and type files\n";
	exit;
}

open (OUTPUT,"> $ARGV[0].mlst");

# gather PubMLST data

chomp(@files=<$ARGV[1]/*>);
foreach $file (@files) {
	open(INPUT,"< $file");
	chomp($line=<INPUT>);			# grab first line of file
	close (INPUT);
	if ($line =~ /^ST/) {						# MLST profiles
		@header = split(/\t/,$line);
		shift(@header);								# remove "ST" at start
		pop(@header);								# remove "clonal_complex" at end
		print OUTPUT "File";
		foreach $allele (@header) {					# alleles as an array
			print OUTPUT "\t$allele";
		}
		print OUTPUT "\tSequence Type\tClonal Complex\n";
		open (INPUT,"< $file");
		while (defined($line=<INPUT>)) {
			chomp($line);
			if ($line =~ /^ST/) { next; }			# header
			@temp = split(/\t/,$line);
			$typeID = shift(@temp);					# extract ST
			undef $typealleles;
			$i = 0;
			foreach $allele (@temp) {				# alleles
				if (defined($header[$i])) {
					$typealleles .= "$header[$i]_$allele";
				} else {
					$MLSTcomplex{$typeID}=$allele;
				}
				$i++;
			}
			$MLSTtype{$typealleles}=$typeID;		# key is e.g. acs_1aro_13gua_12mut_7nuo_72pps_4trp_18, value is ST code
		}
		close (INPUT);
	} elsif ($line =~ /^>/) {					# allele data
		undef $defline;
		undef $sequence;
		open (INPUT,"< $file");
		while (defined($line=<INPUT>)) {
			chomp($line);
			if ($line =~ /^>/) {
				if (defined($defline)) {
					$allele{$defline}=uc($sequence);
				}
				($defline) = $line =~ /^>(.*)/;
				undef $sequence;
			} else {
				$sequence .= $line;
			}
		}
		$allele{$defline}=uc($sequence);
		close (INPUT);
	} else {
		print STDERR "\t$file UNEXPECTED\n";
		exit;
	}
}

# revcomp alleles

foreach $seq (keys %allele) {
	$revallele{$seq} = reversecomplement($allele{$seq});
}

# analyze files

$novelmlstcount = 0;

chomp(@files1=<*.fasta>);
chomp(@files2=<*.fa>);
chomp(@files3=<*.fna>);
undef @files;
foreach $entry (@files1) {
	push(@files,$entry);
}
foreach $entry (@files2) {
	push(@files,$entry);
}
foreach $entry (@files3) {
	push(@files,$entry);
}
foreach $file (@files) {
	print OUTPUT "$file";
	print STDERR "$file\n";
	undef $defline;
	undef $sequence;
	undef %mlsthash;
	undef $typealleles;
	open (INPUT,"< $file");
	while (defined($line=<INPUT>)) {
		chomp($line);
		if ($line =~ /^>/) {
			if (defined($defline)) {
				mlstcheck($sequence);
			}
			($defline) = $line =~ /^>(.*)/;
			undef $sequence;
		} else {
			$sequence .= uc($line);
		}
	}
	mlstcheck($sequence);
	close (INPUT);	
	undef $unresolved;								# tracking type of lack of resolution
	foreach $seq1 (@header) {						# output alleles, allowing for alleles not found & multiple alleles
		undef $status;
		foreach $seq (keys %mlsthash) {
			if ($seq =~ /${seq1}_/) {
				if (defined($status)) {					# check for multiple alleles
					$status .= ",$seq";
					$unresolved = "unresolved allele(s)";
				} else {
					$status = $seq;
				}
			}
		}
		if (defined($status)) {							# allele found
			print OUTPUT "\t$status";
			$typealleles .= $status;					# determine genotype to match to ST
		} else {
			print OUTPUT "\t";									# allele not found
			$unresolved = "unresolved allele(s)";			
		}
	}
	if (defined($typealleles)) {											# allele(s) found
		if (defined($MLSTtype{$typealleles})) {									# type exists
			$finalMLST{$file}=$MLSTtype{$typealleles};
			print OUTPUT "\t$MLSTtype{$typealleles}";
			$typecount{$MLSTtype{$typealleles}}++;
			if (defined($MLSTcomplex{$MLSTtype{$typealleles}})) {					# complex exists
				print OUTPUT "\t$MLSTcomplex{$MLSTtype{$typealleles}}\n";
			} else {																# complex does not exist
				print OUTPUT "\t\n";
			}
		} else {																# type does not exist
			if (defined($unresolved)) {												# does not have all alleles
				print OUTPUT "\t\t\n";
				$typecount{$unresolved}++;
				$finalMLST{$file}="unresolved";		
			} else {																# has all alleles, but not seen before
				$novelmlstcount++;
				$MLSTtype{$typealleles}= "Novel_${novelmlstcount}";
				$typecount{$MLSTtype{$typealleles}}++;
				$finalMLST{$file}=$MLSTtype{$typealleles};
				print OUTPUT "\t$MLSTtype{$typealleles}\t\n";
			}
		}
	} else {																# allele(s) not found
		print OUTPUT "\t\t\n";
		$unresolved = "incorrect pathogen";
		$typecount{$unresolved}++;
		$finalMLST{$file}="incorrect_pathogen";
	}
	foreach $entry (sort keys %mlsthash) {									# generate phylip data
		if (defined($allele{$entry})) {
			$phylip{$file} .= $allele{$entry};
		} else {
			$badphylip{$file}++;
		}
	}
}

close (OUTPUT);

# relabel file

open (OUTPUT2,"> $ARGV[0].mlst.relabel");
foreach $file (keys %finalMLST) {
	print OUTPUT2 "$file\tMLST${finalMLST{$file}}_$file\n";
}
close (OUTPUT2);

# stats

open (OUTPUT,"> $ARGV[0].mlst.log");
print OUTPUT "MLST\tFreq\n";
foreach $typing (sort {$typecount{$b} <=> $typecount{$a}} keys %typecount) {
	print OUTPUT "$typing\t$typecount{$typing}\n";
	print STDERR "$typing\t$typecount{$typing}\n";
}
close (OUTPUT);

# relaxed phylip

foreach $entry (keys %phylip) {
	if (defined($badphylip{$entry})) {
		print STDERR "$file excluded from PHYLIP file due to missing allele sequence(s)\n";
		next;
	}
	$count2 = length($phylip{$entry});
	push(@lengtharray,$count2);
}

@num_sorted = sort {$b <=> $a} @lengtharray;

foreach $entry (sort keys %phylip) {
	if (defined($badphylip{$entry})) {
		next;
	}
	$plength = length($phylip{$entry});
	unless ($plength == $num_sorted[0]) {
		print STDERR "$entry excluded from PHYLIP file due to unexpected allele summed sequence length: $plength\n";
		next;
	}
	$phycount++;
	@name = split(/\./,$entry);
	push(@phyarray,"${name[0]}_MLST_${finalMLST{$entry}} $phylip{$entry}");
}

open (OUTPUT,"> $ARGV[0].mlst.phy");
print OUTPUT "$phycount $num_sorted[0]\n";
foreach $entry (@phyarray) {
	print OUTPUT "$entry\n";
}
close (OUTPUT);

# done

exit;

# subroutines

sub mlstcheck {
	foreach $seq (keys %allele) {
		if ($sequence =~ /$allele{$seq}/) {
			$mlsthash{$seq}++;	
		}
	}
	foreach $seq (keys %revallele) {
		if ($sequence =~ /$revallele{$seq}/) {
			$mlsthash{$seq}++;	
		}
	}
}

sub reversecomplement {

	my($insequence, @sequence, $length, $i, $newsequence, %reverse, $newbase);

	$reverse{"A"}="T";
	$reverse{"C"}="G";
	$reverse{"G"}="C";
	$reverse{"T"}="A";
	$reverse{"R"}="Y";  # R = A or G
	$reverse{"Y"}="R";  # Y = T or C
	$reverse{"M"}="K";  # M = A or C
	$reverse{"K"}="M";  # K = T or G
	$reverse{"S"}="S";  # S = C or G
	$reverse{"W"}="W";  # W = A or T
	$reverse{"H"}="D";  # H = A or C or T
	$reverse{"D"}="H";  # D = T or G or A
	$reverse{"B"}="V";  # B = C or G or T
	$reverse{"V"}="B";  # V = G or C or A
	$reverse{"X"}="N";
	$reverse{"N"}="N";
	
	$insequence = uc($_[0]);
	@sequence=split(//,$insequence);
	$length=@sequence;
	$length=$length-1;
	
	for ($i=$length;  $i >= 0; $i--) {
	    if (defined($reverse{$sequence[$i]})) {
	    	$newsequence .= $reverse{$sequence[$i]};
	    } else {
	    	print STDERR "\tREVCOMP ERROR\tUnrecognized nucleotide \"$sequence[$i]\", using \"N\"\n";	
			$newsequence .= "N";
	    }
	}
	
	return $newsequence;

}

