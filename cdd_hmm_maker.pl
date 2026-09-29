#!/usr/bin/perl 
use strict;
use warnings;

my$FastqFile = $ARGV[0];

open(MYFILE, $FastqFile) or die "$!";
my$Decider = $ARGV[1] - 1;
my$numJobs = $ARGV[2];
my$i=0;
while(defined(my$line = <MYFILE>)){
	chomp($line);
	if($i%$numJobs==$Decider){
	    my$msa = "cdd_fasta/" . $line;
	    $line =~ s/\.FASTA//;
	    my$outFile = "cdd_hmm/" . $line . ".hmm";
		`hmmbuild $outFile $msa`;
	}
	$i++;
}