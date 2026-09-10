#!/usr/bin/perl 
use strict;
use warnings;
use File::Find;

my$Dir = '/scratch/alpine/fhasan1@xsede.org/sra';

my$BLASTNucleotide = $ARGV[0];
my$Catch = $ARGV[1]-1;
my$NumJobs = $ARGV[2];
open(BLASTN,"$BLASTNucleotide") or die "$!";
my$j = 0;
while(defined(my$line=<BLASTN>)){
	if($j%$NumJobs==$Catch){
		chomp($line);
		my@Array = split("\t",$line);
		print "$Array[0] ";
		my$Output = $Dir ;
		my$OutFile = $Output . '/' . $Array[0] . '_1.fastq';
		if(-e($OutFile)){
			next;
		}
		`fasterq-dump $Array[0] --split-files -O $Output`;
	}
	$j++;
}
