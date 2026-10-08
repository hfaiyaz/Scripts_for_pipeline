#!/usr/bin/perl 
use strict;
use warnings;

my$File = $ARGV[0];

open(INPUT, "$File") or die "$!";
open(OUTPUT, ">C_albicans_negative_CircularContigs.fa");
my$Label;
while(defined(my$line = <INPUT>)){
	chomp($line);
	if($line =~ m/\>/){
		$Label = $line;
		next;
	}elsif($line =~ m/\-\-/){
	    next;
	}else{
	    $line =~ s/.+(\:|\-)//;
		$line =~ m/^[acgt]{10}/i;
		my$Begin = $&;
		print "Begin: " . $Begin;
		if($line =~ m/$Begin\w+($Begin\w+)/){
			my$Ending = $1;
			if($line =~ m/^$Ending/){
				print OUTPUT "$Label\n$line\n";
			}
		}
	}
}