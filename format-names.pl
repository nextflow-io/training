#!/usr/bin/perl
use warnings;
use strict;

die "Please specify (1) Input file (2) Post fix name\n" unless(@ARGV==2);

my $file= $ARGV[0];
my $name= $ARGV[1];
my $outfile="$file\.$name";

open(my $emailfile, "<", $file)   or die "Could not open $file \n";
open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";

my %unique;

while (my $line=<$emailfile>){
	chomp $line;
	my @chars = split("@", $line);
	my $newname = $chars[0];
	$newname =~ s/\./-/g;    #Remove dots from name
	$newname =~ s/ /-/g;     #Remove spaces from name
	my $outname = "$newname\_$name";      #Add on session name 
	if ($unique{$outname}){
		#Do nothing, it is a replicate
	}
	else{
		print $outhandle "$line\t$outname\n";
		print "$outname\n";
		$unique{$outname}="Yes";    #Put name in checking hash
	}
	
}

print "For email and Id codes, see $outfile\n";
