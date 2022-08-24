#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# A small script for splitting qid field in the blast output
#

my $usage=      "USAGE: $0 blast.out --header --help\n".
		"\n".
		"blast.out     : blast result table, last col is assumed to be taxid\n".
		"header        : blast.out has header line\n".
		"help          : print this help\n\n";
my $header = !1;
my $help = !1;
GetOptions('header' => \$header,
	   'help' => \$help) or die $usage;
if($help){die $usage; }
unshift(@ARGV,'-') unless @ARGV;



my $blast_file		= shift(@ARGV);

	
open(IN,"<$blast_file") or die "Can\'t open $blast_file: $!\n";
my $ln	= 0; 
my (@headers,@data,@sp,@sp2);
while(my $l=<IN>){
	@headers=();@data=();
	$ln++;
	chomp($l);
	@sp  = split(/\t/,$l,-1);
	if($ln==1 && $header){
		# first line with headers now in @headers_part2
		shift(@sp);
		my @headers_part2 = @sp;
		
		# read next line
		$l=<IN>;
		$ln++;
		chomp($l);
		@sp 	= split(/\t/,$l,-1);
		@sp2	= split(/_/,$sp[0]);
		for(my $j=0; $j<scalar(@sp2); $j++){
			my @pair = split(/=/,$sp2[$j]);
			push(@headers,$pair[0]);
			push(@data,$pair[1]);
		}
		push(@headers,@headers_part2);
		shift(@sp);
		push(@data,@sp);
		
		# print headers + data
		print "",join("\t",@headers),"\n";
		print "",join("\t",@data),"\n";
		
		next;
	}
	
	
	@sp2		= split(/_/,$sp[0]);	# qseqid
	for(my $j=0; $j<scalar(@sp2); $j++){
		my @pair = split(/=/,$sp2[$j]);
		push(@data,$pair[1]);
	}
	shift(@sp);
	push(@data,@sp);
	print "",join("\t",@data),"\n";
}
