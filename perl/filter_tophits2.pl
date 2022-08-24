#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);


# Prints only the top scoring query-subject entry for each query in a given homology search result table.
# This can be e.g. BLAST/SANS/Centrifuge results.
# Entries are ranked based on --bitscol column.
# Comment lines are ignored


my $usage=      "USAGE: $0 --qcol int --bitscol int [--ties -h] db_hits.tab\n".
		"\n".
		"qcol        : query column\n".
		"bitscol     : bitscore column\n".
		"ties        : retain all top-scoring ties [default:NO]\n".
		"db_hits.tab : file to process (if ommited reads STDIN).\n\n".
		"#@ Comment lines are retained\n\n";
		
my $qcol	= 0;
my $bitscol	= 0;
my $retain_ties	= !1;
my $help 	= !1;

GetOptions(	'qcol=i'	=> \$qcol,
		'bitscol=i'	=> \$bitscol,
		'ties'		=> \$retain_ties,
		'help|h' 	=> \$help,) or die $usage;
$qcol -= 1;
$bitscol -=1;
if($help){ die $usage; }
if( ($qcol<0) || ($bitscol<0)){
	print STDERR "ERROR: invalid arguments --qcol $qcol --bitscol $bitscol\n";
	die $usage;
}
unshift(@ARGV,'-') unless @ARGV;

my $file = shift(@ARGV);
open(IN,"<$file") or die "Can\'t open $file: $!\n";
my $l;
my $li 		= 0;
my $score	= -1;
my $best_score	= -1;
my $best_line	= !1;
my $qname	= "";
my $first_ali	= 1;
my @sp;
while($l=<IN>){
	$li++;
	if( $l =~ m/^[#@]/){	# comment lines
		print "$l";
		next;
	}
	
	chomp($l);
	@sp 	= split(/\t/,$l,-1);
	$score	= $sp[$bitscol];
	if($first_ali){
		$qname		= $sp[$qcol];
		$best_line 	= $l;
		$best_score	= $score;
		$first_ali	= !1;
	}
	elsif($sp[$qcol] ne $qname){
		print "$best_line\n";
		$qname = $sp[$qcol];
		$best_line = $l;
		$best_score= $score;
	}
	else{	
		if($score > $best_score){
			$best_score 	= $score;
			$best_line 	= $l;
		}
		elsif(($score == $best_score) && $retain_ties){
			print "$best_line\n";
			$best_line 	= $l;
		}
	}
}
print "$best_line\n";

