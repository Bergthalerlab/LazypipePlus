#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# Prints only the top scoring query-subject entry for each query in a given SAM file
# Entries are ranked based on BWA score in column that is at position >=12 and has format AS:i:num
# Enties with no AS:i:score are ignored


my $usage=      "\nUSAGE: $0 file.sam\n".
		"\n".
		"file.sam     : alignments in SAM format sorted by QNAME\n".
		"		when ommitted reads from STDIN\n".
		"example use:\n".
		"samtools view -F4 file.bam | $0 > file.flt.sam\n".
		"\n";

my $help = !1;
GetOptions('help|h' => \$help) or die $usage;
if($help){ die $usage; }
unshift(@ARGV,'-') unless @ARGV;

my $sam_file = shift(@ARGV);
open(IN,"<$sam_file") or die "Can\'t open $sam_file: $!\n";

my $l;
my $li 		= 0;
my $score	= -1;
my $best_score	= -1;
my $best_line	= !1;
my $qname	= "";
my $first_ali	= 1;
my @sp;

#
# The column to which bwa mem prints AS:i:score varies 
# and can be basically any column starting from SAM OPT-field, i.e. any column starting from 12.
# Thus to make this bulletproof we search columns 12>last for each line

while($l=<IN>){
	$li++;
	if( $l =~ m/^@/){	# SAM HEADER SECTION
		print "$l";
		next;
	}
	
	# SAM ALIGNMENT SECTION
	chomp($l);
	@sp 			= split(/\t/,$l,-1);
	$score			= -1;
	for(my $col=12-1; $col<scalar(@sp); $col++){
		if( $sp[$col] =~ m/AS\:i\:([\d]+)/ ){
			$score  = $1; 
			last;
		}
	}
	if($score< 0){
		next;
	}
	
	#print STDERR "# found AS:i:score on line $li\n";
	
	if($first_ali){
		$qname		= $sp[0];
		$best_line 	= $l;
		$best_score	= $score;
		$first_ali	= !1;
	}
	elsif($sp[0] ne $qname){
		print "$best_line\n";
		$qname = $sp[0];
		$best_line = $l;
		$best_score= $score;
	}
	else{
		if($score > $best_score){
			$best_score 	= $score;
			$best_line 	= $l;
		}
	}
}
print "$best_line\n";

