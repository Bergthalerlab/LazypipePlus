#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# Converts Centrifuge summary output file to taxonomy profile
# Tested for Centrifuge version 1.0.3
#
# input: nodes.dmp from NCBI taxonomy
#	 names.dmp from NCBI taxonomy containing taxnames output by metaphlan2
#	 clade summary report by: centrifuge -S file
# 
#
# output: taxonomy profile in CAMI format (see https://github.com/bioboxes/rfc/blob/master/data-format/profiling.mkd)
#

my $usage=      "USAGE: $0 nodes.dmp names.dmp centrifuge.summary 1> cami_taxonomy_profile\n".
		"\n".
		"nodes.dmp \t\t: node file from NCBI taxonomy\n".
                "names.dmp \t\t: name file from NCBI taxonomy\n".
		"centrifuge.summary\t\t: centrifuge -S file\n".
		"cami_taxonomy_profile\t: taxonomy profile in CAMI format\n\n".
		"NCBI taxonomy *.dmp files:\n".
		"ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz\n".
		"\n";


if(scalar @ARGV<3) { die $usage; }

my $node_dmp_file	= shift(@ARGV);
my $name_dmp_file	= shift(@ARGV);
my $summary_file	= shift(@ARGV);


# read taxonomy nodes
print STDERR "# reading $node_dmp_file\n";
my %taxnodes;
open(IN,"<$node_dmp_file") or die "Can\'t open $node_dmp_file: $!\n";
my $ln=0;
while(my $l=<IN>){
	$ln++;
	chomp($l);
	my @sp= split(/\t\|\t/,$l);
	if(scalar(@sp)<3){
		die "ERROR: missing data at $ln:$node_dmp_file\n";
	}
	$taxnodes{$sp[0]}= join("\t",$sp[0],$sp[1],$sp[2]);
}

# read taxonomy node names
print STDERR "# reading $name_dmp_file\n";
my %taxnames;
open(IN,"<$name_dmp_file") or die "Can\'t open $name_dmp_file: $!\n";
$ln=0;
while(my $l=<IN>){
	$ln++;
	chomp($l);
	my @sp= split(/\t\|\t/,$l);
	for(my $j=0; $j<scalar(@sp); $j++){
		$sp[$j] =~ s/^(\|\t)|(\t\|)$//g; # removing leading and trailing \t|\t characters
	}
	if(scalar(@sp)<4){
		die "ERROR: missing data at $ln:$name_dmp_file\n";
	}
	if( !defined($taxnames{$sp[0]}) ){ # if not defined we add any name available
		$taxnames{$sp[0]}= $sp[1];
	}
	elsif( lc($sp[3]) eq 'scientific name' ){ # else we replace saved name with scientific
		$taxnames{$sp[0]}= $sp[1];
	}
}


my %rankcode_rank = 
	('D' => 'superkingdom',
	 'P' => 'phylum',
	 'C' => 'class',
	 'O' => 'order',
	 'F' => 'family',
	 'G' => 'genus',
	 'S' => 'species');
# root/domain/unclassified and minor ranks are ignored: add here to include
	 #'R' => 'root',
	 #'U' => 'unclassified');
my @rank_list = ('superkingdom','phylum','class','order','family','genus','species');
my %rank_set;
for my $rank(@rank_list){ $rank_set{$rank} = 1;}


print '@',"Version:0.9.1\n";
print '@',"SampleID:$kraken_file\n";
print '@',"Ranks:",join('|',@rank_list),"\n";
print '@',"TaxonomyID:NCBI taxonomy\n";
print '@@',"TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\t_kraken2_fragment_num\n";


open(IN,"<$kraken_file") or die "Can\'t open $kraken_file: $!\n";
$ln=0;
while(my $l=<IN>){
	$ln++;
	chomp($l);
	if($l =~ m/^#/){ next;}
	my ($percentage,$count,$count2,$rankcode,$taxid,$taxname)= split(/\t/,$l,-1);

	if( !defined($rankcode_rank{$rankcode})){ # ignoring root/domain/unclassified and minor ranks (e.g. K1,F2..)
		next;
	}
	
	my $rank	= $rankcode_rank{$rankcode};
	my @taxpath	= reverse get_parent_taxids($taxid,\%taxnodes,\%rank_set, 1);
	my @taxpathsn;
	for my $t(@taxpath){
		my $name = defined($taxnames{$t})? $taxnames{$t} : "";
		push(@taxpathsn,$name);
	}

	
	print join("\t",$taxid,$rank,join('|',@taxpath), join('|',@taxpathsn), $percentage, $count),"\n";
}


# @parent_taxids= get_parent_taxids($taxid,\%taxnodes,\%parent_rank_set);
# 
# Returns parent taxa (including the node itself), limited by the list of ranks
#
# Input:
# $taxid		: taxnode id
# \%taxnodes		: pointer to the %taxnodes hash
# \%parent_rank_set	: set of ranks to report: e.g. {'family' =>1,'genus'=>1,'species'=>1}
#
# Output
# @parent_taxids	: list of parent nodes of $taxid, limited to ranks in @parent_rank_list
		
sub get_parent_taxids{
	my $taxid		= shift(@_);
	my $taxnode_ref		= shift(@_);
	my %parent_rank_set	= %{shift(@_)};
	my @parent_taxids;

	push(@parent_taxids,$taxid);
	while( defined($taxnode_ref->{$taxid}) ){
		my $node	= $taxnode_ref->{$taxid};
		my @sp		= split(/\t/,$node,-1);
		if($sp[0] == $sp[1]){ 		# node is it's own parent
			last;
		}
		$taxid		= $sp[1];	# parent taxid
		push(@parent_taxids,$taxid);
		if($taxid == 1){		# root node
			last;
		}
	}
	
	# limit to rank_list
	my @parent_taxids_pruned;
	for my $taxid(@parent_taxids){
		if(!defined($taxnode_ref->{$taxid})){next;}
		my $node	= $taxnode_ref->{$taxid};
		my @sp		= split(/\t/,$node,-1);
		if(defined($parent_rank_set{$sp[2]})){
			push(@parent_taxids_pruned,$taxid);
		}
	}
	
	return @parent_taxids_pruned;
}



