#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# Converts MetaPhlan2 output to taxonomy profile
#
# input: nodes.dmp from NCBI taxonomy
#	 names.dmp from NCBI taxonomy containing taxnames output by metaphlan2
#	 taxonomy profile output by default metaphlan2.py run
# 
#
# output: taxonomy profile in CAMI format (see https://github.com/bioboxes/rfc/blob/master/data-format/profiling.mkd)
#

my $usage=      "USAGE: $0 nodes.dmp names.dmp metaphlan2.profile\n".
		"\n".
		"nodes.dmp \t\t: node file from NCBI taxonomy\n".
                "names.dmp \t\t: name file from NCBI taxonomy\n".
		"metaphlan2.profile\t: metaphlan2 output with default options\n".
		"nodes.dmp & names.dmp can be downloaded from:\n".
		"ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz\n".
		"\n";


if(scalar @ARGV<3) { die $usage; }

my $node_dmp_file	= shift(@ARGV);
my $name_dmp_file	= shift(@ARGV);
my $metaphlan_file	= shift(@ARGV);
my $ignore_strains	= 1;

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
my %taxid_sciname;	# taxid > scientific taxname
my %taxid_synonym;	# taxid > synonym name
my %taxid_comname;	# taxid > common name
my %taxid_other;	# taxid > other name
my %taxid_allnames;	# taxid > hash pointer: name > 1
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
	
	if( lc($sp[3]) eq 'scientific name'){
		$taxid_sciname{$sp[0]}	= $sp[1];
	}
	elsif( lc($sp[3]) eq 'synonym'){
		$taxid_synonym{$sp[0]}	= $sp[1];
	}
	elsif( lc($sp[3]) eq 'common name'){
		$taxid_comname{$sp[0]}	= $sp[1];
	}
	else{
		$taxid_other{$sp[0]}	= $sp[1];
	}
	
	if(!defined($taxid_allnames{$sp[0]})){ 
		my %tmp;
		$taxid_allnames{$sp[0]} = \%tmp;}
		
	$taxid_allnames{$sp[0]}->{$sp[1]} = 1;
}

my %taxname_id;
for my $id(keys %taxid_allnames){
	for my $name(keys %{$taxid_allnames{$id}}){
		$taxname_id{$name} = $id;
	}
}

my %rankchar_rank = 
	('k' => 'superkingdom',
	 'p' => 'phylum',
	 'c' => 'class',
	 'o' => 'order',
	 'f' => 'family',
	 'g' => 'genus',
	 's' => 'species',
	 't' => 'strain');
		

# EXAMPLE USAGE OF %taxnodes & %taxid_names  hashes:
#my $id= 1656;
#my @sp= split(/\t/, $taxnodes{$id});
#print "hashnode: ",$taxnodes{$id},"\n";
#print "# node $id:\n";
#print "id         :",$sp[0],"\n";
#print "parent id  :",$sp[1],"\n";
#print "rank       :",$sp[2],"\n";
#print "name       :",$taxnames{$sp[0]},"\n";

print STDERR "# reading $metaphlan_file\n";
open(IN,"<$metaphlan_file") or die "Can\'t open $metaphlan_file: $!\n";
my $l = <IN>;chomp($l);
my $sample_id		= $metaphlan_file;
if($l =~ s/^#//){
	my @sp = split(/\t/,$l,-1);
	if(defined($sp[1]) && $sp[1] ne ""){
		$sample_id = $sp[1];
	}
}
close(IN); 
my $taxnames_linked 	= 0;
my $taxnames_noname	= 0;
my $taxnames_unclass	= 0;
my $taxnames_unknown	= 0;
my %taxnames_unknown_h;	#line\ttaxname

print '@',"Version:0.9.1\n";
print '@',"SampleID:$sample_id\n";
print '@',"Ranks:superkingdom|phylum|class|order|family|genus|species\n";
print '@',"TaxonomyID:NCBI taxonomy\n";
print '@@',"TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n";


open(IN,"<$metaphlan_file") or die "Can\'t open $metaphlan_file: $!\n";
$ln=0;
while(my $l=<IN>){
	$ln++;
	chomp($l);
	if($l =~ m/^#/){ next;}
	my ($taxpath_str,$abund)= split(/\t/,$l,-1);
	my @taxpathsn	= split(/\|/,$taxpath_str);
	

	#parsing and removing rank char from taxpathsn
	#translate "_" to " "
	#translate "sp" to "sp."
	#translate "Clostridiales Family XI" to "Clostridiales Family XI."
	my @taxpath_ranks;
	for(my $i=0; $i<scalar(@taxpathsn);$i++){
		$taxpathsn[$i] =~ s/^([kpocfgst])__//;
		my $rank= $rankchar_rank{$1};
		push(@taxpath_ranks,$rank);
		$taxpathsn[$i] =~ s/([A-Za-z])_/$1 /g;
		$taxpathsn[$i] =~ s/ sp / sp\. /g;
		$taxpathsn[$i] =~ s/Clostridiales Family XI/Clostridiales Family XI\./g;
		$taxpathsn[$i] =~ s/ like/-like/g;
		$taxpathsn[$i] =~ s/Tannerella sp\. 6_1_58FAA CT1/Tannerella sp\. 6_1_58FAA_CT1/g;
	}
	
	# deleting strain-level info, when ignore_strains==true
	if($ignore_strains){
		if($taxpath_ranks[$#taxpath_ranks] eq 'strain'){
			pop(@taxpath_ranks);
			pop(@taxpathsn);
		}
	}
	
	#forming taxpath
	my @taxpath;
	for(my $i=0; $i<scalar(@taxpathsn); $i++){
		if(!defined( $taxname_id{$taxpathsn[$i]})){
			if( ($taxpathsn[$i] =~ m/noname/g) ){
				$taxnames_noname++;
			}
			elsif( ($taxpathsn[$i] =~ m/unclassified$/)){
				$taxnames_unclass++;
			}
			else{
				$taxnames_unknown++;
				$taxnames_unknown_h{$taxpathsn[$i]} = 1;
				print STDERR "unknown taxname\tline $ln:\t$taxpathsn[$i]\n";
				$taxpathsn[$i] = "";
			}
			push(@taxpath,"");
		}
		else{
			push(@taxpath,$taxname_id{$taxpathsn[$i]});
		}
	}
	
	if($taxpath[$#taxpath] ne ""){
		print join("\t",$taxpath[$#taxpath], $taxpath_ranks[$#taxpath], join('|',@taxpath), join('|',@taxpathsn), $abund),"\n";
		$taxnames_linked++;
	}
}
$ln--;

print STDERR "# profile lines with all taxnames linked\/lines\t$taxnames_linked\/$ln\n";
print STDERR "# taxnames with \"noname\"-labels\t$taxnames_noname\n";
print STDERR "# taxnames with \"unclassified\"-labels\t$taxnames_unclass\n";
print STDERR "# taxnames unknown\t$taxnames_unknown\n";
#print STDERR "# taxnames unknown:\n";
#for my $t(sort keys %taxnames_unknown_h){
#	print STDERR "$t\n";
#}


