#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# 1) Sorts contig sequences first into main classes
#	1) supkingdomX			: all contig hits belong to superkingdom X, e.g. "viruses"
#	2) supkingdomX-supkingdomY	: all contig hits belong to superkingdom X or Y, e.g. "bacteria-viruses"
#	3) mixture			: contig hits are from a mixture of superkingdoms
#	4) unknown			: no hits for this contig
#       NOTE: contigs that contain sequences from at most two superkingdoms, from which one is from a bacteriophage family are assigned
#             to "bacteriophage" superkingdom. This is convenient for our purposes and we do not claim that this is correct taxonomy.
#	
# 2) In each main class sort into subclasses
#	1) familyX		: all hits belong to the same family
#	2) familyX-familyY	: hits are a mixture of two families
#	3) mixture		: hits are a mixture of tree or more families
#	4) unknown		: all hits belong to an unknown family
# 
# 3) In each subclass sort by genus
#	1) genusX
#	2) genusX-genusY
#	3) mixture
#	4) unknown

# Bacteriophage families according to ICTV classification
#Mc Grath S and van Sinderen D (editors). (2007). Bacteriophage: Genetics and Molecular Biology (1st ed.). Caister Academic Press. ISBN 978-1-904455-14-1. [1].
my @phage_families= ("Myoviridae","Siphoviridae","Podoviridae",
			"Lipothrixviridae","Rudiviridae",
			"Ampullaviridae","Bicaudaviridae","Clavaviridae","Corticoviridae","Cystoviridae","Fuselloviridae","Globuloviridae", "Guttaviridae", "Inoviridae","Leviviridae","Microviridae","Plasmaviridae","Tectiviridae");
my @phage_orders	= ("Caudovirales","Ligamenvirales");

#
# Parameters: references to hashes of the form  'taxonomy_unit_name' -> number_of_occurences_in_the_contig
# 1: \%tax_supking
# 2: \%tax_order	# used for phage classification
# 3: \%tax_family
# 4: \%tax_genus
# 5: \%tax_sciname	# used for phage classification
#
# Returns string:
# SuperkingdomX[-SuperkingdomY]|mixture|unknown tab OrderX[-OrderY]|mixture|unknown tab FamilyX-[FamilyY]|mixture|unknown tab GenusX[-GenusY]|mixture|unknown
#
sub classify_contig{
	my %tax_superking	= %{shift(@_)};
	my %tax_order		= %{shift(@_)};
	my %tax_family		= %{shift(@_)};
	my %tax_genus		= %{shift(@_)};
	my %tax_species		= %{shift(@_)};
	
	my @hash_list = (\%tax_superking,\%tax_order,\%tax_family,\%tax_genus,\%tax_species);
	my $class = '';
	
	for(my $i=0; $i<5; $i++){
		my %h = %{shift(@hash_list)};
		if($i>0){ $class .= "\t";}
		if($h{''}){ delete $h{''};}
		my @k = sort keys %h;
		
		if(scalar(@k)==0){
			$class .= 'unknown';
		}
		elsif(scalar(@k)==1){
			$class .= $k[0];
		}
		elsif( scalar(@k)==2 ){
			$class .= join('+',@k);
		}
		else{
			$class .= 'mixture';
		}
	}
	
	# FILTERING BACTERIOPHAGES BY (1)ORDER (2)FAMILY OR (3)SCPECIES NAME
	my($skingdom,$order,$family,$genus,$species) = split(/\t/,$class,-1);
	
	foreach my $o(@phage_orders){
		if($order =~ m/$o/ig){
			$skingdom = 'bacteriophages';
			return join("\t",$skingdom,$order,$family,$genus,$species);
		}
	}
	
	foreach my $f(@phage_families){
		if($family =~ m/$f/ig){ # this will also match pairs like BACTERIAL_FAMILY-BACTERIOPHAGE_FAMILY
			$skingdom = 'bacteriophages';
			return join("\t",$skingdom,$order,$family,$genus,$species);
		}
	}
	foreach my $s(keys %tax_species){
		if( $species =~ m/\ phage|^phage|\ bacteriophage|^bacteriophage/ig){
			$skingdom = 'bacteriophages';
			return join("\t",$skingdom,$order,$family,$genus,$species);
		}
	}
 	
	return join("\t",$skingdom,$order,$family,$genus,$species);
}


my $usage=  "\nUSAGE: $0 --res|r dir --ann|a contigs.annot.tsv --cont|c contigs.fa [-v]\n".
			"\n".
			"res|r          : result directory, all files will be output here\n".
			"ann|a          : homology search results + taxonomy, *.tsv with named headers\n".
			"                 MUST have these headers:\n".
			"                    contig\n".
			"                    superkingdom\n".
			"                    order\n".
			"                    family\n".
			"                    genus\n".
			"                    species\n". 
			"cont|c         : contigs fasta file\n".
			"-v             : verbal\n\n";
my $res_dir;
my $contigs_annot_file;
my $contigs_file;
my $verb = !1;
GetOptions(
		'res|r=s' 	=> \$res_dir,
    	'ann|a=s' 	=> \$contigs_annot_file,
		'cont|c=s'  => \$contigs_file,
		'v'		=> \$verb) or die $usage;
if(!$res_dir || !$contigs_annot_file || !$contigs_file){
	die $usage;
}


# READING BLAST AND TAXONOMY TO ARRAY
# each contig is classified to main class and a number of subclasses based on taxonomy linked to blast hits
# classification is saved in "seqid" -> "main\tsubclass" hash
my %class_hash;

if($verb){ print STDERR "\treading $contigs_annot_file\n";}
open(IN,"<$contigs_annot_file") or die "Can\'t open $contigs_annot_file: $!\n";

my %header; # header_label > column
my @required= ("contig","superkingdom","order","family","genus","species");
my $l=<IN>; chomp($l);
my @sp		= split(/\t/,$l,-1);
for(my $i=0; $i<scalar(@sp); $i++){
	$header{lc($sp[$i])} = $i;
}
#if( defined($header{'node'}) ){ 
#	$header{'contig'} = $header{'node'}; }
#if( defined($header{'contig_id'}) ){
#	$header{'contig'} = $header{'contig_id'}; }
foreach my $h(@required){
	if( !defined($header{$h})){  die "# ERROR: missing header \"$h\"\n";}
}


$l=<IN>; chomp($l);
@sp = split(/\t/,$l,-1);
my $contig_id	= $sp[$header{'contig'}];
my $new_id;
my (%tax_supking,%tax_order,%tax_family,%tax_genus,%tax_species);
$tax_supking{ 	$sp[$header{'superkingdom'}] }	= 1;
$tax_order{ 	$sp[$header{'order'}] }		= 1;
$tax_family{ 	$sp[$header{'family'}] } 	= 1;
$tax_genus{ 	$sp[$header{'genus'}] }		= 1;
$tax_species{ 	$sp[$header{'species'}] }	= 1;

while($l=<IN>){
	chomp($l);
	@sp 	= split(/\t/,$l,-1);
	$new_id = $sp[$header{'contig'}];
	if( $contig_id ne $new_id ){
		my $class = classify_contig(\%tax_supking,\%tax_order,\%tax_family,\%tax_genus,\%tax_species);
		$class_hash{$contig_id} = $class;
		%tax_supking 	= ();
		%tax_order	= ();
		%tax_family 	= ();
		%tax_genus	= ();
		%tax_species	= ();
	}

	$contig_id = $new_id;
	$tax_supking{$sp[$header{'superkingdom'}]} ? ($tax_supking{$sp[$header{'superkingdom'}]}++) : ($tax_supking{$sp[$header{'superkingdom'}]}=1);
	$tax_order{$sp[$header{'order'}]} ? ($tax_order{$sp[$header{'order'}]}++) : ($tax_order{$sp[$header{'order'}]}=1);
	$tax_family{$sp[$header{'family'}]} ? ($tax_family{$sp[$header{'family'}]}++) : ($tax_family{$sp[$header{'family'}]}=1);
	$tax_genus{$sp[$header{'genus'}]}  ?  ($tax_genus{$sp[$header{'genus'}]}++)  : ($tax_genus{$sp[$header{'genus'}]}=1);
	$tax_species{$sp[$header{'species'}]} ? ($tax_species{$sp[$header{'species'}]}++) : ($tax_species{$sp[$header{'species'}]}=1);
}
close(IN);
# LAST RECORD
if(scalar(keys %tax_supking) > 0){
	my $class = classify_contig(\%tax_supking,\%tax_order,\%tax_family,\%tax_genus,\%tax_species);
	$class_hash{$contig_id} = $class;
}


# READ CONTIG SEQS AND SORT THEM ACCORDING TO TAXONOMY

if($verb){ print STDERR "\tsorting contigs by taxonomy\n";}
system("rm -fR $res_dir");
mkdir("$res_dir"); # this scripts uses file appending, all files are replaced on each call

my %contigs = 	%{readfasta($contigs_file)};

foreach my $seqid(sort keys %contigs){
	
	my $contig_id = $seqid;
	if($contig_id =~ m/^\w+=([A-Za-z0-9\.]+)[_\s]?/g){
		$contig_id = $1;
	}
	else{
		print STDERR "ERROR: unknown contig_id format: contig=\'$contig_id\' in $contigs_file\n";
		exit(1);
	}	
	
	my $class = ($class_hash{$contig_id}) ? ($class_hash{$contig_id}) : "unknown\tunknown\tunknown\tunknown\tunknown";
	my($skingdom,$order,$family,$genus,$species) = split(/\t/,$class,-1);
	
	my $class_dir		= "$res_dir/$skingdom";
	my $subclass_dir 	= "$res_dir/$skingdom/$family";
	my $subclass_file	= "$res_dir/$skingdom/$family".'.fa';
	my $subsubclass_file	= "$res_dir/$skingdom/$family/$genus".'.fa';
	
	if(!(-e "$res_dir/$skingdom")){
		mkdir("$res_dir/$skingdom");
	}
	
	# COLLECTIVE FILE FOR FAMILY SEQS
	my $mod = (-e "$subclass_file")? ">>": ">";
	open(OUT,"$mod$subclass_file") or die "Can\'t open $subclass_file: $!\n";
	print OUT ">$seqid\n";
	print OUT "$contigs{$seqid}\n";
	close(OUT);
	
	
	# separate subsubclass files
	if(!(-e "$subclass_dir")){
		mkdir("$subclass_dir");
	}
	$mod	= (-e "$subsubclass_file")? ">>": ">";
	open(OUT,"$mod$subsubclass_file") or die "Can\'t open $subsubclass_file: $!\n";
	print OUT ">$seqid\n";
	print OUT "$contigs{$seqid}\n";
	close(OUT);
}
	
sub readfasta{
  	my $file		= shift(@_);
	my %sequence;
	my $header;
	my $temp_seq;
	
	#suppose fasta files contains multiple sequences;
	 
	open (IN, "<$file") or die "couldn't open the file $file $!";
	
	while (<IN>){	
		chop;
		next if /^\s*$/; #skip empty line 
		if ($_ =~ s/^>//)  #when see head line
		{	
			$header= $_;
			if ($sequence{$header}){print colored("#CAUTION: SAME FASTA HAS BEEN READ MULTIPLE TIMES.\n#CAUTION: PLEASE CHECK FASTA SEQUENCE:$header\n","red")};
			if ($temp_seq) {$temp_seq=""} # If there is alreay sequence in temp_seq, empty the sequence file
			
		}
		else # when see the sequence line 
		{
		   s/\s+//g;
		   $temp_seq .= $_;
		   $sequence{$header}=$temp_seq; #update the contents
		}
	
	}
	
	return \%sequence;
}
