#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $USAGE = "USAGE: $0 [-h|headers --rgabund int --rgtaxid int --t|cont_score_tail int -w|weights taxacount|bitscore|bitscore2 ] contigid+taxid+score idx_file > readn+contign+taxid\n\n";   

my $headers 		= !1;
my $rgabund 		= !1;
my $rgtaxid 		= !1;
my $weight_model  	= 'taxacount'; # options: taxacount|bitscore|bitscore2
my $cont_score_tail	= 0;

GetOptions(	'headers|h'	=> \$headers,
		'rgabund=i'	=> \$rgabund,
		'rgtaxid=i'	=> \$rgtaxid,
		'weights|w=s'	=> \$weight_model,
		'cont_score_tail|t=i'	=> \$cont_score_tail);
		
$weight_model 	= lc($weight_model);
if($weight_model =~ /taxa|taxon/gi){ $weight_model = 'taxacount';}
if( !(($weight_model =~ 'taxacount') || ($weight_model eq 'bitscore') || ($weight_model eq 'bitscore2')) ){
	print STDERR "ERROR: invalid argument: --weights\n\n";
	print STDERR $USAGE; exit(1);
}
if($cont_score_tail > 10 || $cont_score_tail < 0){
	print STDERR "ERROR: invalid argument: --cont_score_tail $cont_score_tail. Use values in [0,10]\n\n";
	print STDERR $USAGE; exit(1);
}

if( scalar(@ARGV)<2){
	print STDERR "ERROR: missing input files\n\n";
	print STDERR $USAGE; exit(1);
}
my $file1 	= shift(@ARGV);
my $idx_file 	= shift(@ARGV);


my %tax2cont;		# $taxid->$contid->0/1		: taxid mapped to a unique set of contig-ids (%hash->%hash->0/1)
my %cont2tax;		# $contid->$taxid->0/1		: contig-id mapped to a unique set of tax-ids (%hash->%hash->0/1)
my %cont2score;		# $contid->$score		: sum of bitscores for a given contig (parsed from db_hits)
my %cont2tax2score; 	# $contid->$taxid->$score	: sum of bitscores for a given contig and taxid (%hash->%hash->double)
my %tax2score;		# $taxid->$score		: sum of bitscores for a given taxon

print STDERR "# reading $file1\n";
open(IN, "<$file1") or die "failed to open $file1\n";
my $ln  = 0;
while(my $l=<IN>) {
	$ln++;
	chomp($l);
	if($ln==1 && $headers){
		next;
	}
	my ($cont,$tax,$score)	= (0,0,0);
	($cont,$tax,$score) 	= split(/\t/,$l,-1);
	
	if( defined($tax2cont{$tax}) ){
		$tax2cont{$tax}->{$cont} = 1;
	}
	else{
		my %tmp =($cont => 1);
		$tax2cont{$tax} = \%tmp;
	}
	
	if( defined($cont2tax{$cont}) ){
		$cont2tax{$cont}->{$tax} = 1;
	}
	else{
		my %tmp = ($tax => 1);
		$cont2tax{$cont} = \%tmp;
	}
	
	if( defined($cont2score{$cont}) ){
		$cont2score{$cont} += $score;
	}
	else{
		$cont2score{$cont} = $score;
	}
	
	if( defined($cont2tax2score{$cont})  ){
		if( defined($cont2tax2score{$cont}->{$tax}) ){
			$cont2tax2score{$cont}->{$tax} += $score;
		}
		else{
			$cont2tax2score{$cont}->{$tax} = $score;
		}
	}
	else{
		my %tmp = ($tax => $score);
		$cont2tax2score{$cont} = \%tmp;
	}
	
	if( defined($tax2score{$tax}) ){
		$tax2score{$tax} += $score;
	}
	else{
		$tax2score{$tax} = $score;
	}
}
close(IN);


# DELETING TAXA THAT ARE IN THE "TAIL" OF SCORE DIST FOR EACH CONTIG
if($cont_score_tail > 0){
  foreach my $cont(keys %cont2tax2score ){

	my $csum 	= 0;
	my %tmp_hash 	= %{$cont2tax2score{$cont}};
	
	foreach my $tax(sort { $tmp_hash{$a} <=> $tmp_hash{$b} } keys %tmp_hash ){	
		$csum 	+= $tmp_hash{$tax};
		if($csum < $cont2score{$cont}*($cont_score_tail/100.0)){
			#print STDERR "deleting\t$cont\t$tax\n";
			delete($cont2tax{$cont}->{$tax});
			delete($tax2cont{$tax}->{$cont});
			delete($cont2tax2score{$cont}->{$tax});
		}
		else{
			last;
		}
	}	
  }
}


print STDERR "# reading $idx_file\n";
my %cont2rn;
open(IN,"<$idx_file") or die "Can\'t open $idx_file: $!\n";
while(my $l=<IN>){
	chomp($l);
	my @sp = split(/\t/,$l,-1);
	if($sp[0] eq "*"){
		next; 	# unaligned reads
	}
	elsif( $sp[0] =~ m/^contig=([A-Za-z0-9\.]+)[_\s]?/ ){
		$cont2rn{$1} = $sp[2];
	}
	elsif( $sp[0] =~ m/^([A-Za-z0-9\.]+)[_\s]?/ ){
		$cont2rn{$1} = $sp[2];
	}
	else{
		$cont2rn{$sp[0]} = $sp[2];
	}	
}
close(IN);
#foreach my $cont(keys %cont2rn){print STDERR "$cont\t$cont2rn{$cont}\n";}

my %tax2readn;
my %tax2contn;
foreach my $tax(sort keys %tax2cont){
	my @contids		= keys %{$tax2cont{$tax}};
	my $readn		= 0;
	
	foreach my $cont(@contids){
		my $contig_taxa_score 	= 0;
		foreach my $taxon(keys %{$cont2tax{$cont}}){	
			$contig_taxa_score += $tax2score{$taxon};
		}
	
		if(defined($cont2rn{$cont})){
			if($weight_model eq 'taxacount'){
				$readn += $cont2rn{$cont} / scalar(keys %{$cont2tax{$cont}} ); # assign to this taxon readn assigned to linked contig, divided by number of taxa linked to that contig
			}
			if($weight_model eq 'bitscore'){
				$readn += $cont2rn{$cont} * ($cont2tax2score{$cont}->{$tax})/($cont2score{$cont});
			}
			if($weight_model eq 'bitscore2'){
				$readn += $cont2rn{$cont} * ($tax2score{$tax}/$contig_taxa_score);
			}
		};
	}
	$tax2readn{$tax} = $readn;
	$tax2contn{$tax} = scalar(@contids);
	#print STDERR "TAXON=$tax\tREADN=$readn\n";
}

if($rgabund && $rgtaxid){
	if(!defined($tax2readn{$rgtaxid})){
		$tax2readn{$rgtaxid} = 0;
		$tax2contn{$rgtaxid} = 0;
	}
	$tax2readn{$rgtaxid} += $rgabund;
}

print "readn\tcontign\ttaxid\n";
foreach my $tax(sort {$tax2readn{$b}<=>$tax2readn{$a}} keys %tax2readn){
	if($tax2readn{$tax} == 0){
		next;
	}
	print "",join("\t",$tax2readn{$tax},$tax2contn{$tax},$tax),"\n";
}




