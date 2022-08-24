#!/usr/bin/perl
use strict;
use warnings;

# 
# Converts MetaGeneAnnotator predictions to GTF2.2 format
#
my $usage= "\nUSAGE: $0 ORF.mga 1> ORF.gtf\n".
		   "\n".
		   "ORF.mga    \t: ORF predictions by MetaGeneAnnotator (Noguchi 2008)\n".
		   "\n".
		   "converts MetaGeneAnnotator ORF predictions to GTF2.2 format \n\n";

if(scalar(@ARGV)<1) { 
	print "$usage";
	exit(1);
}
my $file_in	= shift(@ARGV);
my $set_frame_zero	= 1;	# shift frame to 0-frame by adjusting start/end positions

# READ IN MGA ORF PREDICTIONS
# MGA format (ref: http://metagene.nig.ac.jp/metagene/download_mga.html):
# # [sequence name]
# # gc = [gc%], rbs = [rbs%]
# # self: [(b)acteria/(a)rchaea/(p)hage/unused(-)]
# [0:gene ID] [1:start pos.] [2:end pos.] [3:strand] [4:frame] [5:complete/partial] [6:gene score] [7:used model] [8:rbs start] [9:rbs end] [10:rbs score]

my $contig_id;
my @gene_pred_list;
my %gene_pred_hash;

open(IN,"<$file_in") or die "Can\'t open $file_in: $!\n";
while(my $l=<IN>){
		chomp($l);
        if( $l =~ m/^[\!#]/){
        	if(scalar(@gene_pred_list) > 0){# genes predictions read for prev contig
	    		my @copy = @gene_pred_list;
				$gene_pred_hash{$contig_id} = \@copy;
	        	@gene_pred_list=();
	    	}
	  
			$l =~ s/^[#\s]+|\s+$//g; # remove leading/trailing spaces
			my @sp= split(/\s+/,$l);
			$contig_id	= $sp[0];
			# DEBUG
			#print STDERR "contig_id:\"$contig_id\"\n";exit(1);	
			# reading two more comment lines
			$l=<IN>;
			$l=<IN>;
			next;
        }
        
		if( $l =~ m/^gene/i ){
			push(@gene_pred_list,$l);
		}
}
if(scalar(@gene_pred_list)>0){
	my @copy = @gene_pred_list;
	$gene_pred_hash{$contig_id} = \@copy;
}

# WRITE GTF2.2: <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments] 

foreach my $seq( sort keys %gene_pred_hash){
	
	my @genepred_list =  @{$gene_pred_hash{$seq}};
	
	foreach my $genepred( @genepred_list){
	
		my ($geneid,$start,$end,$strand,$frame,$complete,$score,$model,$rbs_start,$rbs_end,$rbs_score)
		  	= split(/\t/,$genepred,-1);
		if($set_frame_zero && ($strand eq '+')){
			$start  = $start + $frame;
			$frame	= 0;
		}
		if($set_frame_zero && ($strand eq '-')){
			$end	= $end - $frame;
			$frame	= 0;
		}
		
		print "$seq\tbwa\tORF\t$start\t$end\t$score\t$strand\t$frame\tgene_id \"\"; transcript_id \"\"; complete $complete; model $model;\n";
	}
}




