#! /usr/bin/perl
use strict;
use warnings;
#use YAML::XS 'LoadFile';
use YAML::Tiny;
use Getopt::Long qw(GetOptions);

#
# UPDATE TAXONOMY DATABASE FOR LAZYPIPE
#
# credit: Ilya Plyusnin, University of Helsinki, Ilja.Pljusnin@helsinki.fi
#
my $config_file	= "config.yaml";

my $usage= 	"\n".
			"USAGE: $0 [--force -v]\n\n".
			"       Taxonomy DB will be updated based on settings in $config_file\n\n";

#my %opt	= %{LoadFile( $config_file )};
my $yaml 			= YAML::Tiny->read( $config_file );
my %opt 			= %{$yaml->[0]};
my $taxonomy_nodes	= "$opt{'taxonomy'}/nodes.dmp";
$taxonomy_nodes 	=~ s/\$(\w+)/$ENV{$1}/g;

GetOptions(\%opt,'force','v') or die $usage;

# UPDATING TAXONOMY FILES: this will also load taxonomy on the very first usage
if( !(-e $opt{'taxonomy'}) ){
	system_call("mkdir -p $opt{'taxonomy'}");
}
if( !(-e "$taxonomy_nodes")  
	|| ((-M "$taxonomy_nodes") > $opt{'taxonomy_update_time'}) 
	|| $opt{'force'}){
	if($opt{'v'}){
		print STDERR "\tupdating taxonomy db\n";
	}
	system_call("wget -q $opt{'taxonomy_url'} -O $opt{'taxonomy'}/taxdump.tar.gz", $opt{'v'});
	system_call("tar -xzf $opt{'taxonomy'}/taxdump.tar.gz -C $opt{'taxonomy'}", $opt{'v'});
}
else{
	if($opt{'v'}){
	print STDERR "\ttaxonomy db up to date\n";
	}
}


# system_call($system_call, $verbal=false)
#
# $system_call	String with the system call to excecute
# $verbal	Set true to print the call prior to excecuting the call. Optional, default is false.
# 
sub system_call{
	my $call	= shift;
	my $verbal	= (scalar(@_)>0) ? shift(@_): 0;
	if($verbal){
		print STDERR "\t$call\n";
	}
	my @args= ("bash","-c",$call);
	system(@args) == 0 or die $!;
}
