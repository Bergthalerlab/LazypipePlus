#! /usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long qw(GetOptions);
use YAML::Tiny;
use Data::Dumper;
use Env;
#
# DOWNLOAD AND INSTALL SEQ DATABASES FOR LAZYPIPE
#
# credit: Ilya Plyusnin, University of Helsinki, Ilja.Pljusnin@helsinki.fi
#
my $config_file	= "config.yaml";

my $usage= 	"\n".
			"USAGE: $0 --db str [--force -v]\n\n".
			"Download and install local databases according to $config_file preferences\n\n".
			"db str\n".
			"   blastn_vi  : blastn index for virus sequences\n".
			"   centrifuge : centrifuge index\n".
			"   taxonomy   : NCBI taxonomy\n".
			"force         : force overwrite\n".
			"v             : verbal mode\n".
			"\n";

#my %opt	= %{LoadFile( $config_file )};
my $yaml 			= YAML::Tiny->read( $config_file );
my %opt 			= %{$yaml->[0]};

#print STDERR Dumper(%opt);

# expand paths
$opt{'taxonomy'} 		=~ s/\$(\w+)/$ENV{$1}/g;
$opt{'centrifuge_db'}	=~ s/\$(\w+)/$ENV{$1}/g;
$opt{'blastn_vi_db'}	=~ s/\$(\w+)/$ENV{$1}/g;

GetOptions(\%opt,'db=s','force','v') or die $usage;
if( !defined($opt{'db'}) ){
	print STDERR "ERROR: missing input: --db database\n";
	die $usage;
}

if($opt{'db'} =~ m/^taxonomy/gi){
	$opt{'db'} = 'taxonomy';
}
elsif($opt{'db'} =~ m/^blastn/gi){
	$opt{'db'} = 'blastn_vi';
}
elsif($opt{'db'} =~ m/^cent/gi){
	$opt{'db'} = 'centrifuge';
}
else{
	print STDERR "unknown option -db $opt{'db'}\n";
	die $usage;
}



# INSTALL/UPDATE TAXONOMY FILES
if( $opt{'db'} eq 'taxonomy' ){

	my $taxonomy_nodes		= "$opt{'taxonomy'}/nodes.dmp";
	if($opt{'v'}){
		print STDERR "\tinstall  : NCBI taxonomy\n";
		print STDERR "\turl      : $opt{'taxonomy_url'}\n";
		print STDERR "\tlocal    : $opt{'taxonomy'}\n";
	}
	if( !(-e $opt{'taxonomy'}) ){
		system_call("mkdir -p $opt{'taxonomy'}");
	}
	if( !(-e "$taxonomy_nodes")  
		|| ((-M "$taxonomy_nodes") > $opt{'taxonomy_update_time'}) 
		|| $opt{'force'}){
		system_call("wget -q $opt{'taxonomy_url'} -O $opt{'taxonomy'}/taxdump.tar.gz", $opt{'v'});
		system_call("tar -xzf $opt{'taxonomy'}/taxdump.tar.gz -C $opt{'taxonomy'}", $opt{'v'});
	}
	else{
		if($opt{'v'}){
			print STDERR "\ttaxonomy up to date\n";
		}
	}
}

# INSTALL blastn_vi_db
if( $opt{'db'} eq 'blastn_vi' ){
	if($opt{'v'}){
		print STDERR "\n\tinstall  : blastn index\n";
		print STDERR   "\turl      : $opt{'blastn_vi_db_url'}\n";
		print STDERR   "\tlocation : $opt{'blastn_vi_db'}\n";
	}
	if( !(-e "$opt{'blastn_vi_db'}.00.nsq") 
		|| $opt{'force'}){
		my $dir 	= dirname($opt{'blastn_vi_db'});
		my $file 	= basename($opt{'blastn_vi_db_url'});
		system_call("wget -q $opt{'blastn_vi_db_url'} --directory-prefix=$dir", $opt{'v'});
		system_call("tar -xzf $dir/$file -C $dir", $opt{'v'});
	}
	else{
		if($opt{'v'}){
			print STDERR "\tDB already on disk. Use --force to overwrite\n";
		}
	}
}


# INSTALL centrifuge_db
if( $opt{'db'} eq 'centrifuge' ){
	if($opt{'v'}){
		print STDERR "\n\tinstall   : centrifuge index\n";
		print STDERR   "\turl       : $opt{'centrifuge_db_url'}\n";
		print STDERR   "\tlocation  : $opt{'centrifuge_db'}\n";
	}
	if( !(-e "$opt{'centrifuge_db'}.1.cf") 
		|| $opt{'force'}){
		my $dir 	= dirname($opt{'centrifuge_db'});
		my $file 	= basename($opt{'centrifuge_db_url'});
		system_call("wget -q $opt{'centrifuge_db_url'} --directory-prefix=$dir", $opt{'v'});
		system_call("tar -xzf $dir/$file -C $dir", $opt{'v'});
	}
	else{
		if($opt{'v'}){
			print STDERR "\tDB already on disk. Use --force to overwrite\n";
		}
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
