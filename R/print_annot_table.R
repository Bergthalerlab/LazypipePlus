#
# Interface to NGSlib.R get_abund_tables and get_annot_tables
# 
# prints annotation table to excel
# when append_prefix specified, appends annotation table to the existing excel file
#
source("R/NGSlib.R");
if( !(library("openxlsx", lib.loc="~/Rlib",logical.return=T)) ){
	library("openxlsx");
}

usage	= paste("USAGE: Rscript print_annot_table annot_table.txt annot.xlsx [append_prefix]",sep="\n");
		
args 	= commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
	stop(paste("Missing input files\n",usage,sep="\n"), call.=FALSE);
}

sprintf("# printing annotation table:");
sprintf("#    %s > %s",args[1], args[2]);

if( length(args)<3 ){
	annot<- get_annot_tables(annot_file=args[1], excel_file=args[2]);
}else{
	annot<- get_annot_tables(annot_file=args[1], excel_file=args[2], append_prefix=args[3]);
}


