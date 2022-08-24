#
# Interface to NGSlib.R get_abund_tables
# 
# prints abundance table to excel
#
source("R/NGSlib.R");
if( !(library("openxlsx", lib.loc="~/Rlib",logical.return=T)) ){
	library("openxlsx");
}
 
usage	= paste("USAGE: Rscript print_abund_table abund_table.txt annot_table.txt abund.xlsx [tail host_taxid]",sep="\n");
		
args 	= commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
	stop(paste("Missing input files\n",usage,sep="\n"), call.=FALSE);
}
if(length(args)<4){
	args[4]<- 1;
}
if(length(args)<5){
	args[5]<- 0;
}


sprintf("# printing abundance table:");
sprintf("#   get_abund_tables(abund_table=%s, annot_table=%s, excel_file=%s, tail=%d, host_taxid=%s)",args[1],args[2],args[3],as.numeric(args[4]),args[5]);

abund<- get_abund_tables(abund_table=args[1], annot_table=args[2], excel_file=args[3], tail=as.numeric(args[4]), host_taxid=args[5]);
