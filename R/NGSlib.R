# 
# FUNCTIONS FOR NGS PIPE
#

#
# Sorts files in a list using one of defined rules:
# 'numsuffix'    by numeric suffix, e.g. 'file_B_1' > 'file_A_2' > 'file_A3'
sort.files<- function(filelist, rule='numsuffix'){
    i<- order(as.numeric(gsub("(.*[A-z_-]+)([0-9]+)([A-z]?)$", "\\2", filelist)));
    return(filelist[i]);
}

# Converts filenames to short labels
# by default uses string prefix + the first occuring number
file2label<- function(filelist){
}

# converts data frame to flat format 
# (columnwise catenation with addition of rowname and colname columns)
convert_to_flat<- function(data){
    m<- nrow(data);
    n<- ncol(data);
    data_flat<- as.data.frame(matrix(0, nrow= m*n, ncol=3));
    
    colnames(data_flat)<- c("rowname","colname","value");
    data_rownames<- rownames(data);
    data_colnames<- colnames(data);
    
    for(i in 1:m){
        for(j in 1:n){
            data_flat[(j-1)*m+i,1]<- data_rownames[i];
            data_flat[(j-1)*m+i,2]<- data_colnames[j];
            data_flat[(j-1)*m+i,3]<- data[i,j];
        }
    }
    data_flat$rowname<- factor(data_flat$rowname, levels = unique(data_flat$rowname));
    data_flat$colname<- factor(data_flat$colname, levels = unique(data_flat$colname))
    return(data_flat)
}

# QUALITY CONTROL: READ LENGTH HISTOGRAM
#
seqlen_hist<- function(file_in, file_out = "none", name = "none", width=5.2, res=300, pointsize=9){

    # DEFAULTS
    if(file_out == "none"){ file_out<- gsub(".(\\w*)$",".hist1.pdf",file_in,perl=T);}
    if(name == "none"){ name = gsub(".(\\w*)$","",file_in,perl=T);}
    
    if(file_out!="none"){
        #jpeg(filename = file_out, width=width*res, height=width*res, units="px", res=300, pointsize=pointsize);
        pdf(file = file_out, width=width*res, height=width*res)
    }     
    # HISTOGRAM FROM LENGTH VECTOR FILE
    if(length(grep(".(len|length)$",file_in,perl=T,ignore.case=T)) > 0){
	    readlen<- read.csv(file=file_in, sep="", header=FALSE)
	    readlen<- readlen[,1]
        hist(readlen, breaks=seq(0,max(readlen)+10,by=10), main=name, xlab="seq length");
    	#par(new=T);
	    #plot(density(readlen), main="", ylab="",xlab="",col="red", yaxt="n")
	    axis(4);
	
        # mean+median+quantiles
        q<- quantile(readlen,probs = seq(0,1,0.25))
	    s<- sprintf("mean\t: %1.1f\nmedian\t: %1.1f\nQ0\t: %1.1f\nQ25\t: %1.1f\nQ50\t: %1.1f\nQ75\t: %1.1f\nQ100\t: %1.1f",mean(readlen),median(readlen),q[1],q[2],q[3],q[4],q[5]);
	    mtext(s,side=3,adj = 0.1,padj = 1.1);
    }
    if(file_out!="none"){
        cat("# histogram printed to ",file_out,"\n")
        dev.off();
    }    
}


# Plot read numbers from raw to reads in genes
# filelist      list of assembly.stats.txt files
plot_readcount<- function(filelist, file_out="none",title="",legend=F,labels=basename(dirname(filelist)),
                    plotpc = F, width=5.2, res=300, pointsize=9){
    
    if(require(colorspace, quietly=T) == TRUE){
    	pal<- rainbow_hcl(n=length(filelist));
    }
    else{
    	pal<- rep("#000000",length(filelist));
    }

    
    data<- read.table(file=filelist[1],header=F,sep="\t",row.names=1);
    if(length(filelist) >1){
        for(file in filelist[2:length(filelist)]){
            column<- read.table(file = file,header = F,sep="\t",row.names=1);
            data<- cbind(data,column);
        }
    }
    data<- t(data);
    rownames(data)<- 1:nrow(data);
    data<- as.data.frame(data);
    
    # converting survival rates to pcs
    if(plotpc){
        data$reads.flt    <- (data$reads.flt)/data$reads;
        data$reads.hgflt  <- data$reads.hgflt/data$reads;
        data$reads.contigs <- data$reads.contigs/data$reads;
        data$reads.orfs  <- data$reads.orfs/data$reads;
        data$reads        <- data$reads/data$reads;

        yticks<- seq(0, 1,by = 0.1);
        ylabels<- sprintf("%.0f",yticks*100);
        ylim<- c(0,1);
        ylab<- "Number of reads (%)";
    }
    else{
        #ysteP <- ceiling(max(data$reads)/10);
        yticks<- seq(0, max(data$reads)*1.05,length.out=11);
        ylabels<- sprintf("%iK",yticks/1000);
        ylim<- c(0, max(data$reads)*1.05);
        ylab<- "Number of reads";
    }

    par(mar= c(5,6,4,9)+0.1);
    par(mgp= c(4,1,0), xpd=T);

    if(file_out!="none"){
        #jpeg(filename = file_out, width=width*res, height=width*res, units="px", res=300, pointsize=pointsize);
        pdf(file = file_out, width=width*res, height=width*res);
    }    
    for(i in 1:nrow(data)){
        if(i>1){ par(new=T)};

        plot(c(1:5),data[i,c("reads","reads.flt","reads.hgflt","reads.contigs","reads.orfs")],
            type='b',lty=1,lwd=2,col=pal[i],pch=i,ylab =ylab, xlab="", yaxt="n",xaxt="n",ylim=ylim);
    }
    axis(1, at=c(1:5),labels=c("reads","reads.flt","reads.hgflt","reads.contigs","reads.orfs"), tick=T,las=1)
    axis(2, at=yticks,labels=ylabels, tick=T,las=2)
    if(title != ""){ title(main=title); }
    #if(subtitle != ""){ title(sub=subtitle);}
    if(legend){
        legend("right",legend=labels,lty=1,lwd=2,pch=1:nrow(data),col=pal, bty="n",cex=1,inset=-0.4);
    }

    if(file_out!="none"){
        dev.off();
    }
    par(mar= c(5,4,4,2)+0.1, xpd=F);
}

# input:
# stats     assembly stats, including stats$reads
# timer     timer stats, including timer$cpu
plot_timer1<- function(stats,timer){
    ystep = 100000;
    yticks<- seq(0, ceiling(max(stats$reads)/ystep)*ystep, by=ystep);
    ylabels<- sprintf("%iK",yticks/1000);
    ylim<- c(0, ceiling(max(stats$reads)/ystep)*ystep);    
    
    par(mar= c(5,6,4,2)+0.1);
    
    mod1<- lm(stats$reads ~ timer$cpu)
    plot(timer$cpu, stats$reads, type='p',lty=1,lwd=2,col="blue",pch=15,
        xlab="CPU(s)",ylab="Number of reads",yaxt="n",ylim=ylim)
    axis(2, at=yticks,labels=ylabels, tick=T,las=2)
    title(main="CPU performance")
    abline(mod1, lty=1,lwd=2,col="black")
    
    s<- sprintf("reads/s   : %1.1f\n",mod1$coefficients[[2]]);
	mtext(s,side=3,adj = 0.1,padj = 1.5,cex=1.2);
    
    #pre <- predict(mod1) # plot distances between points and the regression line
    #segments(stats$reads,timer$cpu,stats$reads,pre,col="red")
    
    #mod2<- lm(stats$reads ~ timer$real);
    par(mar= c(5,4,4,2)+0.1, xpd=F);
}

# boxplot for time perf across pipeline parts
plot_timer2<- function(timer2,lwd=1){
    
    tmp<- convert_to_flat(timer2);
    colnames(tmp)<- c("sample","pipeline_step","cpu");
    plot(tmp$pipeline_step,tmp$cpu, main="Stepwise CPU",lwd=lwd, xlab="pipeline step", ylab="CPU(%)");
    labels<- c("1:histogramm","2:Trimmomatic filter","3:refgen filter","4:assemble","5:MetaGene",
                "6:SANSparallel","8:realign reads","9:sort by taxa","10:realign within taxa",
                "11:stats","12:collect results","13:cleanup");
    legend("topright",legend=labels,col="black",cex=1,inset=0.05);
}

# EVALUATES PERFORMANCE ON A GIVEN BENCHMARK
# pred_file     : taxa predictions as data.table with fields: <readn,species,species_id,genus,genus_id,family,family_id,superkingdom, superkingdom_id>
# true_file     : true abundance as data.table (same format)
# excel_file    : print resutls to this file
# filter_phages : filter phages by family/order names
# host_taxid    : reference genome taxid to exclude from predictions and read dist extimation
# tail          : the least abundant taxa that form up to tail% of read distribution will be ignored
#
bmark_stats<- function(pred_file,true_file,excel_file="none",filter_phages=T,host_taxid=9606, tail=1){
    library(reshape);
    library(openxlsx);
    res<- list();
    
    # READ DATA
    pred<- read.table(pred_file,sep="\t",header=T, comment.char="",stringsAsFactors = F);
    true<- read.table(true_file, sep="\t",header=T, comment.char = "",stringsAsFactors = F);
    
    # EXCLUDE HOST TAXID
    ind<- which(pred$taxid == host_taxid);
    if(length(ind)>0){
        pred<- pred[-ind,];
    }
    
        
    # DETECT PHAGES
    if(filter_phages){
    phage_families= c("Myoviridae","Siphoviridae","Podoviridae","Lipothrixviridae","Rudiviridae",
                    "Ampullaviridae","Bicaudaviridae","Clavaviridae","Corticoviridae","Cystoviridae",
                    "Fuselloviridae","Globuloviridae", "Guttaviridae", "Inoviridae","Leviviridae",
                    "Microviridae","Plasmaviridae","Tectiviridae");
    phage_orders = c("Caudovirales","Ligamenvirales");
    ind<- union(which(pred$family %in% phage_families),which(pred$order %in% phage_orders));
    ind<- union(ind,grep("\\<phage\\>|\\<bacteriophage\\>",pred$species,ignore.case = T));
    pred[ind,"superkingdom"]<- "Bacteriophage";
    }
    
    # ADD CSUM COL + exclude CSUM TAIL
    readn_sum<- sum(pred$readn);
    pred$csum<- (cumsum(pred$readn)-pred$readn)/readn_sum; # cumsum of preciding taxa
    ind<- which(pred$csum < (100-tail)/100.0);
    pred<- pred[ind,];


    # DELETE EXTRAC FIELDS: THESE WILL INTERFERE WITH MELT()>RECAST IF NOT NUMERICAL 
    pred<- pred[, c("readn","species","species_id","genus","genus_id","family","family_id","superkingdom","superkingdom_id")];
    true<- true[, c("readn","species","species_id","genus","genus_id","family","family_id","superkingdom","superkingdom_id")];

    tmp<- melt(pred, id=c("species","species_id","genus","genus_id","family","family_id","superkingdom","superkingdom_id"))
    tmp2<- cast(tmp,species+species_id+genus+genus_id+family+family_id+superkingdom+superkingdom_id~variable, sum)
    pred<- tmp2[order(tmp2$readn,decreasing=T),]    
    
    tmp<- melt(true, id=c("species","species_id","genus","genus_id","family","family_id","superkingdom","superkingdom_id"))
    tmp2<- cast(tmp,species+species_id+genus+genus_id+family+family_id+superkingdom+superkingdom_id~variable, sum)
    true<- tmp2[order(tmp2$readn,decreasing=T),]
    
      
    # saving true, pred, TP, FP and FN for Viruses and Bacteria
    res$vi$true<-   true[true$superkingdom=="Viruses",c("species","species_id","readn")];
    res$vi$pred<-   pred[pred$superkingdom=="Viruses",c("species","species_id","readn")];
    res$vi$TP<-     res$vi$pred[
                        match( intersect(res$vi$pred$species,res$vi$true$species), res$vi$pred$species), ];
    res$vi$FP<-     res$vi$pred[
                        match( setdiff(res$vi$pred$species, res$vi$true$species), res$vi$pred$species) , ];
    res$vi$FN<-     res$vi$true[
                        match( setdiff(res$vi$true$species,res$vi$pred$species), res$vi$true$species) , ];
    
    res$ba$true<-   true[true$superkingdom=="Bacteria",c("species","species_id","readn")];
    res$ba$pred<-   pred[pred$superkingdom=="Bacteria",c("species","species_id","readn")];
    res$ba$TP<-     res$ba$pred[
                        match(intersect(res$ba$pred$species,res$ba$true$species), res$ba$pred$species) , ];
    res$ba$FP<-     res$ba$pred[
                        match(setdiff(res$ba$pred$species,res$ba$true$species), res$ba$pred$species) , ];
    res$ba$FN<-     res$ba$true[
                        match(setdiff(res$ba$true$species,res$ba$pred$species), res$ba$true$species) , ];

    
    table<- matrix(0,nrow=3*2,ncol=10);
    table<- data.frame(table);
    colnames(table)<- c("superkingdom","rank","pred","true","TP","FP","FN","Pr","Rc","F1")
    king_list   <- c("Viruses","Bacteria");
    rank_list   <- c("species","genus","family");
    rowi        <- 1;
    for(i in 1:length(king_list)){
        table[rowi,"superkingdom"]<- king_list[i];
        
        for(j in 1:length(rank_list)){
            table[rowi,"rank"]<- rank_list[j];
            
            true_tmp<- unique( true[true$superkingdom==king_list[i],rank_list[j] ] );
            pred_tmp<- unique( pred[pred$superkingdom==king_list[i],rank_list[j] ] );
            table[rowi,"pred"]  <- length(pred_tmp);
            table[rowi,"true"]  <- length(true_tmp);
            table[rowi,"TP"]    <- length(intersect(pred_tmp,true_tmp));
            table[rowi,"FP"]    <- length( setdiff(pred_tmp,true_tmp));
            table[rowi,"FN"]    <- length( setdiff(true_tmp,pred_tmp));
            rowi = rowi+1;
        }
    }
    # add "TOTAL" rows to summary table
    total<- table[1:3,];
    total[1,1]<- "TOTAL";
    total[1:3,3:10]<- table[1:3,3:10]+table[4:6,3:10];
    table<- rbind(table,total);
    table[,"Pr"]<- table[,"TP"]/(table[,"TP"]+table[,"FP"]);
    table[,"Rc"]<- table[,"TP"]/(table[,"TP"]+table[,"FN"]);
    table$F1    <- (2*table$Pr*table$Rc)/(table$Pr+table$Rc);
    res$stats<- table;
    
    # if excel_file is given, write results to excel
    if(excel_file=="none"){ return(res);}

    wb<- createWorkbook();

    addWorksheet(wb, sheetName = "SUMMARY");
    addWorksheet(wb, sheetName = "vi_true");
    addWorksheet(wb, sheetName = "vi_pred");
    addWorksheet(wb, sheetName = "vi_TP");
    addWorksheet(wb, sheetName = "vi_FN");
    addWorksheet(wb, sheetName = "vi_FP");
    addWorksheet(wb, sheetName = "ba_true");
    addWorksheet(wb, sheetName = "ba_pred");
    addWorksheet(wb, sheetName = "ba_TP");
    addWorksheet(wb, sheetName = "ba_FN");
    addWorksheet(wb, sheetName = "ba_FP");

    writeDataTable(wb, sheet = "SUMMARY", x = res$stats, tableStyle = "None", withFilter = F);
    writeDataTable(wb, sheet = "vi_true", x = res$vi$true, tableStyle = "None", withFilter = F);
    writeDataTable(wb, sheet = "vi_pred", x = res$vi$pred, tableStyle = "None", withFilter = F);
    writeDataTable(wb, sheet = "vi_TP",   x = res$vi$TP, tableStyle = "None", withFilter = F);
    writeDataTable(wb, sheet = "vi_FN",   x = res$vi$FN, tableStyle = "None", withFilter = F);
    writeDataTable(wb, sheet = "vi_FP",   x = res$vi$FP, tableStyle = "None", withFilter = F);
    
    writeDataTable(wb, sheet = "ba_true", x = res$ba$true, tableStyle = "None", withFilter = F);
    writeDataTable(wb, sheet = "ba_pred", x = res$ba$pred, tableStyle = "None", withFilter = F);
    writeDataTable(wb, sheet = "ba_TP",   x = res$ba$TP, tableStyle = "None", withFilter = F);
    writeDataTable(wb, sheet = "ba_FN",   x = res$ba$FN, tableStyle = "None", withFilter = F);
    writeDataTable(wb, sheet = "ba_FP",   x = res$ba$FP, tableStyle = "None", withFilter = F);

    # set SUMMARY STYLE
    sty_pc  <- createStyle(numFmt = "0.0%");
    addStyle(wb, sheet = "SUMMARY", style = sty_pc,  rows = 1:(nrow(res$stats)+1), cols = c(8,9,10),gridExpand = T );
    
    # SET STYLE FOR THE REST
    sty_int <- createStyle(numFmt = "0");
    name_list <- names(wb);
    name_list <- setdiff(name_list,"SUMMARY");
    nrow_max  <- max(nrow(res$vi$true),nrow(res$vi$pred),nrow(res$ba$true),nrow(res$ba$pred));
    for(name in name_list){
        addStyle(wb, sheet = name, style = sty_int,  rows = 1:nrow_max, cols = c(2,3),gridExpand = T );
        setColWidths(wb, sheet = name, cols = 1:3, widths = "auto");
    }
    saveWorkbook(wb,excel_file,overwrite = T);

    return(res);
}


get_abund_tables<- function(abund_table,annot_table, excel_file="none", host_taxid=0, tail=1){
    #' @description  Reads raw abundance table: with readn+contign for taxa below species rank.
    #' And converts this to a abundance tables for Viruses/Bacteria at species/genus/family/superkingdom level
    #'
    #' @param abund_table file. Tab delimited file containing taxa abundancies. Format taxid+readn
    #' @param annot_table file. Tab delimited file containing contig annotations.
    #' @param excel_file file. Results are printed here.
    #' @param host_taxid num. Host taxid, host reads are excluded from read distribution and tail cutoff calculation.
    #' @param tail num. Number in percents [0,100], taxa with least abundance summing to this amount will be ignored
    #'

    library(reshape);
    library(openxlsx);
    res<- list();
    pred<-      read.table(file=abund_table,sep="\t",header=T,comment.char="",stringsAsFactors = F);
    
    pred[grep("^$",pred$species),"species"]<- "unknown";    # convert empty taxon labels to "unknown"
    pred[grep("^$",pred$genus),"genus"]<- "unknown";
    pred[grep("^$",pred$family),"family"]<- "unknown";
    pred[grep("^$",pred$superkingdom),"superkingdom"]<- "unknown";
    
    #pred<- replace(pred,is.na(pred),0);
    
    # DETECT PHAGES
    phage_families= c("Myoviridae","Siphoviridae","Podoviridae","Lipothrixviridae","Rudiviridae",
                    "Ampullaviridae","Bicaudaviridae","Clavaviridae","Corticoviridae","Cystoviridae",
                    "Fuselloviridae","Globuloviridae", "Guttaviridae", "Inoviridae","Leviviridae",
                    "Microviridae","Plasmaviridae","Tectiviridae");
    phage_orders = c("Caudovirales","Ligamenvirales");
    ind<- union(which(pred$family %in% phage_families),which(pred$order %in% phage_orders));
    ind<- union(ind,grep("\\<phage\\>|\\<bacteriophage\\>",pred$species,ignore.case = T));
    pred[ind,"superkingdom"]<- "Bacteriophage";
    
    # ADD PERCENTAGE COL
    readn_sum<- sum(pred$readn);
    pred$readn_pc<- pred$readn/readn_sum;
    
    # DELETE EXTRA FIELDS: THESE WILL INTERFERE WITH MELT()>RECAST IF NOT NUMERICAL 
    pred<- pred[, c("readn","readn_pc","contign","species","species_id","genus","genus_id","family","family_id","superkingdom","superkingdom_id")];
    

    # SUM OVER SPECIES/GENUS/FAMILY/SKINGDOM MEMBERS
    tmp<- melt(pred, id=c("species","species_id","genus","genus_id","family","family_id","superkingdom","superkingdom_id"));
    res$species <- cast(tmp,species+species_id+genus+genus_id+family+family_id+superkingdom+superkingdom_id~variable, sum);
    res$genus   <- cast(tmp,genus+genus_id+family+family_id+superkingdom+superkingdom_id~variable, sum);
    res$family  <- cast(tmp,family+family_id+superkingdom+superkingdom_id~variable, sum);
    res$sking   <- cast(tmp,superkingdom+superkingdom_id~variable, sum);
    
    
    # RECOUNT CONTIG NUMBER FROM ANNOT TABLE
    annot<-     read.table(file=annot_table,sep="\t",header=T,comment.char="",stringsAsFactors = F,quote="");
    annot<-     annot[, c("contig","species_id","genus_id","family_id","superkingdom_id")];
    for(i in 1:nrow(res$species)){
        res$species$contign[i] = sum( !is.na(unique( annot$contig[annot$species_id==res$species$species_id[i]]) ));
    }
    for(i in 1:nrow(res$genus)){
        res$genus$contign[i] = sum( !is.na(unique( annot$contig[annot$genus_id==res$genus$genus_id[i]]) ) );
    }
    for(i in 1:nrow(res$family)){
        res$family$contign[i] = sum( !is.na(unique( annot$contig[annot$family_id==res$family$family_id[i]]) ));
    }
    for(i in 1:nrow(res$sking)){
        res$sking$contign[i] = sum( !is.na(unique( annot$contig[annot$superkingdom_id==res$genus$superkingdom_id[i]])) );
    }
    
    # SORT
    name_list<- c("species","genus","family","sking");
    for(name in name_list){
        data<- res[[name]];
        data<- data[order(data$readn,decreasing = T),];
        res[[name]]<- data;
    }

    
    # CALC CUMSUM + DIVIDE CUMSUM INTO INTO CLASSES BY RANGE: 1:[0-95[, 2:[95-99[, 3:[99-100]
    cumsum_bounds<- c(1.0,0.99,0.95);
    cumsum_qlabels<- c(3,2,1);
    name_list<- c("species","genus","family","sking");
    host_ind<- which(pred$species_id==host_taxid);
    host_readn<- 0;
    if( length(host_ind)>0 ){
        host_readn<- pred[host_ind,'readn'];
    }
    for(name in name_list){
        res[[name]]$csum<- (cumsum(res[[name]]$readn) - res[[name]]$readn - host_readn)/ (readn_sum - host_readn); # cumsum of preciding taxa
        res[[name]]$csumq<- cumsum_qlabels[1];
        for(i in 2:length(cumsum_bounds)){
            res[[name]]$csumq[res[[name]]$csum < cumsum_bounds[i]] <- cumsum_qlabels[i];
        }
    }
    
    # REORDER COLS
    cols_pref<- c("readn","readn_pc","csum","csumq","contign");
    res$species <- res$species[,c(cols_pref,"species","species_id","genus","genus_id","family","family_id","superkingdom","superkingdom_id")];
    res$genus   <- res$genus[,c(cols_pref,"genus","genus_id","family","family_id","superkingdom","superkingdom_id")];
    res$family  <- res$family[,c(cols_pref,"family","family_id","superkingdom","superkingdom_id")];
    res$sking   <- res$sking[,c(cols_pref,"superkingdom","superkingdom_id")];
   
    
    # RETRIEVE SUBGROUPS
    res$vi_species  <- res$species[res$species$superkingdom=="Viruses",c(cols_pref,"species","species_id","genus","family")];
    res$vi_genus    <- res$genus[res$genus$superkingdom=="Viruses",c(cols_pref,"genus","genus_id","family")];
    res$vi_family   <- res$family[res$family$superkingdom=="Viruses",c(cols_pref,"family","family_id")];
    
    res$ba_species  <- res$species[res$species$superkingdom=="Bacteria",c(cols_pref,"species","species_id","genus","family")];
    res$ba_genus    <- res$genus[res$genus$superkingdom=="Bacteria",c(cols_pref,"genus","genus_id","family")];
    res$ba_family   <- res$family[res$family$superkingdom=="Bacteria",c(cols_pref,"family","family_id")];
    
    res$ph_species  <- res$species[res$species$superkingdom=="Bacteriophage",c(cols_pref,"species","species_id","genus","family")];
    res$ph_genus    <- res$genus[res$genus$superkingdom=="Bacteriophage",c(cols_pref,"genus","genus_id","family")];
    res$ph_family   <- res$family[res$family$superkingdom=="Bacteriophage",c(cols_pref,"family","family_id")];

    res$eu_species  <- res$species[res$species$superkingdom=="Eukaryota",c(cols_pref,"species","species_id","genus","family")];
    res$eu_genus    <- res$genus[res$genus$superkingdom=="Eukaryota",c(cols_pref,"genus","genus_id","family")];
    res$eu_family   <- res$family[res$family$superkingdom=="Eukaryota",c(cols_pref,"family","family_id")];
  
        
    # ADD ROW NUMS FOR CLARITY (comment out to preserv row numbers from the original pred_file)
    name_list<- names(res);
    for(name in name_list){
        if(nrow(res[[name]]) > 0){
            rownames(res[[name]])<- c(1:nrow(res[[name]]));
        }
    }
 
    # WRITE TO EXCEL FILE  (IF GIVEN)
    if(excel_file=="none"){ return(res);}
    wb<- createWorkbook();
    name_list <- c("vi_family","vi_genus","vi_species","ba_family","ba_genus","ba_species","ph_family","ph_genus","ph_species",
                    "eu_family","eu_genus","eu_species","sking","family","genus");
    sty_int <- createStyle(numFmt = "0");
    sty_pc  <- createStyle(numFmt = "0.0%");
    sty_tail <- createStyle(fontColour = "#808080");

    for(name in name_list){
        if(nrow(res[[name]]) > 0){
            addWorksheet(wb,sheetName=name);
            writeDataTable(wb,sheet=name,x= res[[name]],tableStyle = "TableStyleMedium4");
            # STYLE:
            M<- nrow(res[[name]]) + 1;
            N<- ncol(res[[name]]);
            tail_rows<- which(res[[name]]$csum >= ((100-tail)/100.0) ) + 1;
            addStyle(wb, sheet = name, style = sty_int, rows = 2:M, cols = 1);
            addStyle(wb, sheet = name, style = sty_pc,  rows = 2:M, cols = 2);
            addStyle(wb, sheet = name, style = sty_pc,  rows = 2:M, cols = 3);
            setColWidths(wb, sheet = name, cols = 1:N, widths = "auto");
            
            if( (name!= "sking") && (length(tail_rows) > 0)){ # cutting tail does not work for taxa at highest rank
                addStyle(wb, sheet = name, style = sty_tail,rows = tail_rows, cols = c(1:N),gridExpand = T,stack = T);
            }
        }
    }
    saveWorkbook(wb,excel_file,overwrite = T);
    return(res);
}

# GENERATES SEQUENCE ANNOTATION TABLES
# 
# Parameters:
# annot_file       : path to the annotation file
# excel_file       : Write generated annotation tables as data sheets to this file (overwrite if exists).
# append_prefix    : Append generated annotation tables as data sheets to the existing excel_file.
#                    Appended data sheets are identified with the specified prefix.
#                    If no excel_file exists, create the file and append data sheets as above.
get_annot_tables<- function(annot_file,excel_file="none",append_prefix=""){
    library(openxlsx);
    res<- list();
    pred<- read.table(file=annot_file,sep="\t",header=T,comment.char="",stringsAsFactors = F,quote="");        
    
    # CHANGE HEADERS: MEGAHIT uses NODE/flag/multi/len headers
    headers<- colnames(pred);
    headers<- sub("^NODE$","contig",headers,ignore.case = T);
    headers<- sub("^multi$","coverage",headers,ignore.case = T);
    headers<- sub("^len$","length",headers,ignore.case = T);
    # CHANGE HEADERS FOR VELVET
    headers<- sub("^contig_cov$","coverage",headers,ignore.case = T);
    headers<- sub("^contig_len$","length",headers,ignore.case = T);
    colnames(pred)<- headers;
    
    # convert empty taxon labels to "unknown"
    pred[grep("^$",pred$species),"species"]<- "unknown";    
    pred[grep("^$",pred$genus),"genus"]<- "unknown";
    pred[grep("^$",pred$family),"family"]<- "unknown";
    pred[grep("^$",pred$order),"order"]<-   "unknown";
    pred[grep("^$",pred$superkingdom),"superkingdom"]<- "unknown";
    
    # delete genus_id,family_id,order_id,superkingdom_id cols
    delete_list<- c("genus_id","family_id","order_id","superkingdom_id");
    pred<- pred[, !(colnames(pred) %in% delete_list)];
    
    #pred<- replace(pred,is.na(pred),0);
    
    # DETECT PHAGES
    phage_families= c("Myoviridae","Siphoviridae","Podoviridae","Lipothrixviridae","Rudiviridae",
                    "Ampullaviridae","Bicaudaviridae","Clavaviridae","Corticoviridae","Cystoviridae",
                    "Fuselloviridae","Globuloviridae", "Guttaviridae", "Inoviridae","Leviviridae",
                    "Microviridae","Plasmaviridae","Tectiviridae");
    phage_orders = c("Caudovirales","Ligamenvirales");
    #levels(pred$superkingdom)<- c(levels(pred$superkingdom),"Bacteriophage");
    ind<- union(which(pred$family %in% phage_families),which(pred$order %in% phage_orders));
    ind<- union(ind,grep("\\<phage\\>|\\<bacteriophage\\>",pred$species,ignore.case = T));
    pred[ind,"superkingdom"]<- "Bacteriophage";
    
    
    # SELECT SUBSETS: Viruses/Bacteria/..
    tmp_contigs<- unique( pred[pred$superkingdom=="Viruses","contig"] );
    tmp_ind<- which(pred$contig %in% tmp_contigs);
    res$vi<- pred[tmp_ind,];

    tmp_contigs<- unique( pred[pred$superkingdom=="Bacteria","contig"] );
    tmp_ind<- which(pred$contig %in% tmp_contigs);
    res$ba<- pred[tmp_ind,]; 
    
    tmp_contigs<- unique( pred[pred$superkingdom=="Bacteriophage","contig"] );
    tmp_ind<- which(pred$contig %in% tmp_contigs);
    res$ph<- pred[tmp_ind,];
    
    tmp_contigs<- unique( pred[pred$superkingdom=="Eukaryota","contig"] );
    tmp_ind<- which(pred$contig %in% tmp_contigs);
    res$eu<- pred[tmp_ind,];
    
    tmp_conts<- unique( pred[pred$superkingdom=="unknown","contig"] );
    tmp_conts_flt<- list();
    for(cont in tmp_conts){
        tmp<- pred[pred$contig==cont,];
        annotated<- F;
        for(i in 1:nrow(tmp)){
            if( tmp[i,"superkingdom"] != "unknown"){
                annotated<- T;
                break;
            }
        }
        if(!annotated){
            tmp_conts_flt<- c(tmp_conts_flt,cont);
        }
    }
    tmp_conts_flt<- as.character(tmp_conts_flt);
    tmp_ind<- which(pred$contig %in% tmp_conts_flt);
    res$un<- pred[tmp_ind,];
    
    # SORT
    if(nrow(res$vi) > 0){ res$vi<- res$vi[order(res$vi$length, res$vi$contig, decreasing=T),] }
    if(nrow(res$ba) > 0){ res$ba<- res$ba[order(res$ba$length, res$ba$contig, decreasing=T),] }
    if(nrow(res$ph) > 0){ res$ph<- res$ph[order(res$ph$length, res$ph$contig, decreasing=T),] }
    if(nrow(res$eu) > 0){ res$eu<- res$eu[order(res$eu$length, res$eu$contig, decreasing = T),] }
    if(nrow(res$un) > 0){ res$un<- res$un[order(res$un$length, res$un$contig, decreasing = T),] }
    
    # PRINT
    if(excel_file=="none"){ return(res);}
    
    wb<- createWorkbook();
    sheet_name_prefix = "";
    if( append_prefix !=""){
        sheet_name_prefix = paste(append_prefix,".",sep="");
    }
    if( file.exists(excel_file) && (append_prefix !="") ){
        wb<- loadWorkbook(file=excel_file);
    }
    
    if(nrow(res$vi)> 0){
        sheet_name = paste(sheet_name_prefix,"Viruses",sep="");
        if(sheet_name %in% names(wb)){ removeWorksheet(wb,sheet_name); }
        addWorksheet(wb,sheetName = sheet_name);
        writeDataTable(wb,sheet = sheet_name,x = res$vi, tableStyle = "TableStyleMedium4");}
    if(nrow(res$ba) > 0){
        sheet_name = paste(sheet_name_prefix,"Bacteria",sep="");
        if(sheet_name %in% names(wb)){ removeWorksheet(wb,sheet_name); }
        addWorksheet(wb,sheetName = sheet_name);
        writeDataTable(wb,sheet = sheet_name,x = res$ba, tableStyle = "TableStyleMedium4");}
    if(nrow(res$ph) > 0){
        sheet_name = paste(sheet_name_prefix,"Bacteriophage",sep="");
        if(sheet_name %in% names(wb)){ removeWorksheet(wb,sheet_name); }
        addWorksheet(wb,sheetName = sheet_name);
        writeDataTable(wb,sheet = sheet_name,x = res$ph, tableStyle = "TableStyleMedium4");}
    if(nrow(res$eu) > 0){
        sheet_name = paste(sheet_name_prefix,"Eukaryota",sep="");
        if(sheet_name %in% names(wb)){ removeWorksheet(wb,sheet_name); }
        addWorksheet(wb,sheetName = sheet_name);
        writeDataTable(wb,sheet = sheet_name,x = res$eu, tableStyle = "TableStyleMedium4");}
    if(nrow(res$un) > 0){
        sheet_name = paste(sheet_name_prefix,"Unknown",sep="");
        if(sheet_name %in% names(wb)){ removeWorksheet(wb,sheet_name); }
        addWorksheet(wb,sheetName = sheet_name);
        writeDataTable(wb,sheet = sheet_name,x = res$un, tableStyle = "TableStyleMedium4");}   
    saveWorkbook(wb,excel_file,overwrite = T);
    
    return(res);
}

