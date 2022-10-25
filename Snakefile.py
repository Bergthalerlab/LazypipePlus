#
# LAZYPIPE2: NGS PIPELINE FOR VIRUS DISCOVERY AND METAGENOMICS
#
# SNAKEMAKE INTERFACE
#
# Credit:
# Plyusnin,I., Kant,R., Jaaskelainen,A.J., Sironen,T., Holm,L., Vapalahti,O. and Smura,T. (2020) 
# Novel NGS Pipeline for Virus Discovery from a Wide Spectrum of Hosts and Sample Types. Virus Evolution, veaa091
#
# Contact: grp-lazypipe@helsinki.fi
#

import os
import glob

configfile:     "config.yaml"

# expand env variables
for k in config.keys():
    val = config[k]
    if(type(val) == str):
        val = os.path.expandvars(val)
        config[k] = val
    if(type(val) == dict):
        for k2 in val.keys():
            val2 = val[k2]
            if(type(val2) == str):
                val2 = os.path.expandvars(val2)
                val[k2] = val2
        config[k] = val



def get_all_output_files(wildcards):
    files    = []
    for sample in config["datain"].keys():
        files.append(config["res"]+"/"+sample+".tar.gz")
    files.append(config["res"]+"/"+sample+"/coverage_plot.png")
    return files
        
rule all:
    input:
        get_all_output_files


rule prepro_reads:
    input:
        r1     = lambda wildcards: config["datain"][wildcards.sample],
        r2     = lambda wildcards: config["datain"][wildcards.sample].replace("R1.","R2.")
    output:
        p1     = config["res"]+"/{sample}/trimmed_paired1.fq",
        p2     = config["res"]+"/{sample}/trimmed_paired2.fq",
        up1    = config["res"]+"/{sample}/trimmed_unpaired1.fq",
        up2    = config["res"]+"/{sample}/trimmed_unpaired2.fq"
    threads:
        config["threads_max"]
    params:
        fastp    = config["fastp_par"],
        trimm    = config["trimm_par"]
    log:
        config["logs"]+"/{sample}/prepro_reads.log"
    run:
        if config["pre"]=="fastp":
            shell("fastp --thread {threads} -i {input.r1} -I {input.r2} "
                  "-o {output.p1} -O {output.p2} --unpaired1 {output.up1} --unpaired2 {output.up2} "
                  "{params.fastp} 2> {log}")
        elif config["pre"]=="trimm":
            shell("trimmomatic PE -threads {threads} {input.r1} {input.r2} {output.p1} {output.up1} {output.p2} {output.up2} "
                  "{params.trimm} &> {log}")
        elif config["pre"]=="none":
            shell("cp {input.r1} {output.p1}")
            shell("cp {input.r2} {output.p2}")
        else:
            shell("fastp --thread {threads} -i {input.r1} -I {input.r2} "
                  "-o {output.p1} -O {output.p2} --unpaired1 {output.up1} --unpaired2 {output.up2} "
                  "{params.fastp} 2> {log}")


rule filter_hostgen:
    input:
        hostgen     = lambda wildcards: config["hostgen_sm"][wildcards.sample],
        read1       = config["res"]+"/{sample}/trimmed_paired1.fq",
        read2       = config["res"]+"/{sample}/trimmed_paired2.fq"
    output:
        sam         = config["res"]+"/{sample}/hostgen.sam",
        sam_flt     = config["res"]+"/{sample}/hostgen.sam.flt",
        bam         = config["res"]+"/{sample}/hostgen.bam",
        bam_sorted  = config["res"]+"/{sample}/hostgen.sorted.bam",
        readids     = config["res"]+"/{sample}/hostgen.readids",
        read1       = config["res"]+"/{sample}/trimmed_paired1_hostflt.fq",
        read2       = config["res"]+"/{sample}/trimmed_paired2_hostflt.fq"
    log:
        config["logs"]+"/{sample}/filter_hostgen.log"
    threads:
        config["threads_max"]
    params:
        bwa         = "-T {}".format(config["hostgen_flt_th"]),
        bwaindex    = "-a bwtsw -b 50000000",
        wrkdir      = config["wrkdir"]
    run:
        hostgen     = config["hostgen_sm"][wildcards.sample]
        if not os.path.isfile( hostgen+".amb" ):    # index hostgen if needed
            shell("bwa index {params.bwaindex} {input.hostgen} &> {log}")
        shell("bwa mem -t {threads} {params.bwa} {input.hostgen} {input.read1} {input.read2} 1> {output.sam} 2> {log}")
        shell("samtools view -@ {threads} -F4 -Shb {output.sam} 1> {output.bam} 2>> {log}")
        shell("samtools sort -@ {threads} -n {output.bam} -T {params.wrkdir} -o {output.bam_sorted} 2>> {log}")
        shell("samtools view -F4 {output.bam_sorted} | perl perl/filter_tophits.pl 1> {output.sam_flt} 2>> {log}")
        shell("cat {output.sam_flt} | cut -f1 | sort -T {params.wrkdir} | uniq 1> {output.readids} 2>> {log}")
        shell("bin/filtfq -1 {input.read1} -2 {input.read2} -o {output.read1} -O {output.read2} -f {output.readids} -m filter -v &>> {log}")

# this rule aligns reads to the virousaurus reference and produces coverage plots
rule coverage_plots:
	input:
		fasta1       = config["res"]+"/{sample}/trimmed_paired1_hostflt.fq",
		fasta2       = config["res"]+"/{sample}/trimmed_paired2_hostflt.fq"
	output:
		sam           = temp(config["res"]+"/{sample}/virosaurus_mapped.sam"),
		bam_sorted    = temp(config["res"]+"/{sample}/virosaurus_mapped_sorted.bam"),
		uniqmap_bam   = config["res"]+"/{sample}/virosaurus_mapped_sorted_uniquemap.bam",
		uniqmap_tsv   = config["res"]+"/{sample}/virosaurus_mapped_sorted_uniquemap.tsv",
		allmap_tsv   = config["res"]+"/{sample}/virosaurus_mapped_sorted_allmap.tsv",
		coverage_plot_genes = config["res"]+"/{sample}/genes_coverage_plot.png"
		coverage_plot_genomes = config["res"]+"/{sample}/genomes_coverage_plot.png"
	threads:
		config["threads_max"]
	log:
		config["logs"]+"/{sample}/coverage_plots.log"
	params:
		bwa         = "-T {}".format(config["hostgen_flt_th"]),
		bwaindex    = "-a bwtsw -b 50000000",
		wrkdir      = config["wrkdir"],
		sample_name = "{sample}"
		output_path = config["res"]+"/{sample}/"
	run:
		virosaurus_fasta = config["virosaurus_fasta"]
		if not os.path.isfile(virosaurus_fasta+".amb" ):   # index virosaurus if needed
			shell("bwa index {params.bwaindex} {virosaurus_fasta} &> {log}")
		shell("bwa mem -t {threads} {virosaurus_fasta} {input.fasta1} {input.fasta2} 1> {output.sam} 2>> {log}")
		shell("samtools sort -O bam {output.sam} 1> {output.bam_sorted} 2>> {log}")
		shell("samtools index {output.bam_sorted} 2>> {log}")
		shell("samtools view -h {output.bam_sorted} | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b 1> {output.uniqmap_bam} 2>> {log}")
		shell("samtools depth -Q 20 {output.uniqmap_bam} 1> {output.uniqmap_tsv} 2>> {log}")
		shell("samtools depth -Q 20 {output.bam_sorted} 1> {output.allmap_tsv} 2>> {log}")
		shell("python python/plot.py {output.uniqmap_tsv} {output.allmap_tsv} GENOMES {params.output_path} 2>> {log}")
		shell("python python/plot.py {output.uniqmap_tsv} {output.allmap_tsv} GENES {params.output_path} 2>> {log}")

        
rule assemble:
    input:
        r1 = lambda wildcards: config["res"]+"/{sample}/trimmed_paired1_hostflt.fq" if bool(config["hostgen_sm"][wildcards.sample]) else config["res"]+"/{sample}/trimmed_paired1.fq",
        r2 = lambda wildcards: config["res"]+"/{sample}/trimmed_paired2_hostflt.fq" if bool(config["hostgen_sm"][wildcards.sample]) else config["res"]+"/{sample}/trimmed_paired2.fq"
    output:
        assembler_out   = directory(config["res"]+"/{sample}/assembler_out"),
        contigs         = config["res"]+"/{sample}/contigs.fa"
    log:
        config["logs"]+"/{sample}/assemble.log"
    threads:
        config["threads_max"]
    params:
        megahit   = "",
        spades    = ""
    run:
        #shell("rm -fr {output.assembler_out}")
        if config["ass"]=="megahit":
            shell("megahit -t {threads} -1 {input.r1} -2 {input.r2} --out-dir {output.assembler_out} &> {log}")
            shell("cat {output.assembler_out}/final.contigs.fa | "
                  "sed 's/>\\([[:alnum:]]\\+\\)_\\([0-9]\\+\\)/>contig=\\1.\\2/' | "
                  "sed 's/multi/coverage/' | "
                  "sed 's/len/length/' | "
                  "sed 's/\\s\\+/_/g' 1> {output.contigs} 2>> {log}")
        elif config["ass"]=="spades":
            shell("spades.py --meta -t {threads} -1 {input.r1} -2 {input.r2} -o {output.assembler_out} &> {log}")
            shell("cat {output.assembler_out}/contigs.fasta | "
                  "sed 's/>NODE_/>contig=/' | "
                  "sed 's/length_/length=/' | "
                  "sed 's/cov_/coverage=/' 1> {output.contigs} 2>> {log} ")
            shell("cp {output.assembler_out}/scaffolds.fasta results/{wildcards.sample}/scaffolds.fa")
        else:
            os.write(2,"ERROR: unknown assembler: "+config["ass"]+"\n")
            
            
rule realign_reads:
    input:
        contigs   = config["res"]+"/{sample}/contigs.fa",
        r1        = lambda wildcards: config["res"]+"/{sample}/trimmed_paired1_hostflt.fq" if bool(config["hostgen_sm"][wildcards.sample]) else config["res"]+"/{sample}/trimmed_paired1.fq",
        r2        = lambda wildcards: config["res"]+"/{sample}/trimmed_paired2_hostflt.fq" if bool(config["hostgen_sm"][wildcards.sample]) else config["res"]+"/{sample}/trimmed_paired2.fq",
    output:
        sam       = temp(config["res"]+"/{sample}/contigs.bwa.sam"),
        bam_top   = temp(config["res"]+"/{sample}/contigs.top.bam"),
        idxstats  = config["res"]+"/{sample}/contigs.idxstats",
        idmap     = config["res"]+"/{sample}/readid_contigid.tsv"
    log:
        config["logs"]+"/{sample}/realign_reads.log"
    threads:
        config["threads_max"]        
    params:
        bwa         = "-T {}".format(config["realign_read_th"]),
        bwaindex     = "",        # "-a bwtsw -b 50000000"
        wrkdir      = config["wrkdir"]
    run:
        shell("bwa index {params.bwaindex} {input.contigs} &> {log}")
        shell("bwa mem -t {threads} {params.bwa} {input.contigs}  {input.r1} {input.r2} 1> {output.sam} 2>> {log}")
        shell("samtools sort -@ {threads} -n -T {params.wrkdir} {output.sam} | "
              "samtools view -F4 -h | perl perl/filter_tophits.pl | "
              "samtools sort -@ {threads} -T {params.wrkdir} -o {output.bam_top} 2>> {log}")
        shell("samtools index -@ {threads} {output.bam_top} 2>> {log}")
        shell("samtools idxstats {output.bam_top} 1> {output.idxstats} 2>> {log}")   
        shell("samtools view  {output.bam_top} | cut -f1,3 1> {output.idmap} 2>> {log}") 
        
        
rule detect_orfs:
    input:
        contigs     = config["res"]+"/{sample}/contigs.fa"
    output:
        orfs_raw    = temp(config["res"]+"/{sample}/ORFs.raw"),
        orfs_gtf    = config["res"]+"/{sample}/ORFs.gtf",
        orfs_nt     = config["res"]+"/{sample}/ORFs.nt.fa",
        orfs_aa     = config["res"]+"/{sample}/ORFs.aa.fa"
    log:
        config["logs"]+"/{sample}/detect_orfs.log"
    threads:
        config["threads_max"]
    params:
        "-w 90 --min-len {}".format(config["min_gene_length"])
    run:
        shell("mga {input.contigs} -m 1> {output.orfs_raw} 2> {log}")
        shell("perl perl/mga2gtf.pl {output.orfs_raw} 1> {output.orfs_gtf} 2> {log}")
        shell("rm -f {input.contigs}.seqkit.fai")
        shell("seqkit subseq --gtf {output.orfs_gtf} {input.contigs} | "
	      "seqkit seq {params} | "
	      "seqkit replace -p '_([\w\:+\-]+)\s*$' -r '_ORF=$1' 1> {output.orfs_nt} 2>> {log}")
	shell("seqkit translate -w 90 -f 1 {output.orfs_nt} 1> {output.orfs_aa} 2> {log}")

rule sans_orfs:
    input:
        orfs          = config["res"]+"/{sample}/ORFs.aa.fa"
    output:
        dbhits        = config["res"]+"/{sample}/dbhits.sans.tsv"
    log:
        config["logs"]+"/{sample}/sans_orfs.log"
    threads:
        config["threads_max"]
    params:
        sans          = "-m SANStopHtaxid --SANS_H 5 -R "
    shell:
        "python3.10 SANSPANZ.3/runsanspanz.py {params.sans} -i {input.orfs} -o {output.dbhits} &> {log}"
        
rule blastp_orfs:
    input:
        config["res"]+"/{sample}/ORFs.aa.fa",
    output:
        config["res"]+"/{sample}/dbhits.blastp.tsv"
    log:
        config["logs"]+"/{sample}/blastp_orfs.log"
    threads:
        config["threads_max"]
    params:
        "-evalue 10 -use_sw_tback -max_target_seqs 5 -outfmt '6 qseqid sacc bitscore pident length stitle sscinames staxid ' "+
        "-db "+config["blastp_db"]
    shell:
        "blastp -num_threads {threads} {params} -query {input} | "+
        "csvtk add-header -I -t -n qseqid,sacc,bitscore,pident,length,stitle,sscinames,staxid 1> {output} 2> {log}"
        
rule blastn_contigs_vidb:
    input:
        config["res"]+"/{sample}/{contigs}.fa"
    output:
        config["res"]+"/{sample}/{contigs}.vidb.hits.tsv"
    log:
        config["logs"]+"/{sample}/blastn_{contigs}.vidb.log"
    threads:
        config["threads_max"]
    params:
        "-evalue 10 -max_target_seqs 5  -outfmt '6 qseqid sseqid bitscore pident length stitle sscinames staxid ' "+
        "-db "+config["blastn_vi_db"]
    shell:
        "blastn -num_threads {threads} {params} -query {input} | "+
        "csvtk add-header -I -t -n qseqid,sseqid,bitscore,pident,length,stitle,sscinames,staxid 1> {output} 2> {log}"    


rule centrifuge_contigs_habv:
    input:
        config["res"]+"/{sample}/contigs.fa"
    output:
        dbhits    = config["res"]+"/{sample}/dbhits.cent.tsv",
        report    = temp(config["res"]+"/{sample}/dbhits.cent.report")
    log:
        config["logs"]+"/{sample}/centrifuge_contigs_habv.log"
    threads:
        config["threads_max"]
    params:
        "-f --mm -x "+config["centrifuge_db"] 
    shell:
        "centrifuge --threads {threads} {params} -U {input} --report-file {output.report} -S {output.dbhits} &> {log}"

        
def select_dbhits_file(wildcards):
    if config["ann"]=="sans":
        return config["res"]+"/{sample}/dbhits.sans.tsv"
    elif config["ann"]=="blastp":
        return config["res"]+"/{sample}/dbhits.blastp.tsv"
    elif config["ann"]=="centrifuge":
        return config["res"]+"/{sample}/dbhits.cent.tsv"

rule format_hsearch_main:
    input:
        hits_raw        = select_dbhits_file
    output:
        hits_srt        = config["res"]+"/{sample}/dbhits.srt.tsv",            # db search results: sorted/edited
        hits_top        = config["res"]+"/{sample}/dbhits.top.tsv",            # db search results: top hits
        contigs_annot   = config["res"]+"/{sample}/contig_taxid_score.tsv"      # tab-separated table with fields: contigid + taxid + score
    log:
        config["logs"]+"/{sample}/format_hsearch_main.log"
    params:
        min_sans_bits    = config["min_sans_bits"],
        min_blasp_bits    = config["min_blastp_bits"],
        min_cent_bits    = config["min_cent_bits"]
    run:
        if config["ann"]== "sans":
            shell("cat {input.hits_raw} | "
                  "csvtk rename -t -f qpid,spid -n qseqid,sseqid | "
                  "csvtk rename -t -f bits,taxid -n bitscore,staxid | "
                  "csvtk filter -t -f 'bitscore>{params.min_sans_bits}' | "
                  "csvtk cut -t -f qseqid,sseqid,bitscore,qcov,scov,pide,lali,desc,staxid 1> {output.hits_srt} 2>> {log}")
            shell("perl perl/filter_tophits2.pl --qcol 1 --bitscol 3 --ties {output.hits_srt} 1> {output.hits_top} 2>> {log}")
            shell("cat {output.hits_top} | "
                  "csvtk cut -t -f qseqid,staxid,bitscore | "
                  "csvtk replace -t -f qseqid -p '^contig=([A-Za-z0-9\.]+)_(.+)' -r'$1' 1> {output.contigs_annot} 2>> {log}")
            
        if config["ann"]== "blastp":
            shell("cat {input.hits_raw} | "
                  "csvtk filter -t -f 'bitscore>{params.min_blastp_bits}' | "
                  "csvtk sort -t -k qid:N -k bitscore:nr 1> {output.hits_str} 2> {log}") 
            shell("perl perl/filter_tophits2.pl --qcol 1 --bitscol 3 --ties {output.hits_srt} 1> {output.hits_top} 2>> {log}")
            shell("cat {output.hits_top} | "
                  "csvtk cut -t -f qseqid,staxid,bitscore | "
                  "csvtk replace -t -f qseqid -p '^contig=([A-Za-z0-9\.]+)_(.+)' -r'$1' 1> {output.contigs_annot} 2>> {log}")
        
        if config["ann"]=="centrifuge":
            shell("cat {input.hits_raw} | "
                  "csvtk rename -t -f readID,seqID,taxID,score,hitLength,queryLength  -n qseqid,sseqid,staxid,bitscore,slen,qlen | "
                  "csvtk filter -t -f 'bitscore>={params.min_cent_bits}' | "
                  "csvtk sort -t -k qseqid:N -k bitscore:nr | "
                  "csvtk cut -t -f qseqid,sseqid,bitscore,qlen,slen,staxid 1> {output.hits_srt} 2>> {log}")
            shell("perl perl/filter_tophits2.pl --qcol 1 --bitscol 3 --ties {output.hits_srt} 1> {output.hits_top} 2>> {log}")
            shell("cat {output.hits_top} | "
                  "csvtk cut -t -f qseqid,staxid,bitscore | "
                  "csvtk replace -t -f qseqid -p '^contig=([A-Za-z0-9\.]+)_(.+)' -r'$1' 1> {output.contigs_annot} 2>> {log}")
        

rule create_reports:
    input:
        contigs             = config["res"]+"/{sample}/contigs.fa",
        contig_taxid_score  = config["res"]+"/{sample}/contig_taxid_score.tsv",
        contigs_stats       = config["res"]+"/{sample}/contigs.idxstats",
        hits_top            = config["res"]+"/{sample}/dbhits.top.tsv"
    output:
        readn_taxid         = config["res"]+"/{sample}/readn_taxid.tsv",
        abund_table         = config["res"]+"/{sample}/abund_table.tsv",
        annot_table         = config["res"]+"/{sample}/contigs.annot.tsv",
        taxprofile          = config["res"]+"/{sample}/taxprofile.txt",
        contigs_dir         = directory(config["res"]+"/{sample}/contigs") # sorted contigs
    log:
        config["logs"]+"/{sample}/generate_reports.log"
    params:
        abund_table     = "-h -w "+ config["weights"] +" --cont_score_tail "+str(config["cont_score_tail"]),
        taxprofile      = "--sample {sample} --tail "+str(config["tail"]),
        hosttaxid       = "0",
        taxonomy        = config["taxonomy"]
    threads:
        min(8,config["threads_max"])
    run:
        # UPDATE TAXONOMY
        if config["taxonomy_update"]:
            shell("perl perl/update_taxonomy.pl")
	
        #shell("rm -fr {output.contigs_dir}")		
    
        # ABUNDANCE AND CONTIG ANNOT TABLES
        shell("perl perl/get_abund_table.pl {params.abund_table} {input.contig_taxid_score} {input.contigs_stats} 1> {output.readn_taxid} 2> {log}")
        shell("csvtk del-header -t {output.readn_taxid} | "
              "taxonkit reformat --data-dir {params.taxonomy} -j {threads} -I 3 -t -f '{{s}}\\t{{g}}\\t{{f}}\\t{{o}}\\t{{c}}\\t{{p}}\\t{{k}}' | "
              "csvtk add-header -I -t -n readn,contign,taxid,species,genus,family,order,class,phylum,superkingdom,species_id,genus_id,family_id,order_id,class_id,phylum_id,superkingdom_id "
              "1> {output.abund_table} 2>> {log}")
        
        shell("csvtk del-header -t {input.hits_top} | "
              "taxonkit reformat --data-dir {params.taxonomy} -j {threads} -I `csvtk ncol -t {input.hits_top}` -t -f '{{s}}\\t{{g}}\\t{{f}}\\t{{o}}\\t{{c}}\\t{{p}}\\t{{k}}' | "
              "csvtk add-header -I -t -n `head -n 1 {input.hits_top} | tr '\\t' ','`,species,genus,family,order,class,phylum,superkingdom,species_id,genus_id,family_id,order_id,class_id,phylum_id,superkingdom_id | "
              "perl perl/split_qid.pl --header 1> {output.annot_table} 2>> {log}");

        # TAXONOMIC PROFILE
        shell("perl perl/abundtable2taxprofile.pl {params.taxprofile} {output.abund_table} 1> {output.taxprofile} 2>> {log}")
        
        # SORT CONTIGS TO DIRS USING TAXONOMY CLASSIFICATION
        shell("perl perl/sort_contigs_bytaxa.pl --res {output.contigs_dir} --ann {output.annot_table} --cont {input.contigs} &>> {log}")


rule format_hsearch_contigs_vidb:
    input:
        config["res"]+"/{sample}/{contigs}.vidb.hits.tsv"
    output:
        contigs_annot    = config["res"]+"/{sample}/{contigs,[A-Za-z0-9]+_[A-Za-z0-9]+}.annot.tsv" 
    log:
        config["logs"]+"/{sample}/format_hsearch_{contigs}.log"
    threads:
        config["threads_max"]
    params:
        taxonomy        = config["taxonomy"]
    run:
        # UPDATE TAXONOMY
        if config["taxonomy_update"]:
            shell("perl perl/update_taxonomy.pl")
            
        #shell("perl perl/filter_tophits2.pl --qcol 1 --bitscol 4 --ties {input} 1> {output.hits_top} 2> {log}")
        shell("csvtk del-header -t {input} | "
              "taxonkit reformat --data-dir {params.taxonomy} -j {threads} -I `csvtk ncol -t {input}` -f '{{s}}\\t{{g}}\\t{{f}}\\t{{o}}\\t{{c}}\\t{{p}}\\t{{k}}' | "
              "csvtk add-header -I -t -n `head -n 1 {input} | tr '\\t' ','`,species,genus,family,order,class,phylum,superkingdom | "
              "perl perl/split_qid.pl --header 1> {output.contigs_annot} 2>> {log}")
            
            

rule extract_vicontigs:
    input:
        contigs         = config["res"]+"/{sample}/contigs.fa",
        contigs_annot   = config["res"]+"/{sample}/contigs.annot.tsv"
    output:
        contigids       = temp(config["res"]+"/{sample}/contigs_vi.ids"),
        contigs         = config["res"]+"/{sample}/contigs_vi.fa"
    log:
        config["logs"]+"/{sample}/extract_vicontigs.log"
    run:
        shell("cat {input.contigs_annot} | "
              "csvtk filter2 -t -f '$superkingdom==\"Viruses\"' | "
              "csvtk cut -t -f contig | uniq 1> {output.contigs}.tmp1 2> {log}")
        shell("sed -n \"/^>/s/^>// p\" {input.contigs} | "
              "csvtk add-header -I -t -n qseqid | csvtk mutate -t -f qseqid -n contig -p '^contig=([A-Za-z0-9\.]+)_(.+)' 1> {output.contigs}.tmp2 2>> {log}")
        shell("csvtk join -t -f contig {output.contigs}.tmp1 {output.contigs}.tmp2 | "
              "csvtk cut -t -f qseqid | csvtk del-header 1> {output.contigids} 2>> {log}")
        shell("bin/filtfa -i {input.contigs} -o {output.contigs} -f {output.contigids} -m select -v &> {log}")
        shell("rm -f {output.contigs}.tmp1 {output.contigs}.tmp2")


# TODO: Extracts contigs with no dbhits as "unknown"
rule extract_uncontigs:
    input:
        contigs          = config["res"]+"/{sample}/contigs.fa",
        contigs_annot    = config["res"]+"/{sample}/contigs.annot.tsv"
    output:
        contigflt        = temp(config["res"]+"/{sample}/contigs_un.flt"),
        contigids        = config["res"]+"/{sample}/contigs_un.ids",
        contigs          = config["res"]+"/{sample}/contigs_un.fa"
    log:
        config["logs"]+"/{sample}/extract_uncontigs.log"
    run:
        shell("cat {input.contigs_annot} | "
              "csvtk cut -t -f contig | uniq 1> {output.contigs}.tmp1 2> {log}")
        shell("sed -n \"/^>/s/^>// p\" {input.contigs} | "
              "csvtk add-header -I -t -n qseqid | csvtk mutate -t -f qseqid -n contig -p '^contig=([A-Za-z0-9\.]+)_(.+)' 1> {output.contigs}.tmp2 2>> {log}")
        shell("csvtk join -t -f contig {output.contigs}.tmp1 {output.contigs}.tmp2 | "
              "csvtk cut -t -f qseqid | csvtk del-header 1> {output.contigflt} 2>> {log}")
        shell("bin/filtfa -i {input.contigs} -o {output.contigs} -f {output.contigflt} -m filter -v &> {log}")
        shell("rm -f {output.contigs}.tmp1 {output.contigs}.tmp2")
        shell("sed -n \"/^>/s/^>// p\" {output.contigs} 1> {output.contigids}")
    


rule summary_excel:
    input:
        abund_table         = config["res"]+"/{sample}/abund_table.tsv",
        annot_table         = config["res"]+"/{sample}/contigs.annot.tsv"
    output:
        abund_excel         = config["res"]+"/{sample}/abund_table.xlsx",
        annot_excel         = config["res"]+"/{sample}/contigs.annot.xlsx"
    log:
        config["logs"]+"/{sample}/generate_reports_excel.log"
    params:
        tail            = config["tail"],
        hosttaxid       = lambda wildcards: config["hostgen_taxid_sm"][wildcards.sample]
    run:
        shell(config["R_call"]+" R/print_abund_table.R {input.abund_table} {input.annot_table} {output.abund_excel} {params.tail} {params.hosttaxid}")
        shell(config["R_call"]+" R/print_annot_table.R {input.annot_table} {output.annot_excel} ")

# Simple excel tables for contigs_vi and contigs_un
rule annot_excel:
    input:
        annot_table         = config["res"]+"/{sample}/{contigs}.annot.tsv"
    output:
        annot_excel         = config["res"]+"/{sample}/{contigs,[A-Za-z0-9]+_[A-Za-z0-9]+}.annot.xlsx"
    log:
        config["logs"]+"/{sample}/generate_annot_excel_{contigs}.log"
    run:
        shell(config["R_call"]+" R/print_annot_table.R {input.annot_table} {output.annot_excel} ")


rule krona_graph:
    input:
        abund_table         = config["res"]+"/{sample}/abund_table.tsv"
    output:
        krona_table         = temp(config["res"]+"/{sample}/krona_table.tsv"),
        krona_graph         = config["res"]+"/{sample}/krona_graph.html"
    log:
        config["logs"]+"/{sample}/krona_graph.log"
    params:
        tail                = config["tail"],
        hosttaxid           = lambda wildcards: config["hostgen_taxid_sm"][wildcards.sample]            
    run:
        shell("perl perl/abundtable2krona.pl --tail {params.tail} {input.abund_table} 1> {output.krona_table} 2> {log}")
        shell("ktImportText {output.krona_table} -o {output.krona_graph} -u \"http://krona.sourceforge.net\" 2>> {log}")



#rule stats:
rule stats:
    input:
        r1          = lambda wildcards: config["datain"][wildcards.sample],
        r1_flt      = config["res"]+"/{sample}/trimmed_paired1.fq",
        r1_hgflt    = lambda wildcards: config["res"]+"/{sample}/trimmed_paired1_hostflt.fq" if bool(config["hostgen_sm"][wildcards.sample]) else config["res"]+"/{sample}/trimmed_paired1.fq",
        contigs     = config["res"]+"/{sample}/contigs.fa",
        idx         = config["res"]+"/{sample}/contigs.idxstats",
        orfs       = config["res"]+"/{sample}/ORFs.nt.fa"
    output:
        stats       = config["res"]+"/{sample}/assembly.stats.tsv",
    log:
        config["logs"]+"/{sample}/stats.log"
    run:
        shell("perl perl/assembly_stats2.pl --format col,names "
              "--reads {input.r1} --reads_flt {input.r1_flt} --reads_hgflt {input.r1_hgflt} "
              "--cont {input.contigs} --idx {input.idx} --orfs {input.orfs} 1> {output.stats} 2> {log} ")
        shell("perl perl/assembly_stats.pl {input.contigs} {input.orfs} mean,sum,N50,LN500,Lbp500,LN1000,Lbp1000 col,names 1>> {output.stats} 2>> {log}")



rule qc_plots:
    input:
        r1           = lambda wildcards: config["datain"][wildcards.sample],
        r2           = lambda wildcards: config["datain"][wildcards.sample].replace("R1.","R2."),
        contigs      = config["res"]+"/{sample}/contigs.fa",
        stats        = config["res"]+"/{sample}/assembly.stats.tsv"
    output:
        r1_len       = temp(config["res"]+"/{sample}/r1.readlen"),
        r2_len       = temp(config["res"]+"/{sample}/r2.readlen"),
        cont_len     = temp(config["res"]+"/{sample}/contigs.len"),
        r1_fig       = config["res"]+"/{sample}/qc.read1.pdf",
        r2_fig       = config["res"]+"/{sample}/qc.read2.pdf",
        cont_fig     = config["res"]+"/{sample}/qc.contigs.pdf",
        surv_fig     = config["res"]+"/{sample}/qc.readsurv.pdf"
    params:
        wrkdir       = config["wrkdir"]
    run:
        if input.r1.endswith('.gz'):
            shell("gunzip -c {input.r1} | awk '{{ if(NR%4==2) print length($1) }}' | sort -T {params.wrkdir} -n 1> {output.r1_len}")
        else:
            shell("cat {input.r1} | awk '{{ if(NR%4==2) print length($1) }}' | sort -T {params.wrkdir} -n 1> {output.r1_len}")
        if input.r2.endswith('.gz'):
            shell("gunzip -c {input.r2} | awk '{{ if(NR%4==2) print length($1) }}' | sort -T {params.wrkdir} -n 1> {output.r2_len}")		
        else:
            shell("cat {input.r2} | awk '{{ if(NR%4==2) print length($1) }}' | sort -T {params.wrkdir} -n 1> {output.r2_len}")
        shell("seqkit seq -w0 {input.contigs} | awk '{{ if(NR%2==0) print length($1)  }}' | sort -T {params.wrkdir} -n 1> {output.cont_len}")

        shell(config["R_call"] + " R/qc_plots.R seqlen {output.r1_len} {output.r1_fig}")
        shell(config["R_call"] + " R/qc_plots.R seqlen {output.r2_len} {output.r2_fig}")
        shell(config["R_call"] + " R/qc_plots.R seqlen {output.cont_len} {output.cont_fig}")
        shell(config["R_call"] + " R/qc_plots.R readn  {input.stats} {output.surv_fig}")

def get_required_files(wildcards):
    files = [config["res"]+"/"+wildcards.sample+"/abund_table.tsv",
            config["res"]+"/"+wildcards.sample+"/abund_table.xlsx",
            config["res"]+"/"+wildcards.sample+"/contigs.annot.tsv",        
            config["res"]+"/"+wildcards.sample+"/contigs.annot.xlsx",
            config["res"]+"/"+wildcards.sample+"/krona_graph.html",
			config["res"]+"/"+wildcards.sample+"/taxprofile.txt",
            config["res"]+"/"+wildcards.sample+"/contigs.fa",
            config["res"]+"/"+wildcards.sample+"/contigs",
            config["res"]+"/"+wildcards.sample+"/ORFs.nt.fa",
            config["res"]+"/"+wildcards.sample+"/ORFs.aa.fa",
            config["res"]+"/"+wildcards.sample+"/qc.read1.pdf",
            config["res"]+"/"+wildcards.sample+"/qc.read2.pdf",
            config["res"]+"/"+wildcards.sample+"/qc.contigs.pdf",
            config["res"]+"/"+wildcards.sample+"/assembly.stats.tsv"]
    if config["blastv"]:
        files.append(config["res"]+"/"+wildcards.sample+"/contigs_vi.annot.tsv")
        files.append(config["res"]+"/"+wildcards.sample+"/contigs_vi.annot.xlsx")
    if config["blastu"]:
        files.append(config["res"]+"/"+wildcards.sample+"/contigs_un.annot.tsv")
        files.append(config["res"]+"/"+wildcards.sample+"/contigs_un.annot.xlsx")
    return files
        
rule pack:
    input:
        get_required_files
    output:
        config["res"]+"/{sample}.tar.gz"
    log:
        config["logs"]+"/{sample}/pack.log"
    run:
        # add additional non-required files:
        input        += glob.glob(config["res"]+"/"+wildcards.sample+"/*.html")
        input        += glob.glob(config["res"]+"/"+wildcards.sample+"/*.fa")
        input        += glob.glob(config["res"]+"/"+wildcards.sample+"/*.jpeg")
        input        += glob.glob(config["res"]+"/"+wildcards.sample+"/*.xlsx")
        input        = set(input)    # delete duplicates
                
        shell("mkdir -p results/"+wildcards.sample+".tar")
        shell("cp -r {input} results/"+wildcards.sample+".tar/")
        shell("tar -czf {output} -C results/"+wildcards.sample+".tar .")
        shell("rm -fR results/"+wildcards.sample+".tar")
        
