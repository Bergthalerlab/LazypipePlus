#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
//#include <cstring>

#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <getopt.h>
#include "bioio.h"


/* LAZYPIPE PROJECT: C++ CODE FOR NGS PIPELINE (2022)
 *
 * Retrieve reads for contig or taxid based on contig taxonomy assignments in resdir
 *
 * Author: Ilya Plyusnin, University of Helsinki (2022)
 * Creadit: Lazypipe project, https://doi.org/10.1093/ve/veaa091
 */
using namespace std;

void print_usage(const string name){
	cerr << "USAGE: "<< name <<" -t taxid(s) -c contid [-r resdir -1 read1 -2 read2 -w wrkdir -p prefix -v]\n"
		<<"\n"
		<<"Retrieve reads for contig or taxid based on taxonomy assignments in resdir\n"
		<<"\n\n"
		<<"-t str            : taxon id (or comma-separated list of ids)\n"
		<<"-c str            : short contig id (the id part in contid=id_[\\w+] string).\n"
		<<"-r dir            : directory with pipeline results. Default: results. MUST include files: contig_taxid_score.tsv + readid_contigid.tsv.\n"
		<<"-1 file           : PE reads forward. Default: results/trimmed_paired1_hostflt.fq|trimmed_paired1.fq)\n"
		<<"-2 file           : PE reads reverse. Default: results/trimmed_paired2_hostflt.fq|trimmed_paired2.fq)\n"
		<<"-w dir            : work directory. Default: .\n"
		<<"-p str            : prefix for output read files. Default: taxid or contid\n"
		<<"-v                : verbal mode. Default: false\n"
		<<"output            : $resdir/$prefix_R1.fq and $resdir/$prefix_R2.fq\n"
		<<"\n\n"
		<<"Credit: Lazypipe project, https://doi.org/10.1093/ve/veaa091\n\n"
		<<"\n";
}

inline std::vector<std::string> split(const std::string& str, const char delim)
{
    std::stringstream ss {str};
    std::string item;
    std::vector<std::string> result {};
    while (std::getline(ss, item, delim)){
    	if( isspace(item.back()) ){
		item.erase(item.length()-1,1);
	}
    	result.emplace_back(item);
	}

    return result;
}


// Test file for accessibility/existance
inline bool exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

void readto_map_map(const string file, unsigned int keyi, unsigned int vali, unordered_map<string,unordered_map<string,bool> > &map_map){

	ifstream in(file, ios::in);
	if(!in.is_open()){
		cerr << "\nERROR: failed to open \'"<<file<<"\'\n"; exit(1);
	}
	char buffer[10000];
	vector<string> sp;
	while(in.getline(buffer,10000)){
		if(buffer[0] == '@'){
			continue;
		}
		sp = split(string(buffer),'\t');
		if(  keyi<sp.size() && vali<sp.size()){
			if( map_map.count(sp[keyi]) == 0){
				unordered_map<string,bool> map;
				map_map[sp[keyi]] = map;
			}
			map_map[sp[keyi]][sp[vali]] = true;
		}	
	}
	in.close();
}

int main(int argc, char** argv) {
	
	string prog_name= string(argv[0]);	
	// paramenters
	string taxid		= "";
	string contid		= "";
	string resdir		= "results";
	string wrkdir		= ".";
	string prefix		= "";
	string contid_taxid_file= "";
	string read1		= "";
	string read2		= "";
	string read1_out	= "";
	string read2_out	= "";
	bool verbal	= false;
	
	// PARSING PARAMETERS
	int option= 0; // -t taxid -c contid -r resdir [-1 read1 -2 read2 -w wrkdir -p prefix]\n"
	while ((option = getopt(argc, argv,"t:c:r:w:p:1:2:v")) != -1) {
        switch (option) {
		case 't':
			taxid	= string(optarg);
			break;
		case 'c':
			contid	= string(optarg);
			break;
		case 'r':
			resdir	= string(optarg);
			break;
		case '1':
			read1	= string(optarg);
			break;
		case '2':
			read2	= string(optarg);
			break;
		case 'w':
			wrkdir	= string(optarg);
			break;
		case 'p':
			prefix	= string(optarg);
			break;
		case 'v':
			verbal= true;
			break;	
		case '?':
        		cerr << "option -"<<optopt<<" requires an argument\n\n";
        		exit(1);
             	default:
	     		print_usage(prog_name); 
                	exit(1);
        	}
    	}
	if(taxid=="" && contid==""){ 	cerr<<"ERROR: please specify -t taxid or -c contid\n\n"; print_usage(prog_name); exit(1);}
	if(resdir==""){ 		cerr<<"ERROR: please specify -r resdir\n\n"; print_usage(prog_name); exit(1);}
	if(read1 == ""){
		if( exists(resdir + "/trimmed_paired1_hostflt.fq") &&  exists(resdir + "/trimmed_paired2_hostflt.fq")){
			read1	= resdir + "/trimmed_paired1_hostflt.fq";
			read2	= resdir + "/trimmed_paired2_hostflt.fq";}
		else if( exists(resdir+"/trimmed_paired1.fq") && exists(resdir+"/trimmed_paired1.fq") ){
			read1	= resdir + "/trimmed_paired1.fq";
			read2	= resdir + "/trimmed_paired1.fq";}
		else{
			cerr << "\nERROR: missing read files in "<<resdir<<": trimmed_paired*.fq or trimmed_paired*_hostflt.fq \n\n";
			print_usage(prog_name);
			exit(1);}
	}
	if(read2 == ""){
		read2	= read1;
		size_t found = read2.find("_R1");
		if(found != string::npos){
			read2.replace(found,3,"_R2");
		}
		else{
			cerr << "\nERROR: could not inferr read2 from read1\n\n";
			print_usage(prog_name);
			exit(1);
		}
		if( !exists(read2) ){
			cerr << "\nERROR: missing read2 file: "<< read2 << "\n\n",
			print_usage(prog_name);
			exit(1);
		}
	}
	if( prefix == "" ){
		if( contid!=""){
			prefix = contid; }
		else if(taxid !=""){
			prefix = taxid; }
		else{
			prefix = "tmp"; }
	}
	read1_out= resdir+"/"+prefix+"_R1.fq";
	read2_out= resdir+"/"+prefix+"_R2.fq";

	// contid_taxid_[score] file
	if( exists(resdir+"/contig_taxid.tsv") ){
		contid_taxid_file = resdir+"/contid_taxid.tsv";
	}
	else if(exists(resdir+"/contig_taxid_score.tsv")){
		contid_taxid_file = resdir+"/contig_taxid_score.tsv";
	}
	else{
		cerr << "\nERROR: missing file: "<< (resdir+"/contig_taxid[_score].tsv") << "\n\n";
	}
		

	// Listing taxid(s)
	vector<string> taxid_list = split(taxid,',');

	// Listing contid(s)
	vector<string> contid_list;
	if(contid !=""){
		contid_list.push_back( contid );
	}
	else if( taxid != ""){
	
	
		if(verbal){ cerr << "\t# reading taxid-contid map from "<< contid_taxid_file << "\n";}
		
		unordered_map<string,unordered_map<string,bool> > taxid_contid_map;
		readto_map_map(contid_taxid_file, 1, 0, taxid_contid_map);

		unordered_map<string,bool> contid_map;
		for(unsigned int i=0; i<taxid_list.size(); i++){
			if( taxid_contid_map.count(taxid_list[i]) > 0){
				contid_map.insert( taxid_contid_map[taxid_list[i]].begin(), taxid_contid_map[taxid_list[i]].end() );
			}
		}
		if( contid_map.size() == 0){
			cerr << "\tWARNING: no contigs found for taxid(s) "+taxid+"\n";
			exit(0);
		}
		for(auto it = contid_map.begin(); it != contid_map.end(); ++it ){
			contid_list.push_back( it->first );
		}
	}
	if(verbal){ cerr << "\t# contids found: "<<contid_list.size()<<"\n";}
	
	
	// List reads for contids
    string readid_contid_file = resdir+"/readid_contigid.tsv";
	if(verbal){ cerr << "\t# reading "<< readid_contid_file <<"\n";}
	unordered_map<string,unordered_map<string,bool> > contid_readid_map_tmp;
	unordered_map<string,unordered_map<string,bool> > contid_readid_map;
	readto_map_map(readid_contid_file, 1, 0, contid_readid_map_tmp);

	// converting contids to short format
	vector<string> sp;
	for(auto it = contid_readid_map_tmp.begin(); it != contid_readid_map_tmp.end(); ++it){
		string contid = it->first;
		sp	= split(contid,'_');
		contid	= sp[0];
		sp	= split(contid,'=');
		if(sp.size() > 1){
			contid = sp[1]; // assuming contig id of form "contid|contig=id_.*"
		}
		unordered_map<string,bool> map_tmp	= (it->second);
		contid_readid_map[contid] 		= map_tmp;
		//DEBUG
		//cout << "key: "<< it->first <<"\n";
		//cout << "contid: "<< contid<<"\n";
	}
	
	unordered_map<string,bool> readid_map;
	
	for(unsigned int i=0; i<contid_list.size(); i++){
		if(contid_readid_map.count(contid_list[i]) == 0){
			fprintf(stderr,"\t# scipping contig with no reads: %s\n",contid_list[i].c_str());
			continue;
		}
		
		unordered_map<string,bool> map_tmp	= contid_readid_map[contid_list[i]];
		
		for(auto it = map_tmp.begin(); it != map_tmp.end(); ++it)
			readid_map[ it->first ] = true;
	}
	if(verbal){ cerr << "\t# readids found: "<<readid_map.size()<<"\n";}


	// Filtering reads
	vector<string> read_files;
	read_files.push_back(read1);
	read_files.push_back(read2);
	vector<string> read_files_out;
	read_files_out.push_back(read1_out);
	read_files_out.push_back(read2_out);
	
	for(unsigned int k=0; k<read_files.size(); k++){
		if(verbal){ cerr << "\t# processing "<<read_files[k]<<"\n";}
	
		std::ifstream fastq {read_files[k], std::ios::binary};
		std::ofstream fout (read_files_out[k], std::ofstream::out);
		unsigned int seqn	=  bioio::count_fastq_records(fastq);
		unsigned int seqn_corr	= 0;
		unsigned int sel_num	= 0;
		vector<string>	sp;
		
		for(unsigned int i=0; i<seqn; i++,seqn_corr++){
			auto seq = bioio::detail::read_fastq_record<string,string,string>(fastq);
			if(seq.name.length()== 0){
			// count_fastq_records counts "^@" lines in fastq, which can be larger than rec num due to matches in quality string
				break;
			}
			string name = bioio::detail::split(seq.name,' ')[0];
			//cerr << "name:"<<name<<"\n"; exit(1);
			if( name[0]=='@' ){
				name = name.substr(1);
			}
			if( name.substr(name.length()-2,name.length()) == "/1" || name.substr(name.length()-2,name.length()) == "/2"){
				name = name.substr(0,name.length()-2);
			}
			//cerr<< "seqname:\""<<name<<"\"\n";exit(1);
			
			if( readid_map.count(name)>0 ){
				fout << seq.name << "\n";
				fout << seq.seq << "\n";
				fout << "+\n";
				fout << seq.qual<< "\n";
				sel_num++;
			}
		}
		fastq.close();
		fout.close();
		if(verbal){
			fprintf(stderr,"\t# selected %u/%u (%2.2f%%) reads to %s\n",
				sel_num,seqn_corr,((sel_num+0.1)/seqn_corr)*100.0,read_files_out[k].c_str());
		}
	}
		
	
}




