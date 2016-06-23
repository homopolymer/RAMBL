#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <list>
#include <tuple>
#include <map>
using namespace std;



long double thresh_identity = 95.;
long double thresh_evalue = 1e-10;
string blast_outfile = "";

void help_msg(){
    cerr << "usage: blastout2abundance [-I identity_thresh] [-E evalue_thresh] blast_output" << endl
         << "" << endl
         << "Options" << endl
         << "  -I    Threshold for sequence identity [95.0]" << endl
         << "  -E    Threshold for E-value [1e-10]" << endl
         << "  -h    Print this message" << endl
         << "" << endl;
}


void parse_cmd_line(int argc, char** argv, string& blast_outfile, long double& thresh_identity, long double& thresh_evalue){
    for (int i=0; i<argc; ++i)
    {
        string op = string(argv[i]);
        if (op[0]=='-')
        {
            if (op=="-h" or op=="--help" or op=="-h")
            {
                help_msg();
                exit(0);
            }else if (op=="-I")
            {
                thresh_identity = stold(string(argv[++i]));
            }else if (op=="-E")
            {
                thresh_evalue = stold(string(argv[++i]));
            }
        }else
        {
            blast_outfile = string(argv[i]);
        }
    }
}


class ContigHit{
public:
    ContigHit(){id="";len=0;evalue=0;}
    ~ContigHit(){}
public:
    string id;
    long double len;
    long double evalue;
};

typedef long double Evalue;
typedef tuple<list<ContigHit>,Evalue> ContigHits;
map<string,ContigHits> ReadInfo;
map<string,long double> ContigRawAbundance;
map<tuple<string,string>,int> ReadContigCount;


void read_blast_outfile(string blast_outfile){
    ifstream data(blast_outfile);
    
    string line;
    while (getline(data,line)){
        stringstream line_stream(line);
        string cell;

        int i = 0;

        string read_segment_id;
        long double identity;
        long double evalue;
        long double align_len;
        string contig_id;
        long double read_start;
        long double read_end;
        long double contig_start;
        long double contig_end;
        long double len;
        long double read_len;
        while (getline(line_stream,cell,',')){
            if (i==0) {read_segment_id = cell;}
            else if (i==1) {contig_id = cell;}
            else if (i==2) {identity = stold(cell);}
            else if (i==3) {align_len = stold(cell);}
            else if (i==4) {read_start = stold(cell);}
            else if (i==5) {read_end = stold(cell);}
            else if (i==6) {contig_start = stold(cell);}
            else if (i==7) {contig_end = stold(cell);}
            else if (i==8) {evalue = stold(cell);}
            else if (i==9) {read_len = stold(cell);}
            i++;
        }


        long double perc_identity = identity*100/align_len;

        if (perc_identity<thresh_identity || evalue>thresh_evalue){
            continue;
        }

        //len = (read_start>read_end)?(read_start-read_end+1):(read_end-read_start+1);
        //if (len<0.7*read_len) {continue;}

        // get read name, remove /1 and /2
        string read_id = read_segment_id;
        bool left_segment = true;
        if (read_segment_id.substr(read_segment_id.length()-2,2)=="/1" ||
            read_segment_id.substr(read_segment_id.length()-2,2)==".1" ){
            read_id = read_segment_id.substr(0,read_segment_id.length()-2);
        }else if (read_segment_id.substr(read_segment_id.length()-2,2)=="/2" ||
            read_segment_id.substr(read_segment_id.length()-2,2)==".2"){
            read_id = read_segment_id.substr(0,read_segment_id.length()-2);
            left_segment = false;
        }

        // collect read information
        ContigHit contig_hit;
        contig_hit.id = contig_id;
        contig_hit.evalue = evalue;
        contig_hit.len = ((contig_start>contig_end)?contig_start-contig_end:contig_end-contig_start)+1;

        // check duplicate
        auto key = make_tuple(read_segment_id,contig_id);
        if (ReadContigCount.count(key)>0) {continue;}
        ReadContigCount[key] = 1;

        // update read information
        auto found = ReadInfo.find(read_id);
        if (found==ReadInfo.end()){
            list<ContigHit> contig_hits {contig_hit};
            ReadInfo[read_id] = make_tuple(contig_hits,evalue);
        }else{
            long double evalue0 = get<1>(ReadInfo[read_id]);
            if (evalue0>evalue){
                list<ContigHit> contig_hits {contig_hit};
                ReadInfo[read_id] = make_tuple(contig_hits,evalue);
            }else if (evalue0==evalue){
                get<0>(ReadInfo[read_id]).emplace_back(contig_hit);
            }
        }
  
    }
}

void compute_contig_raw_abundance(){
    // iterate over reads
    auto iter1 = ReadInfo.begin();
    for (;iter1!=ReadInfo.end();iter1++){
        long double mean_len = 0;
        int n = 0, mxn = 0;
        map<string,int> contig_counts;
        map<string,int> contig_counts2;
        auto iter2 = get<0>(iter1->second).begin();
        for (;iter2!=get<0>(iter1->second).end();iter2++){
            n += 1;
            mean_len += iter2->len;
            auto found = contig_counts.find(iter2->id);
            if (found==contig_counts.end()) {contig_counts[iter2->id]=1;}
            else {contig_counts[iter2->id]+=1;}
        }
        mean_len /= n; 

        // 
        auto iter4 = contig_counts.begin();
        for (;iter4!=contig_counts.end();iter4++){
            if (mxn<iter4->second) {mxn=iter4->second;}
        }
        auto iter5 = contig_counts.begin();
        for (;iter5!=contig_counts.end();iter5++){
            if (iter5->second==mxn) {contig_counts2[iter5->first] = iter5->second;}
        }

        //
        n = contig_counts2.size();
        auto iter3 = contig_counts2.begin();
        for (;iter3!=contig_counts2.end();iter3++){
            auto found = ContigRawAbundance.find(iter3->first);
            if (found==ContigRawAbundance.end()){
                ContigRawAbundance[iter3->first] = iter3->second/(n+0.);
            }else{
                ContigRawAbundance[iter3->first] += iter3->second/(n+0.);
            }
        }
    }
}

int main(int argc, char** argv){
    if (argc==1) {
        help_msg();
        exit(0);
    }
    parse_cmd_line(--argc,++argv,blast_outfile,thresh_identity,thresh_evalue);

    // load blast results
    read_blast_outfile(blast_outfile);

    // compute raw abundance
    compute_contig_raw_abundance();

    // output results
    auto iter=ContigRawAbundance.begin();
    for (;iter!=ContigRawAbundance.end();iter++){
        cout << iter->first << "\t" << iter->second << endl;
    }

}
