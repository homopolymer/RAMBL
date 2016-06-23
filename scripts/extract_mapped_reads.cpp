#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <map>
using namespace std;


class Read{
public:
    Read(){id="";first_segment="";second_segment="";}
    ~Read(){}

public:
    string id;
    string first_segment;
    string second_segment;
};

map<string,Read> ReadPool;

void help_msg(){
    cerr << "usage: extract_mapped_reads bam_file prefix" << endl
         << "" << endl;
}


void load_reads(string samfile){
    ifstream data(samfile);

    string line;
    while (getline(data,line)){
        stringstream line_stream(line);
        string cell;
        int i = 0;
        string read_id;
        string read_seq;
        int read_flag;
        while (getline(line_stream,cell,'\t')){
            if (i==0) {read_id = cell;}
            else if (i==1) {read_flag = stoi(cell);}
            else if (i==9) {read_seq = cell;}

            if (read_id.substr(read_id.length()-2,2)=="/1" ||
                read_id.substr(read_id.length()-2,2)=="/2"){
                read_id = read_id.substr(0,read_id.length()-2);
            }

            i ++;
        }

        // collect read information
        auto found = ReadPool.find(read_id);
        if (found==ReadPool.end()){
            Read read;
            read.id = read_id;
            if (read_flag&0x40) {read.first_segment = read_seq;}
            else if (read_flag&0x80) {read.second_segment = read_seq;}
            else {read.first_segment = read_seq;}
            ReadPool[read_id] = read;
        }else{
            if (read_flag&0x40) {found->second.first_segment = read_seq;}
            else if (read_flag&0x80) {found->second.second_segment = read_seq;}
            else {found->second.first_segment = read_seq;}
        }
    }
}


int main(int argc, char** argv){
    if (argc==1) {help_msg(); exit(0);}

    string bamfile = string(argv[1]);
    string prefix = string(argv[2]);

    // convert bam to sam
    string samfile = prefix+"_mapped_reads.tsv";
    string cmd = "samtools view -F4 "+bamfile+" > "+samfile+" 2>/dev/null";
    system(cmd.c_str());

    // load reads
    load_reads(samfile);

    // write to disk files
    ofstream pair1_out(prefix+".1.fa");
    ofstream pair2_out(prefix+".2.fa");
    ofstream single_out(prefix+".single.fa");

    auto iter = ReadPool.begin();
    for (; iter!=ReadPool.end(); iter++){
        auto read = iter->second;
        if (read.first_segment.length()>0 && read.second_segment.length()>0){
            pair1_out << (">"+read.id+"/1") << endl
                      << read.first_segment << endl;
            pair2_out << (">"+read.id+"/2") << endl
                      << read.second_segment << endl;
        }else if (read.first_segment.length()>0){
            single_out << (">"+read.id) << endl
                       << read.first_segment << endl;
        }else{
            single_out << (">"+read.id) << endl
                       << read.second_segment << endl;
        }
    }

    pair1_out.close();
    pair2_out.close();
    single_out.close();

    // remove samfile
    cmd = "rm -f "+samfile;
    system(cmd.c_str());
}

