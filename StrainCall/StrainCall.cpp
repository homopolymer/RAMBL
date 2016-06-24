// StrainCall 
// ----------
// It is used to call strains for a marker gene.
//
//
// Example
// -------
//
//
// Reference
// ---------
// Zeng F and Chen T. Full-length 16S rRNA gene reconstruction at strain level.
//
// Last Changed
// ------------
// July 30, 2015    Feng Zeng    Create it.
// October 3, 2015  Feng Zeng    Add a count filter to assembly
//
//
#include "PartialOrderGraph.hpp"
#include <sys/stat.h>
#include <cstdlib>
#include <ctime>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>
#include <algorithm>
using namespace std;

void help()
{
    cerr << "StrainCall marker_gene read_mapping" << endl;
    cerr << "           [-r gn:p0-p1] [-w window_size]" << endl;
    cerr << "           [-e error_rate] [-q map_qual]" << endl;
    cerr << "" << endl;
    cerr << "Options" << endl;
    cerr << "-r,--roi           region of interesting, gn is gene name," << endl;
    cerr << "                   p0 is starting position, p1 is ending position (inclusive)" << endl;
    cerr << "-w,--window        the size of scanning window [500]" << endl;
    cerr << "-o,--overlap       the size of window-window overlap [100]" << endl;
    cerr << "-e,--error-rate    sequencing error rate [0.01]" << endl;
    cerr << "-D,--max-depth     downsample data to the specified depth [800]" << endl;
    cerr << "-q,--map-qual      only include reads with mapping quality >= INT [3]" << endl;
    cerr << "-I,--max-ins       only include reads with insertions <= INT [10]" << endl;
    cerr << "-l,--read-len      only include reads with length >=INT [80]" << endl;
    cerr << "-t,--tau           only include strains with abundance level >=FLT [0.02]" << endl;
    cerr << "-d,--diff-rate     only include strains with difference rate >=FLT [0.01]" << endl;
    cerr << "-G,--plot-graph    print graph" << endl;
    cerr << "-h,--help          print this message" << endl;
    cerr << "" << endl;
}


struct sc_parameter
{
    sc_parameter()
    {
        gene_file = "";
        mapping_file = "";
        roi = "";
        window_size = 500;
        overlap_size = 100;
        error_rate = 0.01;
        mapping_qual = 3;
        max_ins = 10;
        read_len = 80;
        print_help = false;
        d0 = 0;
        d1 = 0;
        tau = 0.02;
        diff_rate = 0.01;
        max_depth = 800;
        plot_graph = false;
    }
    string gene_file;
    string mapping_file;
    string roi;
    int    window_size;
    int    overlap_size;
    float  error_rate;
    int    mapping_qual;
    int    max_ins;
    int    read_len;
    bool   print_help;
    float  tau;
    float  diff_rate;
    int    max_depth;
    int    d0;
    int    d1;
    bool   plot_graph;
};


void sc_parse_cmd_line(int argc, char** argv, sc_parameter* parameters)
{
    for (int i=0,p=0; i<argc; ++i)
    {
        string op = string(argv[i]);
        if (op[0]=='-')
        {
            if (op=="-h" or op=="--help" or op=="-h")
            {
                parameters->print_help = true;
            }else if (op=="-r" or op=="--roi" or op=="-roi")
            {
                parameters->roi = string(argv[++i]);
            }else if (op=="-w" or op=="--window" or op=="-window")
            {
                parameters->window_size = stoi(string(argv[++i]));
            }else if (op=="-e" or op=="--error-rate" or op=="-error-rate")
            {
                parameters->error_rate = stof(string(argv[++i]));
            }else if (op=="-q" or op=="--map-qual" or op=="-map-qual")
            {
                parameters->mapping_qual = stoi(string(argv[++i]));
            }else if (op=="-o" or op=="--overlap" or op=="-overlap")
            {
                parameters->overlap_size = stoi(string(argv[++i]));
            }else if (op=="-l" or op=="--read-len" or op=="-read-len")
            {
                parameters->read_len = stoi(string(argv[++i]));
            }else if (op=="-t" or op=="--tau" or op=="-tau")
            {
                parameters->tau = stof(string(argv[++i]));
            }else if (op=="-d" or op=="--diff-rate" or op=="-diff-rate")
            {
                parameters->diff_rate = stof(string(argv[++i]));
            }else if (op=="-D" or op=="--max-depth" or op=="-max-depth")
            {
                parameters->max_depth = stoi(string(argv[++i]));
            }else if (op=="-I" or op=="--max-ins" or op=="-max-ins")
            {
                parameters->max_ins = stoi(string(argv[++i]));
            }else if (op=="-G" or op=="--plot-graph" or op=="-plot-graph")
            {
                parameters->plot_graph = true;
            }
        }else
        {
            if (p==0)
            {
                parameters->gene_file = string(argv[i]);
            }else
            {
                parameters->mapping_file = string(argv[i]);
            }
            p++;
        }
    }
}


void load_gene_seq(string& gene_file, string& gene_roi, GenomeSeq& gene_seq)
{
    time_t t = time(0);
    srand(t);
    int r = rand();
    string tmp_file = gene_roi+"_"+to_string(t)+"_"+to_string(r);
    
    string cmd;
    
    // extract seq to temporary file
    cmd = "samtools faidx " + gene_file + " " + gene_roi + " 2>/dev/null 1>" + tmp_file;
    system(cmd.c_str());

    // read temporary file
    gene_seq = "";
    ifstream in(tmp_file);
    for (string line; getline(in,line);)
    {
        if (line[0]!='>')
        {
            gene_seq += line;
        }
    }
    in.close();

    // remove temporary file
    cmd  = "rm " + tmp_file;
    system(cmd.c_str());
}


string gene_roi_name(string& gene_roi)
{
    size_t x = gene_roi.find_first_of(':');
    return gene_roi.substr(0,x);
}

int gene_roi_start_pos(string& gene_roi)
{
    string pos="";
    size_t x = gene_roi.find_first_of(':');
    auto it = next(gene_roi.begin(),x+1);
    for (;it!=gene_roi.end();++it)
    {
        if (*it=='-')
        {
            break;
        }
        pos += string(1,*it);
    }
    return stoi(pos);
}

int gene_roi_end_pos(string& gene_roi)
{
    string pos="";
    size_t x = gene_roi.find_first_of('-');
    auto it = next(gene_roi.begin(),x+1);
    for (;it!=gene_roi.end();++it)
    {
        pos += string(1,*it);
    }
    return stoi(pos);
}

string gene_name(string& gene_file)
{
    string gn;
    string cmd;
    struct stat buffer;
    if (!(stat(gene_file.c_str(),&buffer)==0))
    {
        cmd = "samtools faidx " + gene_file + " 2>/dev/null";
        system(cmd.c_str());
    }

    ifstream in(gene_file+".fai");
    for (string line; getline(in,line);)
    {
        stringstream ss(line);
        string f1,f2,f3,f4,f5;
        ss >> f1 >> f2 >> f3 >> f4 >> f5;
        if (!f1.empty())
        {
            gn = f1;
        }
    }
    in.close();

    return gn;
}

int gene_length(string& gene_file, string gene_name)
{
    int len=0;
    string cmd;
    struct stat buffer;
    if (!(stat(gene_file.c_str(),&buffer)==0))
    {
        cmd = "samtools faidx " + gene_file + " 2>/dev/null";
        system(cmd.c_str());
    }

    ifstream in(gene_file+".fai");
    for (string line; getline(in,line);)
    {
        stringstream ss(line);
        string f1,f2,f3,f4,f5;
        ss >> f1 >> f2 >> f3 >> f4 >> f5;
        if (f1==gene_name)
        {
            len = stoi(f2);
        }
    }
    in.close();
    return len;
}


int read_align_end_pos(int p0, vector<CigarRecord> cigars)
{
    for (auto it=cigars.begin(); it!=cigars.end(); ++it)
    {
        CigarOp op = get<0>(*it);
        CigarOpLen opl = get<1>(*it);
        if (op=="M" or op=="D")
        {
            p0 += opl;
        }
    }

    return (p0-1);
}

void crop_read_within_window(int wp0,int wp1,
                             string& seq,string& qual,
                             vector<CigarRecord>& cigars,
                             int rp0,int rp1,
                             string& crop_seq,string& crop_qual,string& crop_cigar)
{
    int i=0,j=0,ki=0,kj=0,opl; 
    string op;
    vector<CigarRecord> crop_cigars;

    auto it = cigars.begin();
    op = get<0>(*it);
    opl = get<1>(*it);
    if (op=="S")
    {
        i += opl;
        ++it;
    }

    op = get<0>(*it);
    opl = get<1>(*it);
    if (rp0<wp0 and rp0<wp1)
    {
        while (rp0<wp0 and rp0<wp1)
        {
            ki = 0;
            op = get<0>(*it);
            opl = get<1>(*it);
            if (op=="M")
            {
                for (; ki<opl; ++ki,++i,++rp0)
                {
                    if (rp0==wp0)
                    {
                        break;
                    }
                }
            }else if (op=="D")
            {
                for (; ki<opl; ++ki,++rp0)
                {
                    if (rp0==wp0)
                    {
                        break;
                    }
                }
            }else if (op=="I")
            {
                i += opl;
            }
            ++it;
        }
    }else
    {
        ++it;
    }
    if (ki<opl)
    {
        crop_cigars.push_back(CigarRecord(op,opl-ki));
    }
    for (;it!=cigars.end();++it)
    {
        crop_cigars.push_back(*it);
    }

    auto rit = cigars.rbegin();
    op = get<0>(*rit);
    opl = get<1>(*rit);
    if (op=="S")
    {
        j+=opl;
        ++rit;
        crop_cigars.pop_back();
    }
    op = get<0>(*rit);
    opl = get<1>(*rit);
    while (rp1>wp1 and rp1>wp0)
    {
        kj = 0;
        op = get<0>(*rit);
        opl = get<1>(*rit);
        if (op=="M")
        {
            for (; kj<opl;++kj,++j,--rp1)
            {
                if (rp1==wp1)
                {
                    break;
                }
            }
        }else if (op=="D")
        {
            for (; kj<opl;++kj,--rp1)
            {
                if (rp1==wp1)
                {
                    break;
                }
            }
        }else if (op=="I")
        {
            j += opl;
        }
        ++rit;
        if (kj==opl or op=="I")
        {
            crop_cigars.pop_back();
        }else
        {
            get<1>(*crop_cigars.rbegin()) -= kj;
        }
    }

    // crop
    crop_seq = seq.substr(i,seq.length()-i-j);
    crop_qual = qual.substr(j,qual.length()-i-j);
    crop_cigar = "";
    for (auto cit=crop_cigars.begin(); cit!=crop_cigars.end(); ++cit)
    {
        op = get<0>(*cit);
        opl = get<1>(*cit);
        crop_cigar += to_string(opl)+op;
    }
}

int number_of_ambiguous_base(string& read)
{
    int n = 0;
    for_each(read.begin(),read.end(),[&n](char x)
            {
                if (x=='N' or x=='n')
                    n += 1;
            });
    return n;
}

int max_insert_size(string& cigar)
{
    vector<CigarRecord> cigars;
    parse_cigar(cigar,cigars);
    
    int ins = 0;
    for_each(cigars.begin(),cigars.end(),[&ins](CigarRecord c)
            {
                if (get<1>(c)>ins and get<0>(c)=="I")
                    ins = get<1>(c);
            });

    return ins;
}

double read_error_frac(string& gene_seq, int g0, string& read_seq, string& cigar)
{
    double erf = 0;
    int i=g0,j=0,len=0;

    vector<CigarRecord> cigars;
    parse_cigar(cigar,cigars);

    string op;
    int opl;

    for (auto c=cigars.begin(); c!=cigars.end(); ++c)
    {
        op = get<0>(*c);
        opl = get<1>(*c);
 
        if (op=="I") j+=opl;
        else if (op=="D") 
        {
            i+=opl;
            len += opl;
        }
        else
        {
            len += opl;
            for (int k=0; k<opl; k++,i++,j++)
            {
                if (gene_seq[i]!=read_seq[j])
                {
                    erf += 1;
                }
            }
        }
    }

    return erf/len;
}

void load_mapping_reads(string& gene_seq, string& mapping_file, int mq, int rl, int max_ins, int max_depth,
                        string& gene_roi, vector<AlignRead>& mapping_reads, vector<ReadAux>& aux_reads,
                        ReadPairs& read_pairs)
{
    time_t t = time(0);
    srand(t);
    int r = rand();
    string tmp_file = gene_roi+"_"+to_string(t) + "_" + to_string(r);

    string cmd;
    DoubleL rho;
    std::random_device rd;
    std::mt19937 gen(1234);
    std::uniform_real_distribution<> dicer(0, 1);
    
    // extract reads to temporary file
    cmd = "samtools view " + mapping_file + " -q " + to_string(mq) + " -F 1804 " + gene_roi + " 2>/dev/null 1>" + tmp_file;
    system(cmd.c_str());

    // starting position at gene
    int p0 = gene_roi_start_pos(gene_roi);
    int p1 = gene_roi_end_pos(gene_roi);

    ifstream in(tmp_file);
    // sequencing depth
    int depth = 0;
    for (string line;getline(in,line);)
    {
        stringstream ss(line);
        string f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11;
        ss >> f1 >> f2 >> f3 >> f4 >> f5 >> f6 >> f7 >> f8 >> f9 >> f10 >> f11;

        vector<CigarRecord> cigars;
        parse_cigar(f6,cigars);
        int len = 0;
        int r0 = stoi(f4);
        for_each(cigars.begin(),cigars.end(),[&len](CigarRecord cr)
                 {
                     if (get<0>(cr)=="M" || get<0>(cr)=="D")
                         len += get<1>(cr);
                 });
        int r1 = r0+len-1;
        if (p0<=r0 and p1>r1) depth += r1-r0+1;
        else if (p0<=r0 and p1<=r1) depth += p1-r0+1;
        else if (p0>r0 and p1<=r1) depth += p1-p0+1;
        else if (p0>r0 and p1>r1) depth += r1-p0+1;
    }
    in.close();
    depth /= p1-p0+1;
    rho = min(1.0,max_depth/(depth+0.));

    // parse temporary file
    map<AlignRead,vector<string>> tmp_duplicates;
    int id=0;
    in.open(tmp_file);
    for (string line;getline(in,line);++id)
    {
        stringstream ss(line);
        string f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11;
        string rn;
        int flag;
        ss >> f1 >> f2 >> f3 >> f4 >> f5 >> f6 >> f7 >> f8 >> f9 >> f10 >> f11;

        if (f10.length()<rl)
        {
            continue;
        }
        if (number_of_ambiguous_base(f10)>0)
        {
            continue;
        }

        // read name
        rn = f1;
        flag = stoi(f2);
        if ((flag&65)==65)
        {
            rn += "/1";

        }else if ((flag&129)==129)
        {
            rn += "/2";
        }

        vector<CigarRecord> cigars;
        parse_cigar(f6,cigars);

        int read_p0 = stoi(f4);
        int read_p1 = read_align_end_pos(read_p0,cigars);

        int relative_pos = read_p0-p0;
        if (relative_pos<0)
        {
            relative_pos = 0;
        }

        string seq,qual,cigar;
        crop_read_within_window(p0,p1,f10,f11,cigars,read_p0,read_p1,seq,qual,cigar);   

        // error rate
        double erf = read_error_frac(gene_seq,relative_pos,seq,cigar);

        // max insertion size
        int maxins = max_insert_size(cigar);

        // save to temporary pool
        if (seq.length()>rl and maxins<max_ins)
        {
            // sampling
            if (dicer(gen)>rho)
            {
                continue;
            }

            AlignRead ar(relative_pos,cigar,seq,"",1);
            auto it = tmp_duplicates.find(ar);
            if (it==tmp_duplicates.end())
            {
                tmp_duplicates[ar] = vector<string>(1,rn);
            }else
            {
                it->second.push_back(rn);
            }
        }
    }
    in.close();

    // move to result
    map<string,int> tmp_read_uids;
    id = 0;
    for (auto it=tmp_duplicates.begin(); it!=tmp_duplicates.end(); ++it,++id)
    {
        AlignRead ar = it->first;
        get<4>(ar) = it->second.size();
        mapping_reads.push_back(ar);

        // save to aux
        for (auto tt=it->second.begin(); tt!=it->second.end(); ++tt)
        {
            aux_reads.push_back(ReadAux(id,*tt));
        }

        // save to read uid
        for (auto tt=it->second.begin(); tt!=it->second.end(); ++tt)
        {
            tmp_read_uids[*tt] = id;
        }
    }

    // save to read pairs
    for (auto it=tmp_read_uids.begin(); it!=tmp_read_uids.end(); ++it)
    {
        string rn1 = it->first;
        int uid = it->second;
        string rn2;
        bool paired = false;

        if (rn1.substr(rn1.length()-2,2)=="/1")
        {
            paired = true;
            rn2 = rn1.substr(0,rn1.length()-2)+"/2";
        }else if (rn1.substr(rn1.length()-2,2)=="/2")
        {
            paired = true;
            rn2 = rn1.substr(0,rn1.length()-2)+"/1";
        }
   
        auto tt = read_pairs.find(uid);
        if (paired)
        {
            auto iz = tmp_read_uids.find(rn2);
            if (iz==tmp_read_uids.end())
            {
                if (tt==read_pairs.end()) read_pairs[uid] = vector<int>(1,-1);
                else read_pairs[uid].push_back(-1);
            }else
            {
                if (tt==read_pairs.end()) read_pairs[uid] = vector<int>(1,iz->second);
                else read_pairs[uid].push_back(iz->second);
            }
        }else
        {
            if (tt==read_pairs.end()) read_pairs[uid] = vector<int>(1,-1);
            else read_pairs[uid].push_back(-1);
        }
    }

    // remove temporary file
    cmd = "rm " + tmp_file;
    system(cmd.c_str());
}


void window_adjust(string& mapping_file, int mq, string gn, int p0, int p1, int z, int L, int& d0, int& d1)
{
    int P,Q;
    string cmd;

    if (p0-z<1)
    {
        z = p0-1;
    } 
    P = p0-z;

    Q = p1+z;
    if (Q>L)
    {
        Q = L;
    }
    
    // extract pileup to temporary file
    size_t t = time(0);
    srand(t);
    int r = rand();
    string tmp_file = gn+":"+to_string(p0)+"-"+to_string(p1)+"_"+to_string(t) + "_" + to_string(r);
    
    cmd = "samtools mpileup -q " + to_string(mq) + " -Q0 " + " -A "
          + " -r " + gn + ":" + to_string(P) + "-" + to_string(Q) + " " + mapping_file
          + " 2>/dev/null 1>" + tmp_file;
    system(cmd.c_str());

    // parse temporary file
    typedef tuple<bool,bool> mapping_state; // <has_insert,has_delete>
    map<int,mapping_state> mapping_info;

    ifstream in(tmp_file);
    for (string line; getline(in,line);)
    {
        size_t pp;
        mapping_state ms(false,false);

        stringstream ss(line);
        string f1,f2,f3,f4,f5,f6;
        ss >> f1 >> f2 >> f3 >> f4 >> f5 >> f6;

        // check insert
        pp = f5.find_first_of('+');
        if (pp!=f5.npos)
        {
            get<0>(ms) = true;
        }

        // check delete
        pp = f5.find_first_of('-');
        if (pp!=f5.npos)
        {
            get<1>(ms) = true;
        }
        pp = f5.find_first_of('*');
        if (pp!=f5.npos)
        {
            get<1>(ms) = true;
        }
        
        // save mapping state
        mapping_info[stoi(f2)] = ms;
    }
    in.close();

    // find a good start position
    auto it0 = mapping_info.find(p0);
    if (it0 == mapping_info.end())
    {
        P = mapping_info.begin()->first;
    }else
    {
        P = p0;
        while (get<0>(it0->second) or get<1>(it0->second))
        {
            --it0;
            if (it0==mapping_info.end())
            {
                break;
            }
            --P;
        }
    }

    // find a good end position
    auto it1 = mapping_info.find(p1);
    if (it1 == mapping_info.end())
    {
        Q = mapping_info.rbegin()->first;
    }else
    {
        Q = p1;
        while (get<0>(it1->second) or get<1>(it1->second))
        {
            ++it1;
            if (it1==mapping_info.end())
            {
                break;
            }
            ++Q;
        }
    }
    
    d0 = p0-P;
    d1 = Q-p1;

    // remove temporary file
    cmd = "rm " + tmp_file;
    system(cmd.c_str());
}

struct sc_window
{
    sc_window(string gn_,int p0_,int p1_)
    {
        gn = gn_;
        p0 = p0_;
        p1 = p1_;
    }
    string gn;
    int p0;
    int p1;
};

void make_scan_window(sc_parameter* parameters, vector<sc_window>& windows)
{
    int p0,p1,d0=0,d1=0;
    int z=50,l,L,LL;
    string gene_file = parameters->gene_file;
    string mapping_file = parameters->mapping_file;
    string gene_roi = parameters->roi;
    int mq = parameters->mapping_qual;
    int wsize = parameters->window_size;
    int osize = parameters->overlap_size;
    string gn;

    if (gene_roi.empty())
    {
        gn = gene_name(gene_file);
        l = 1;
        L = gene_length(gene_file,gn);
        LL = gene_length(gene_file,gn);
    }else
    {
        gn = gene_roi_name(gene_roi);
        l = gene_roi_start_pos(gene_roi);
        L = gene_roi_end_pos(gene_roi);
        LL = gene_length(gene_file,gn);
    }
    
    set<int> visited;
    for (p0=l,p1=l; p1<L; p0+=wsize-osize)
    {
        p1 = p0+wsize-1;
        if (p1>L)
        {
            p1 = L;
        }
        window_adjust(mapping_file,mq,gn,p0,p1,z,LL,d0,d1);
        if (p0==l)
        {
            parameters->d0 = d0;
        }

        if (visited.count(p1+d1)>0)
        {
            continue;
        }

        windows.push_back(sc_window(gn,p0-d0,p1+d1));
        visited.emplace(p1+d1);
    }
    parameters->d1 = d1;

}

void strain_identity(vector<Strain>& strains)
{
    for (int i=0; i<strains.size(); i++)
    {
        string si = strains[i].strain_seq();
        for (int j=i+1; j<strains.size(); j++)
        {
            string sj = strains[j].strain_seq();
            
            double iden=0, len=0;
            for (int k=0; k<si.length(); k++)
            {
                if ((si[k]=='-' or si[k]=='=') and 
                    (sj[k]=='-' or sj[k]=='='))
                {
                    continue;
                }
                if (si[k]==sj[k])
                {
                    iden += 1;
                }
                len += 1;
            }

            cout << "strain pair (" << i << "," << j << "), identity=" << iden/len << endl;
        }
    }
}

struct LocalResult
{
    LocalResult(string _gn, int _p0, int _p1, vector<Strain>& _strains)
    {
        gn = _gn;
        p0 = _p0;
        p1 = _p1;
        strains = _strains;
    }

    string gn;
    int p0;
    int p1;
    vector<Strain> strains;
};

bool read_is_paired(string& rn)
{
    if (rn.substr(rn.length()-2,2)=="/1" or
        rn.substr(rn.length()-2,2)=="/2")
    {
        return true;
    }
    return false;
}

bool read_is_pair1(string& rn)
{
    if (rn.substr(rn.length()-2,2)=="/1")
    {
        return true;
    }
    return false;
}

bool read_is_pair2(string& rn)
{
    if (rn.substr(rn.length()-2,2)=="/2")
    {
        return true;
    }
    return false;
}

string read_name(string& rn)
{
    if (read_is_paired(rn))
    {
        return rn.substr(0,rn.length()-2);
    }else
    {
        return rn;
    }
}

int number_of_reads_cover(Strain& a, Strain& b)
{
    int n = 0;
    set<string> visited;
    for (int i=0; i<a.assign_reads.size(); ++i)
    {
        for (int j=0; j<b.assign_reads.size(); ++j)
        {
            string rni0 = get<0>(a.assign_reads[i]);
            string rnj0 = get<0>(b.assign_reads[j]);
            string rni = read_name(rni0);
            string rnj = read_name(rnj0);
            DoubleL di = get<1>(a.assign_reads[i]);
            DoubleL dj = get<1>(b.assign_reads[j]);
            if (visited.count(rni)>0) continue;
            if (rni==rnj and !read_is_paired(rni0))
            {
                n += 1;
                visited.emplace(rni);
            }else if (rni==rnj and read_is_pair1(rni0) and read_is_pair2(rnj0))
            {
                n += 1;
                visited.emplace(rni);
            }else if (rni==rnj and read_is_pair2(rni0) and read_is_pair1(rnj0))
            {
                n += 1;
                visited.emplace(rni);
            }else if (rni==rnj)
            {
                n += 1;
                visited.emplace(rni);
            }
        }
    }

    return n;
}

int main(int argc, char** argv)
{
    string cmd;

    sc_parameter *parameters = new sc_parameter;

    // parse command line
    sc_parse_cmd_line(--argc,++argv,parameters);

    // print help message or not
    if (parameters->print_help or argc==0)
    {
        help();
        exit(0);
    }

    // make scan window
    vector<sc_window> windows;
    make_scan_window(parameters,windows);

    // scan through window
    for (auto it=windows.begin(); it!=windows.end(); ++it)
    {
        string gene_roi = it->gn + ":" + to_string(it->p0) + "-" + to_string(it->p1);

        // load gene file
        GenomeSeq gene_seq;
        load_gene_seq(parameters->gene_file,gene_roi,gene_seq);

        // load mapping reads
        vector<AlignRead> reads;
        vector<ReadAux> aux_reads;
        ReadPairs read_pairs;
        load_mapping_reads(gene_seq,parameters->mapping_file,parameters->mapping_qual,parameters->read_len,
                           parameters->max_ins,parameters->max_depth,gene_roi,reads,aux_reads,read_pairs);

        // skip if the number of reads less than 1
        if (reads.empty())
        {
            continue;
        }

        // construct pog
        vector<Strain> strains;

        PartialOrderGraph* pog = new PartialOrderGraph(gene_seq,reads);

        if (!parameters->plot_graph)
        {
            pog->infer_strains(strains,read_pairs,5000,parameters->error_rate,parameters->tau,
                               parameters->diff_rate);

            pog->read_assign(strains,reads,read_pairs,5000);

            // sort according to abundance
            sort(strains.begin(),strains.end(),[](Strain& a, Strain& b){return a.abundance>b.abundance;});

            // output local results
            GenomeSeq gene_seq0;
            load_gene_seq(parameters->gene_file,it->gn,gene_seq0);
            int si = 0;
            for (auto s=strains.begin(); s!=strains.end(); ++s,++si)
            {
                if (s->abundance>=parameters->tau)
                {   
                    // change by Feng Zeng, Jan 22, 2016
                    GenomeSeq out_seq = gene_seq0;
                    //out_seq.replace(it->p0-1, it->p1-it->p0+1, s->plain_seq());
                    out_seq = s->plain_seq();

                    cout << ">contig" << it->gn << "" << it->p0 << "" << it->p1 << "" << si 
                         << endl;
                    cout << out_seq << endl;
                }
            }
        }else
        {

            pog->output_edge(cout);
        }

        delete pog;
    }

    delete parameters;

    return 0;
}

