#include "PartialOrderGraph.hpp"
#include <map>
#include <vector>
#include <cmath>
using namespace std;

map<int,string> A = {{0,"A"},{1,"C"},{2,"G"},{3,"T"},{4,"-"},{5,"="}};

Strain::Strain()
{
    int i,j;
    for (i=0; i<6; i++)
    {
        for (j=0; j<6; j++)
        {
            if (i==j)
            {
                sub_count[Substitution(A[i],A[j])] = 100;
            }else
            {
                sub_count[Substitution(A[i],A[j])] = 1;
            }
        }
    }

    Z = 0;
    for (i=0; i<6; i++)
    {
        comp_count[A[i]] = 0;
        for (j=0; j<6; j++)
        {
            comp_count[A[i]] += sub_count[Substitution(A[i],A[j])];
        }

        Z += comp_count[A[i]];
    }

    abundance = 0;
}

Strain::Strain(int N,DoubleL e)
{
    int i,j;
    for (i=0; i<6; i++)
    {
        for (j=0; j<6; j++)
        {
            if (i==j)
            {
                sub_count[Substitution(A[i],A[j])] = N*(1-e);
            }else
            {
                sub_count[Substitution(A[i],A[j])] = N*e;
            }
        }
    }

    Z = 0;
    for (i=0; i<6; i++)
    {
        comp_count[A[i]] = 0;
        for (j=0; j<6; j++)
        {
            comp_count[A[i]] += sub_count[Substitution(A[i],A[j])];
        }

        Z += comp_count[A[i]];
    }

    abundance = 0;
}

Strain::Strain(const Strain& other)
{
    abundance = other.abundance;
    sub_count = other.sub_count;
    comp_count = other.comp_count;
    read_loglik = other.read_loglik;   
    read_info = other.read_info;
    assign_reads = other.assign_reads;
    path = other.path;
    Z = other.Z;
}//*/

void Strain::update_read_loglik(int id, DoubleL loglik)
{
    auto it = read_loglik.find(id);
    if (it==read_loglik.end())
    {
        read_loglik[id] = loglik;
    }else
    {
        read_loglik[id] += loglik;
    }
}

void Strain::erase_read(int id)
{
    auto it = read_loglik.find(id);
    if (it!=read_loglik.end())
    {
        read_loglik.erase(it);
    }
}

void Strain::update_model(DoubleL al, SubstitutionCount& sc)
{
    abundance += al;
    
    for (auto it=sc.begin(); it!=sc.end(); ++it)
    {
        sub_count[it->first] += it->second;
    }

    Z = 0;
    for (int i=0; i<6; i++)
    {
        comp_count[A[i]] = 0;
        for (int j=0; j<6; j++)
        {
            comp_count[A[i]] += sub_count[Substitution(A[i],A[j])];
        }
        Z += comp_count[A[i]];
    }
}

DoubleL Strain::logprob(string a)
{
    return (log(comp_count[a])-log(Z));
}

DoubleL Strain::logprob(string a,string b)
{
    return (log(sub_count[Substitution(a,b)])-log(comp_count[a]));
}

DoubleL Strain::operator()(string a)
{
    return exp(logprob(a));
}

DoubleL Strain::operator()(string a, string b)
{
    return exp(logprob(a,b));
}

DoubleL Strain::logprob(int id)
{
    return read_loglik[id];
}

DoubleL Strain::operator()(int id)
{
    return exp(logprob(id));
}

void Strain::path_extend(PartialOrderGraphNode* n)
{
    path.push_back(n);
}

string Strain::strain_seq()
{
    string seq = "";
    for (auto it=path.begin(); it!=path.end(); ++it)
    {
        seq += get<1>((*it)->state);
    }
    return seq;
}

void Strain::clean_read(vector<ReadBase>& current_reads)
{
    typedef map<int,DoubleL>::iterator ll_iter;
    vector<int> delete_reads;
    for (ll_iter it=read_loglik.begin(); it!=read_loglik.end(); ++it)
    {
        auto tt = current_reads.begin();
        for (; tt!=current_reads.end(); ++tt)
        {
            if (it->first==get<0>(*tt))
            {
                break;
            }
        }
        if (tt==current_reads.end())
        {
            delete_reads.push_back(it->first);
        }
    }

    for (auto it=delete_reads.begin(); it!=delete_reads.end(); ++it)
    {
        erase_read(*it);
    }
}

Strain& Strain::operator=(const Strain& other)
{
    abundance = other.abundance;
    sub_count = other.sub_count;
    comp_count = other.comp_count;
    read_loglik = other.read_loglik;
    read_info = other.read_info;
    assign_reads = other.assign_reads;
    path = other.path;
    Z = other.Z;
    return *this;
}//*/

string Strain::plain_seq()
{
    string seq="";
    for (auto p=path.begin(); p!=path.end(); ++p)
    {
        string pl = get<1>((*p)->state);
        if (pl!="^" and pl!="$" and pl!="-" and pl!="=")
        {
            seq += pl;
        }
    }
    return seq;
}

void Strain::update_read_info(int id, SubstitutionCount& c)
{
    auto it = read_info.find(id);
    if (it==read_info.end())
    {
        for (auto tt=c.begin(); tt!=c.end(); ++tt)
        {
            read_info[id][tt->first] = tt->second;
        }
    }else
    {
        for (auto tt=c.begin(); tt!=c.end(); ++tt)
        {
            if (it->second.find(tt->first)==it->second.end())
                it->second[tt->first] = tt->second;
            else it->second[tt->first] += tt->second;
        }
    }

    // update read loglik
    DoubleL loglik = 0;
    for (auto it=read_info[id].begin(); it!=read_info[id].end(); ++it)
    {
        loglik += logprob(get<0>(it->first),get<1>(it->first))*it->second;
    }
    read_loglik[id] = loglik;
}
