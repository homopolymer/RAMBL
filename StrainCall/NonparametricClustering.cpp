#include "PartialOrderGraph.hpp"
#include <random>
#include <algorithm>
#include <cmath>
#include <queue>
#include <vector>
#include <set>
using namespace std;

void normalize(vector<DoubleL>& f)
{
    DoubleL z = 0;
    for_each(f.begin(),f.end(),[&z](DoubleL x){z+=x;});
    for_each(f.begin(),f.end(),[z](DoubleL& x){x/=z;});
}

void PartialOrderGraph::hard_clustering(vector<Strain>& strains,
                                        vector<ReadBase>& reads,
                                        vector<int>& new_reads,
 					ReadPairs& read_pairs
                                        )
{
    int c,i,j,m,n,id,uid,ri,si,cn;

    vector<DoubleL> abundance(strains.size(),0);

    vector<SubstitutionCount> substitute(strains.size());

    // assign read to strain
    ri = 0;
    for (auto r=reads.begin(); r!=reads.end(); ++r,++ri)
    {
        id = get<0>(*r);
        cn = get<2>(*r);
        // loop over reads
        for (;cn>0;cn--)
        {
            // paired id
            uid = read_pairs[id][cn-1];

            // compute posterior
            vector<DoubleL> p(strains.size(),0);
            // strain prior
            si=0;
            for (auto s=strains.begin(); s!=strains.end(); ++s,++si)
            {
                p[si] = s->abundance;
            }
            normalize(p);
            // read likelihood
            si = 0;
            for (auto s=strains.begin(); s!=strains.end(); ++s,++si)
            {
                p[si] = log(p[si]) + s->logprob(id);
                if (uid>=0)
                {
                    p[si] += s->logprob(uid);
                }
                p[si] = exp(p[si]);
            }
            normalize(p);
        
            // update counting
            for (si=0; si<strains.size(); ++si)
            {
                // abundance
                abundance[si] += p[si];

                // substitution
                if (get<1>(*r).length()==1)
                {
                    Substitution sb(get<1>((*strains[si].path.rbegin())->state),
                                    get<1>(*r));
                    if (substitute[si].count(sb)==0)
                    {
                        substitute[si][sb] = p[si];
                    }else
                    {
                        substitute[si][sb] += p[si];
                    }
                }else if (new_reads[ri])
                {
                    i = get<1>((*strains[si].path.rbegin())->state).length();
                    j = get<1>(*r).length();
                    while (i>0 and j>0)
                    {
                        Substitution sb(string(1,get<1>((*strains[si].path.rbegin())->state)[--i]),
                                        string(1,get<1>(*r)[--j]));
                        if (substitute[si].count(sb)==0)
                        {
                            substitute[si][sb] = p[si];
                        }else
                        {
                            substitute[si][sb] += p[si];
                        }
                    }
                }else
                {
                    m = get<1>((*strains[si].path.rbegin())->state).length();
                    n = get<1>(*r).length();
                    i = 0;
                    j = 0;
                    while (i<m and j<n)
                    {
                        Substitution sb(string(1,get<1>((*strains[si].path.rbegin())->state)[i++]),
                                        string(1,get<1>(*r)[j++]));
                        if (substitute[si].count(sb)==0)
                        {
                            substitute[si][sb] = p[si];
                        }else
                        {
                            substitute[si][sb] += p[si];
                        }
                    }
                }
            }
        }
    }

    // update strain model
    for (c=0; c<strains.size(); ++c)
    {
        strains[c].update_model(abundance[c],substitute[c]);
    }
}

typedef discrete_distribution<> MultinomialDist;
void PartialOrderGraph::np_bayes_clustering(vector<Strain>& strains,
                                            vector<ReadBase>& reads, ReadPairs& read_pairs, int n, 
                                            vector<DoubleL>& abundance,
                                            vector<DoubleL>& posterior)
{
    int i,j,s,c,id,cn,k,m=reads.size(),S=strains.size();
    vector<DoubleL> a0,a;
    vector<DoubleL> p;
    set<int> clusters;
    vector<SubstitutionCount> sc;
    vector<DoubleL> posterior1(S,0);

    std::random_device rd;
    //std::mt19937 gen(rd());
    std::mt19937 gen(1234);

    for (auto it=strains.begin(); it!=strains.end(); ++it)
    {
        a0.push_back(it->abundance);
        a.push_back(it->abundance);
        p.push_back(0);
    }

    sc.resize(S);

    int read_size = 0;
    for (auto it=reads.begin(); it!=reads.end(); ++it)
    {
        read_size += get<2>(*it);
    }

    // sampling n times
    n = min(n,40000/read_size);
    for (i=0; i<n; i++)
    {
        // scan through reads
        for (j=0; j<m; j++)
        {
            id = get<0>(reads[j]);
            cn = get<2>(reads[j]);

            for (;cn>0;cn--)
            {
                // strain prior
                for (s=0; s<S; s++)
                {
                    p[s] = a[s];
                }
                normalize(p);
                // strain posterior
                for (s=0; s<S; s++)
                {
                    p[s] = log(p[s])+strains[s].logprob(id);
                    // paired read
                    if (read_pairs[id][cn-1]>=0)
                    {
                        int uid = read_pairs[id][cn-1];
                        if (strains[s].read_loglik.count(uid)>0)
                        {
                            p[s] += strains[s].logprob(uid);
                        }
                    }
                    p[s] = exp(p[s]);
                }
                // multinomial sampling
                MultinomialDist sampler(p.begin(),p.end());
                c = sampler(gen);
                // update a
                a[c] += 1;
                // update sc
                Substitution sub(get<1>((*strains[c].path.rbegin())->state),
                                 get<1>(reads[j]));
                if (sc[c].count(sub)==0)
                {
                    sc[c][sub] = 1;
                }else
                {
                    sc[c][sub] += 1;
                }
                // update cluster
                clusters.emplace(c);
            }   
        }
        // update k posterior
        k = clusters.size();
        posterior1[k-1] += 1;
        clusters.clear();
    }

    // divide by n
    normalize(a);
    for_each(a.begin(),a.end(),[read_size](DoubleL& x){x*=read_size;});
    for_each(sc.begin(),sc.end(),[n](SubstitutionCount& sc)
            {
                for(auto it=sc.begin(); it!=sc.end(); ++it)
                {
                    it->second /= n;
                }
            });
    for_each(posterior1.begin(),posterior1.end(),[n](DoubleL& x){x/=n;});

    // copy a to abundance
    abundance.resize(a.size());
    copy(a.begin(),a.end(),abundance.begin());

    // update posterior
    for (s=0; s<S; ++s)
    {
        posterior[s] += posterior1[s];
    }

    // update strain model
    for (s=0; s<S; ++s)
    {
        strains[s].update_model(a[s],sc[s]);
    }
}

DoubleL Qx(vector<DoubleL> a, int n)
{
    sort(a.begin(),a.end(),[](DoubleL x, DoubleL y){return x>y;});
    if (n>=a.size())
    {
        return *a.rbegin();
    }
    return a[n];
}

int seq_diff_num(string a, string b);
int seq_len(string a, string b);
DoubleL seq_identity(string a, string b);
void merge_strains(vector<Strain>& strains, DoubleL diff);
void read_reassign(vector<Strain>& strains, set<tuple<int,int>>& reads);

void PartialOrderGraph::streaming_clustering(vector<Strain>& strains, ReadPairs& read_pairs,
                                             int n, DoubleL e, DoubleL tau, DoubleL diff)
{
    int level = 0;
    bool branching = false;
    PartialOrderGraphNode *u,*v;
    vector<Strain> level_strains;
    vector<Strain> sublevel_strains;
    queue<PartialOrderGraphNode*> level_node;
    queue<PartialOrderGraphNode*> sublevel_node;
    set<PartialOrderGraphNode*> has_visited;
    vector<ReadBase> level_reads;
    int level_read_count;
    vector<DoubleL> abundance;
    vector<DoubleL> posterior(100,0);
    set<tuple<int,int>> total_reads;
    vector<int> new_reads;
    DoubleL Z,Z0,Zt,Zt0;

    level_strains.push_back(Strain(100,e));
    level_node.push(nodes[0]);
    level_read_count = 0;
    while (!level_node.empty())
    {
        // debug
        if (0){
        if (!level_strains.empty()){
            cerr << "------------------------------" << endl
                 << "before clustering" << endl
                 << "level: " << level << endl;

            for (auto ls_iter = level_strains.begin(); ls_iter!=level_strains.end(); ls_iter++)
            {
                cerr << ls_iter->strain_seq() << "\t" << ls_iter->abundance << endl;
            }
        }
        }

        // ---------------------------------
        // this part is level-order visiting
        // ---------------------------------
        u = level_node.front(); level_node.pop();
        if (u==nodes[0])
        {
            // initialize a strain
            level_strains[0].path_extend(u);
            level_strains[0].abundance = 1;
        }else if (get<1>(u->state)=="$")
        {
            read_reassign(level_strains,total_reads);
            merge_strains(level_strains,diff);
            strains.resize(level_strains.size());
            copy(level_strains.begin(),level_strains.end(),strains.begin());
        }else 
        {
            // collect reads at the level
            for (auto it=u->read_pool.begin(); it!=u->read_pool.end(); ++it)
            {
                level_reads.push_back(*it);
                level_read_count += get<2>(*it);
            }
        }

        // collect nodes at next level
        for (auto it=u->out.begin(); it!=u->out.end(); ++it)
        {
            v = *it;
            if (has_visited.count(v)==0)
            {
                sublevel_node.push(v);
                has_visited.emplace(v);
            }
        }

        if (level_node.empty())
        {
            // -----------------------------
            // this part is strain inference
            // -----------------------------

            // update strain read likelihood
            if (!level_reads.empty())
            {
                new_reads.resize(level_reads.size());
                fill(new_reads.begin(),new_reads.end(),0);
                for (auto s=level_strains.begin(); s!=level_strains.end(); ++s)
                {
                    int ri = 0;
                    for (auto r=level_reads.begin(); r!=level_reads.end(); ++r,++ri)
                    {
                        int rid = get<0>(*r);
                        string rb = get<1>(*r);
                        string sb = get<1>((*(s->path.rbegin()))->state);
                        DoubleL loglik;
                        if (sb.length()==1)
                        {
                            if (sb=="N") sb=rb;
                            loglik = s->logprob(sb,rb);
                        }else
                        {
                            loglik = 0;
                            int ii,jj;
                            if(s->read_loglik.count(rid)==0) // start at this level
                            {
                                ii=sb.length();
                                jj=rb.length();
                                while(ii>0 and jj>0)
                                {
                                    string ssb = string(1,sb[--ii]);
                                    string rrb = string(1,rb[--jj]);
                                    if (ssb=="N") ssb = rrb;
                                    loglik += s->logprob(ssb,rrb);
                                }
                                new_reads[ri] = 1;
                            }else // start before this level
                            {
                                ii=0;
                                jj=0;
                                while(ii<sb.length() and jj<rb.length())
                                {
                                    string ssb = string(1,sb[ii++]);
                                    string rrb = string(1,rb[jj++]);
                                    if (ssb=="N") ssb = rrb;
                                    loglik += s->logprob(ssb,rrb);
                                }
                            }
                        }
                        s->update_read_loglik(rid,loglik);
                    }
                }
                // collect reads
                for (auto r=level_reads.begin(); r!=level_reads.end(); ++r)
                {
                    tuple<int,int> rr(get<0>(*r),get<2>(*r));
                    if (total_reads.count(rr)==0)
                        total_reads.emplace(rr);
                }
            }
            
            // clustering if having branches
            if (branching and !level_reads.empty())
            {
                map<string,DoubleL> A_prior,A_posterior;
                // abundance before clustering
                for (auto lsiter = level_strains.begin(); lsiter!=level_strains.end(); lsiter++)
                {
                    A_prior[lsiter->strain_seq()] = lsiter->abundance;
                }


                // probabilistic clustering
                np_bayes_clustering(level_strains,level_reads,read_pairs,n,abundance,posterior);


                // abundance after clustering
                for (auto lsiter = level_strains.begin(); lsiter!=level_strains.end(); lsiter++)
                {
                    A_posterior[lsiter->strain_seq()] = lsiter->abundance;
                }

                // abundance delta
                DoubleL A_delta_max = 0;
                for (auto lsiter = level_strains.begin(); lsiter!=level_strains.end(); lsiter++)
                {
                    DoubleL delta = A_posterior[lsiter->strain_seq()] - A_prior[lsiter->strain_seq()];
                    if (A_delta_max < delta) A_delta_max = delta;
                }

                // terminate a strain if its abundance level is lower than tau
                Z = 0;
                for_each(abundance.begin(),abundance.end(),[&Z](DoubleL z){Z+=z;});

                Zt = Z*tau;
 
                vector<vector<Strain>::iterator> delete_strains;

                // collect strains to be deleted
                int si=0;
                for (auto s=level_strains.begin(); s!=level_strains.end(); ++s,++si)
                {

                    DoubleL delta_abundance = A_posterior[s->strain_seq()] - A_prior[s->strain_seq()];

                    if (abundance[si]<Zt or delta_abundance<0.01*A_delta_max)
                    {
                        delete_strains.push_back(s);
                    }
                }

                for (auto s=delete_strains.rbegin(); s!=delete_strains.rend(); ++s)
                {
                    level_strains.erase(*s);
                }//*/
            }else if (!level_reads.empty())
            {
                hard_clustering(level_strains,level_reads,new_reads,read_pairs);
            }

if (0){
        if (!level_strains.empty()){
            cerr << "------------------------------" << endl
                 << "after clustering" << endl
                 << "level: " << level << endl;

            for (auto ls_iter = level_strains.begin(); ls_iter!=level_strains.end(); ls_iter++)
            {
                cerr << ls_iter->strain_seq() << "\t" << ls_iter->abundance << endl;
            }
        }
        }

            // generate strains at next level
            branching = false;
            int si = 0;
            for (auto s=level_strains.begin(); s!=level_strains.end(); ++s,++si)
            {
                v = *(s->path.rbegin());

                DoubleL oz = 0,zz,tt;
                vector<DoubleL> oc;
                DoubleL moc = 0;
                for (auto it=v->out.begin(); it!=v->out.end(); ++it)
                {
                    int oc0 = number_of_reads_cover_nodes(v,*it);
                    oc.push_back(oc0);
                    oz += oc0;
                    if (moc<oc0) moc = oc0;
                }

                int oi = 0, dd = 0;
                for (auto it=v->out.begin(); it!=v->out.end(); ++it,++oi)
                {
                    zz = max(oz,(DoubleL)level_read_count);
                    tt = max((DoubleL)0.05,tau);

                    if (get<1>((*it)->state)!="$" and oz>0)
                    {
                        if (oc[oi]<=1. and oc[oi]<moc) 
                        {
                            dd += 1;
                            continue;
                        }

                        Strain ns = *s;
                        ns.path_extend(*it);
                        
                        if (oc[oi]>0)
                        {
                            ns.abundance = s->abundance*oc[oi]/oz;
                        }
                        else
                        {
                            ns.abundance = oz*min(0.01,(double)tau);
                        }
                        sublevel_strains.push_back(ns);
                    }else
                    {
                        Strain ns = *s;
                        ns.path_extend(*it);
                        ns.abundance = s->abundance;
                        sublevel_strains.push_back(ns);
                    }
                }
                if (outdegree(v)>1+dd)
                {
                    branching = true;
                }
            }

            // maker sure the candidate set not too large
            if (sublevel_strains.size()>80)
            {
                vector<DoubleL> ssa;
                for_each(sublevel_strains.begin(),sublevel_strains.end(),[&ssa](Strain& x){ssa.push_back(x.abundance);});
                Zt0 = Qx(ssa,80);

                vector<vector<Strain>::iterator> delete_strains;
                for (auto s=sublevel_strains.begin(); s!=sublevel_strains.end(); ++s)
                {
                    if (s->abundance<Zt0)
                    {
                        delete_strains.push_back(s);
                    }
                }

                for (auto s=delete_strains.rbegin(); s!=delete_strains.rend(); ++s)
                {
                    sublevel_strains.erase(*s);
                }
            }

            // ------------------------------------------
            // this part is to update graph walking state
            // ------------------------------------------
            
            // update level
            level += 1;
            // update level node
            while(!sublevel_node.empty())
            {
                v = sublevel_node.front(); sublevel_node.pop();
                level_node.push(v);
            }
            // update level strains
            level_strains.clear();
            level_strains.resize(sublevel_strains.size());
            copy(sublevel_strains.begin(),sublevel_strains.end(),level_strains.begin());
            sublevel_strains.clear();
            sublevel_strains.resize(0);
            has_visited.clear();
            // clean level reads
            level_reads.clear();
            level_reads.resize(0);
            level_read_count = 0;
        }
    }

//    strains.resize(level_strains.size());
//    copy(level_strains.begin(),level_strains.end(),strains.begin());
//
}

DoubleL seq_identity(string a, string b)
{
    int iden = 0, len = 0;
    for (int i=0; i<a.length(); ++i)
    {
        if (a[i]=='-' and b[i]=='-')
        {
            continue;
        }else if (a[i]=='=' and b[i]=='=')
        {
            continue;
        }else if (a[i]=='=' and b[i]=='-')
        {
            continue;
        }else if (a[i]=='-' and b[i]=='=')
        {
            continue;
        }else if (a[i]=='^' and b[i]=='^')
        {
            continue;
        }else if (a[i]==b[i])
        {
            iden += 1;
        }
        len += 1;
    }

    return (iden+0.0)/len;
}

int seq_len(string a, string b)
{
    int len = 0;
    for (int i=0; i<a.length(); ++i)
    {
        if ((a[i]=='-' or a[i]=='=') and (b[i]=='-' or b[i]=='='))
        {
            continue;
        }
        len += 1;
    }
    return len;
}

int seq_diff_num(string a, string b)
{
    int diff = 0;
    for (int i=0; i<a.length(); ++i)
    {
        if ((a[i]=='-' or a[i]=='=') and (b[i]=='-' or b[i]=='='))
        {
            continue;
        }
        if (a[i]!=b[i])
        {
            diff += 1;
        }
    }
    return diff;
}

void merge_strains(vector<Strain>& strains, DoubleL diff)
{
    sort(strains.begin(),strains.end(),[](Strain& a, Strain& b){return a.abundance>b.abundance;});
    
    int i,j;
    vector<Strain> merged(1,strains[0]);   
    for (i=1; i<strains.size(); i++)
    {
        for (j=0; j<merged.size(); j++)
        {
            if (seq_identity(strains[i].strain_seq(),merged[j].strain_seq())>1-diff)
            {
                merged[j].abundance += strains[i].abundance;
                break;
            }
        }
        
        if (j==merged.size())
        {
            merged.push_back(strains[i]);
        }
    }

    strains.resize(merged.size());
    copy(merged.begin(),merged.end(),strains.begin());
}

void read_reassign(vector<Strain>& strains,set<tuple<int,int>>& reads)
{
    int c;
    sort(strains.begin(),strains.end(),[](Strain& a,Strain& b){return a.abundance>b.abundance;});
    vector<int> count(strains.size(),0);
    
    for (auto it=reads.begin();
              it!=reads.end(); ++it)
    {
        DoubleL max_loglik = numeric_limits<DoubleL>::lowest(), loglik;
        for (auto s=strains.begin(); s!=strains.end(); ++s)
        {
            loglik = s->logprob(get<0>(*it));
            if (max_loglik<loglik)
            {
                max_loglik = loglik;
            }
        }

        c=0;
        for (auto s=strains.begin(); s!=strains.end(); ++s,++c)
        {
            loglik = s->logprob(get<0>(*it));
            if (max_loglik==loglik)
            {
                count[c] += get<1>(*it);
                break;
            }
        }
    }
}

void PartialOrderGraph::infer_strains(vector<Strain>& strains, ReadPairs& read_pairs,
                                      int n, DoubleL e, DoubleL tau, DoubleL diff)
{
    streaming_clustering(strains,read_pairs,n,e,tau,diff);
}


void PartialOrderGraph::read_assign(vector<Strain>& strains, vector<ReadAux>& reads, int n)
{
    n = min(n,40000/(int)reads.size());
    vector<vector<DoubleL>> assign(reads.size(),vector<DoubleL>(strains.size(),0));

    std::random_device rd;
    //std::mt19937 gen(rd());
    std::mt19937 gen(1234);

    int id,i,j,c,m;
    vector<DoubleL> a,p;
    for (auto s=strains.begin(); s!=strains.end(); ++s)
    {
        a.push_back(s->abundance);
        p.push_back(0);
    }

    for (;n>=0;--n)
    {
        j = 0;
        for (auto r=reads.begin(); r!=reads.end(); ++r,++j)
        {
            int id = get<0>(*r);

            // prior
            copy(a.begin(),a.end(),p.begin());
            normalize(p);

            // posterior
            for (i=0; i<strains.size(); i++)
            {
                p[i] = exp(log(p[i])+strains[i].logprob(id));
            }                      

            // Multinomial sampling
            MultinomialDist sampler(p.begin(),p.end());
            c = sampler(gen);

            a[c] += 1;
            assign[j][c] += 1;
        }
    }

    for (j=0; j<reads.size(); j++)
    {
        normalize(assign[j]);
        auto it = max_element(assign[j].begin(),assign[j].end());
        for (i=0; i<strains.size(); i++)
        {  
            if (assign[j][i]==*it)
            {
                tuple<string,DoubleL> res(get<1>(reads[j]),assign[j][i]);
                strains[i].assign_reads.push_back(res);
                break;
            }
        }
    }

    normalize(a);
    for (i=0; i<strains.size(); i++)
    {
        strains[i].abundance = a[i];
    }
}

void PartialOrderGraph::read_assign(vector<Strain>& strains, vector<AlignRead>& reads, ReadPairs& read_pairs, int n)
{
    int read_size = 0;
    for_each(reads.begin(),reads.end(),[&read_size](AlignRead& ar){read_size+=get<4>(ar);});

    n = min(n,40000/(int)read_size);

    std::random_device rd;
    //std::mt19937 gen(rd());
    std::mt19937 gen(1234);

    int id,i,j,c,m,cn,uid;
    vector<DoubleL> a,p;
    for (auto s=strains.begin(); s!=strains.end(); ++s)
    {
        a.push_back(s->abundance);
        p.push_back(0);
    }

    // sampling
    for (;n>0;n--)
    {
        j=0;
        id=0;
        for (auto r=reads.begin(); r!=reads.end(); ++r,++id)
        {
            cn = get<4>(*r);
            // scan through reads within a group
            for (;cn>0;--cn,++j)
            {
                uid = read_pairs[id][cn-1];
                // prior
                copy(a.begin(),a.end(),p.begin());
                normalize(p);
                // posterior
                i = 0;
                for (auto s=strains.begin(); s!=strains.end(); ++s,++i)
                {
                     p[i] = log(p[i]) + s->logprob(id);
                     if (uid>=0)
                     {
                          p[i] += s->logprob(uid);
                     }
                     p[i] = exp(p[i]);
                }

                // Multinomial sampling
                MultinomialDist sampler(p.begin(),p.end());
                c = sampler(gen);

                a[c] += 1;
            }
        } 
    }

    normalize(a);
    for (i=0; i<strains.size(); i++)
    {
        strains[i].abundance = a[i];
    }
}
