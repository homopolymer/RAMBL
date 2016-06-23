#include "PartialOrderGraph.hpp"
#include "MultipleSequenceAlignment.hpp"
#include <algorithm>
#include <iterator>
#include <random>
using namespace std;


//extern typedef string CigarOp;
//typedef int CigarOpLen;
//typedef tuple<CigarOp,CigarOpLen> CigarRecord;

void parse_cigar(AlignCigar& cigar, vector<CigarRecord>& records)
{
    string op_l = "";
    for (auto it=cigar.begin(); it!=cigar.end(); ++it)
    {
        if (*it=='M')
        {
            records.push_back(CigarRecord("M",stoi(op_l)));
            op_l = "";
        }else if (*it=='I')
        {
            records.push_back(CigarRecord("I",stoi(op_l)));
            op_l = "";
        }else if (*it=='D')
        {
            records.push_back(CigarRecord("D",stoi(op_l)));
            op_l = "";
        }else if (*it=='N')
        {
            records.push_back(CigarRecord("N",stoi(op_l)));
            op_l = "";
        }else if (*it=='S')
        {
            records.push_back(CigarRecord("S",stoi(op_l)));
            op_l = "";
        }else if (*it=='H')
        {
            records.push_back(CigarRecord("H",stoi(op_l)));
            op_l = "";
        }else if (*it=='P')
        {
            records.push_back(CigarRecord("P",stoi(op_l)));
            op_l = "";
        }else if (*it=='=')
        {
            records.push_back(CigarRecord("M",stoi(op_l)));
            op_l = "";
        }else if (*it=='X')
        {
            records.push_back(CigarRecord("M",stoi(op_l)));
            op_l = "";
        }else
        {
            op_l += string(1,*it);
        }
    }
}

PartialOrderGraph::PartialOrderGraph(GenomeSeq& G, vector<AlignRead>& R)
{
    N = 0;
    build(G,R);
}

void PartialOrderGraph::build(GenomeSeq& G, vector<AlignRead>& R)
{
    PartialOrderGraphNode* u;
    PartialOrderGraphNode* v;
    PartialOrderGraphNode* w;

    // begin node
    PartialOrderGraphNode* B = new PartialOrderGraphNode(N,PartialOrderGraphNodeState(mat,"^"));
    add_node(B);

    // build the backbone
    int i = 0;
    u = B;
    for (auto it=G.begin(); it!=G.end(); ++it,++i)
    {
        w = new PartialOrderGraphNode(N,PartialOrderGraphNodeState(mat,string(1,*it)));
        add_node(w);
        w->pos = i;
        add_edge(u,w);
        u = w;
    }

    // end node
    PartialOrderGraphNode* E = new PartialOrderGraphNode(N,PartialOrderGraphNodeState(mat,"$"));
    add_node(E);
    add_edge(u,E);

    // scan through alignment
    int rid = 0;
    for (auto it=R.begin(); it!=R.end(); ++it,++rid)
    {

        // preceding node
        u = nodes[get<0>(*it)];
        // current node
        v = nodes[get<0>(*it)+1];
        // pointer on genome 
        int i = get<0>(*it);

        // pointer on read
        int j = 0; 

        int dl = 0;

        // read sequence
        ReadSeq r = get<2>(*it);

        // parse cigar
        vector<CigarRecord> cigar;
        parse_cigar(get<1>(*it),cigar);
        
        // scan through cigar
        for (auto c_it=cigar.begin(); c_it!=cigar.end(); ++c_it)
        {
            auto op  = get<0>(*c_it);
            auto opl = get<1>(*c_it);

            if (op=="S")
            {
                j += j+opl;
                dl = 0;
                continue;
            }else if (op=="M")
            {
                for (int k=0; k<opl+dl; k++,j++)
                {
                    PartialOrderGraphNodeState s;
                    if (G[i]==r[j]) get<0>(s) = mat;
                    else get<0>(s) = mis;
                    get<1>(s) = string(1,r[j]);
                    

                    if (v->state == s)
                    {
                        if (!linking(u,v))
                            add_edge(u,v);

                        v->read_pool.push_back(ReadBase(rid,get<1>(s),get<4>(*it)));

                        // update
                        u = v;
                        v = nodes[++i+1];
                    }else
                    {
                        auto s_it = v->find_sibling(s);
                        if (s_it == v->sibling.end())
                        {
                            w = new PartialOrderGraphNode(N,s);
                            add_node(w);
                            add_edge(u,w);
                            w->read_pool.push_back(ReadBase(rid,get<1>(s),get<4>(*it)));
 
                            v->insert_sibling(w);

                            // update
                            u = w;
                            v = nodes[++i+1];
                        }else
                        {
                            if (!linking(u,*s_it))
                                add_edge(u,*s_it);

                            (*s_it)->read_pool.push_back(ReadBase(rid,get<1>(s),get<4>(*it)));
                    
                            // update
                            u = *s_it;
                            v = nodes[++i+1];
                        }
                    }
                }
                dl = 0;
            }else if (op=="I")
            {
                Gap gap;
                for (int k=0;k<opl;k++,j++)
                {
                    PartialOrderGraphNodeState s(ins,string(1,r[j]));
                    w = new PartialOrderGraphNode(N,s);
                    add_node(w);
                    w->read_pool.push_back(ReadBase(rid,get<1>(s),get<4>(*it)));
                    gap.push_back(w);
                }
                //add_edge(u,v,gap);
                add_edge(u,gap);

                // update
                u = gap.back();
                dl = 0;
            }else if (op=="D")
            {
                Gap gap;
                for (int k=0;k<opl;k++)
                {
                    PartialOrderGraphNodeState s(del,"=");
                    w = new PartialOrderGraphNode(N,s);
                    add_node(w);
                    w->read_pool.push_back(ReadBase(rid,get<1>(s),get<4>(*it)));
                    gap.push_back(w);
                    //update
                    v = nodes[++i+1];
                }
 
                if (get<1>(v->state)=="$")
                {
                    add_edge(u,v,gap);
                    // update
                    u = v;
                    continue;
                }
                
                add_edge(u,gap);
                /*
                PartialOrderGraphNodeState s;
                if (G[i]==r[j]) get<0>(s) = mat;
                else get<0>(s) = mis;
                get<1>(s) = r[j];

                if (s==v->state)
                {
                    add_edge(u,v,gap);
                    v->read_pool.push_back(ReadBase(rid,get<1>(s),get<4>(*it)));
                    // update
                    u = v;
                    v = nodes[++i+1];
                }else
                {
                    w = new PartialOrderGraphNode(N,s);
                    add_node(w);
                    add_edge(u,w,gap);
                    w->read_pool.push_back(ReadBase(rid,get<1>(s),get<4>(*it)));

                    v->insert_sibling(w);

                    // update
                    u = w;
                    v = nodes[++i+1];
                }//*/

                u = gap.back();
                //j++;
                //dl = -1;
                dl = 0;
            }
        }

        if (!linking(u,v) and u!=v)
            add_edge(u,v);

    }

    // canonization
    canonize_graph();

    // path collapse
    path_collapse();

    // compute node level
    node_level();
}

void PartialOrderGraph::add_node(PartialOrderGraphNode* w)
{
    N += 1;
    nodes.push_back(w);
}

void PartialOrderGraph::add_edge(PartialOrderGraphNode* u,
                                 PartialOrderGraphNode* w)
{
    u->insert_out(w);
    w->insert_in(u);
}


void PartialOrderGraph::add_edge(PartialOrderGraphNode* u,
                                 Gap& gap)
{
    PartialOrderGraphNode* a = u;

    auto it = gap.begin();
    for (;it!=gap.end();++it)
    {
        add_edge(a,*it);
        a = *it;
    }
}

void PartialOrderGraph::add_edge(PartialOrderGraphNode* u,
                                 PartialOrderGraphNode* v,
                                 Gap& gap)
{
    PartialOrderGraphNode* a = u;
    PartialOrderGraphNode* b = v;

    auto it = gap.begin();
    for (;it!=gap.end();++it)
    {
        add_edge(a,*it);
        a = *it;
    }

    add_edge(gap.back(),b);
}

void PartialOrderGraph::delete_edge(PartialOrderGraphNode* u,
                                    PartialOrderGraphNode* v)
{
    u->delete_out(v);
    v->delete_in(u);
}

void PartialOrderGraph::output_edge(ostream& f)
{
    for (auto it=nodes.begin();it!=nodes.end();++it)
    {
        int rc = 0;
        for_each((*it)->read_pool.begin(),(*it)->read_pool.end(),[&rc](ReadBase rb)
                {rc+=get<2>(rb);});
        f << "#" << "\t" << (*it)->id << "\t" << (*it)->level << "\t" << get<1>((*it)->state) 
          << "\t" << rc << endl;
    }
    for (auto it=nodes.begin();it!=nodes.end();++it)
    {
        auto o_it = (*it)->out.begin();
        for (;o_it!=(*it)->out.end();++o_it)
        {
            int nr = number_of_reads_cover_nodes(*it,*o_it);
            f << (*it)->id << "\t" << (*o_it)->id << "\t" << nr << endl;
        }
    }
}

bool PartialOrderGraph::linking(PartialOrderGraphNode* u, PartialOrderGraphNode* v)
{
    for (auto it=u->out.begin(); it!=u->out.end(); ++it)
    {
        if (v==(*it))
        {
            return true;
        }
    }
    return false;
//    auto it = u->find_out(v->state);
//    if (it==u->out.end())
//        return false;
//    return true;
}

void PartialOrderGraph::find_insert_from(PartialOrderGraphNode* u, vector<GapEx>& inserts)
{
    int c;
    Gap g;
    PartialOrderGraphNode* v;
    stack<PartialOrderGraphNode*> to_visit;

    to_visit.push(u);
    while (!to_visit.empty())
    {
        v = to_visit.top(); to_visit.pop();

        if (v==u)
        {
            for (auto it=v->out.begin(); it!=v->out.end(); ++it)
            {
                if (get<0>((*it)->state)==ins)
                {
                    to_visit.push(*it);
                }
            }
        }else if (get<0>(v->state)==mat || get<0>(v->state)==mis)
        {
            inserts.push_back(GapEx(u,v,g));
            g.clear();
            g.resize(0);
        }else 
        {
            g.push_back(v);
            for (auto it=v->out.begin(); it!=v->out.end(); ++it)
            {
                //if (get<0>((*it)->state)==ins)
                {
                    to_visit.push(*it);
                }
            }
        }
    }
}

void PartialOrderGraph::find_insert_at_level(int i, vector<GapEx>& inserts)
{
    PartialOrderGraphNode* u = nodes[i];
    find_insert_from(u,inserts);
    for (auto it=u->sibling.begin(); it!=u->sibling.end(); ++it)
    {
        find_insert_from(*it,inserts);
    }
}


void PartialOrderGraph::delete_node(PartialOrderGraphNode* w, bool bridging)
{
    for (auto it=w->in.begin(); it!=w->in.end(); ++it)
    {
        for (auto itt=w->out.begin(); itt!=w->out.end(); ++itt)
        {
            if (!linking(*it,*itt) and bridging)
            {
                add_edge(*it,*itt);
            }
            (*itt)->delete_in(w);
        }
        (*it)->delete_out(w);
    }   

    int w_id = w->id;

    auto it = nodes.begin();
    for (; it!=nodes.end(); ++it)
    {
        if ((*it)->id==w_id)
        {
            break;
        }
    }
    if (it!=nodes.end())
    {
        nodes.erase(it);
    }

    N -= 1;

    for (int i=w_id; i<N; ++i)
    {
        nodes[i]->id -= 1;
    }

    deleted_nodes.push_back(w);
}

void PartialOrderGraph::canonize_insert_at_level(int i)
{

    MultipleSequenceAlignmentSP
    <Index2D,
     SimpleScoreModel,
     vector,
     string,
     char> msa;
    MSA<vector,char> result;

    // collect sequences
    int n=0,l=0,k=1e9,d,t;
    vector<GapEx> inserts;
    find_insert_at_level(i,inserts);

    if (inserts.empty())
        return;

    // sort inserts
    sort(inserts.begin(), inserts.end(), [](GapEx& a, GapEx& b){return get<2>(a).size()>get<2>(b).size();});

    vector<string> seqs;
    for (auto it=inserts.begin(); it!=inserts.end(); ++it)
    {
        Gap g = get<2>(*it);
        string seq="";
        for (auto itt=g.begin(); itt!=g.end(); ++itt)
        {
            seq += get<1>((*itt)->state);
        }
        seqs.push_back(seq);

        n += 1;
        if (seq.length()>l)
        {
            l = seq.length();
        }
        if (seq.length()<k)
        {
            k = seq.length();
        }
    }
    d = l-k;

    PartialOrderGraphNode *w;
    if (n==1)
    {
        add_edge(i,l);
        delete_edge(i);
    }else if (d==0)
    {
        add_edge(i,l);
        delete_edge(i);
    }else
    {
        // WIRED: if temp vector were removed, msa runs wrong and slow
        //        Something wrong!!!
        vector<int> temp;

        // run msa
        msa.align(seqs,result);
        
        // scan through inserts
        for (t=0; t<n; ++t)
        {
            vector<char> res;
            string new_seq="";
            result.get(t,res);
            for (auto it=res.begin(); it!=res.end(); ++it)
            {
                new_seq += string(1,*it);
            }

            if (new_seq!=seqs[t])
            {

                int rid,rcn;
                // delete old inserted nodes
                for (auto it=get<2>(inserts[t]).begin();
                     it!=get<2>(inserts[t]).end(); ++it)
                {
                    rid = get<0>((*it)->read_pool[0]);
                    rcn = get<2>((*it)->read_pool[0]);
                    delete_node(*it);
                }

                // add new inserted nodes
                Gap g;
                for (auto it=res.begin(); it!=res.end(); ++it)
                {
                    w = new PartialOrderGraphNode(N,PartialOrderGraphNodeState(ins,string(1,*it)));
                    add_node(w);
                    w->read_pool.push_back(ReadBase(rid,string(1,*it),rcn));
                    g.push_back(w);
                }
                add_edge(get<0>(inserts[t]),get<1>(inserts[t]),g);
            }
        }

        l = result.size();
        add_edge(i,l);
        delete_edge(i);       
    }
}


void PartialOrderGraph::canonize_insert()
{
    int i;
    for (i=0; i<N; ++i)
    {
        if (get<1>(nodes[i]->state)=="$")
        {
            break;
        }
        canonize_insert_at_level(i);
    }
}

int PartialOrderGraph::distance(PartialOrderGraphNode* u,
                                PartialOrderGraphNode* v)
{
}

int PartialOrderGraph::node_level_exclude_delete(PartialOrderGraphNode* w)
{
    stack<PartialOrderGraphNode*> level_node;
    stack<PartialOrderGraphNode*> sublevel_node;
    set<PartialOrderGraphNode*> visited_node;

    int level = 0;
    PartialOrderGraphNode *u,*v;
    level_node.push(nodes[0]);
    while (!level_node.empty())
    {
        u = level_node.top(); level_node.pop();
        if (u==w)
        {
            break;
        }

        for (auto it=u->out.begin(); it!=u->out.end(); ++it)
        {
            if (get<0>((*it)->state)==del)
            {
                continue;
            }
            sublevel_node.push(*it);
        
            // avoid the case: mat->del->mis
            v = *it;
            for (auto itt=v->sibling.begin(); itt!=v->sibling.end(); ++itt)
            {
                sublevel_node.push(*itt);
            }
        }

        if (level_node.empty())
        {
            while (!sublevel_node.empty())
            {
                v = sublevel_node.top(); sublevel_node.pop();
                if (visited_node.find(v) != visited_node.end())
                {
                    continue;
                }               
                level_node.push(v);
                visited_node.emplace(v);
            }
            level += 1;
            visited_node.clear();
        }
    }

    return level;
}

void PartialOrderGraph::find_delete_from(PartialOrderGraphNode* w, vector<GapEx>& deletes)
{
    int c;
    Gap gap;
    PartialOrderGraphNode* u;
    stack<PartialOrderGraphNode*> to_visit;
    stack<int> visited_count;
    to_visit.push(w);
    visited_count.push(0);
    while(!to_visit.empty())
    {
        u = to_visit.top(); to_visit.pop();
        c = visited_count.top(); visited_count.pop();

        if (u == w)
        {
            for (auto it=u->out.begin(); it!=u->out.end(); ++it)
            {
                if (get<0>((*it)->state)==del)
                {
                    to_visit.push(*it);
                    visited_count.push(0);
                }
            }
        }else if (get<0>(u->state)==del)
        {
            if (c==0)
            {
                to_visit.push(u);
                visited_count.push(1);
                gap.push_back(u);
                for (auto it=u->out.begin(); it!=u->out.end(); ++it)
                {
                    //if (get<0>((*it)->state)==del)
                    {
                        to_visit.push(*it);
                        visited_count.push(0);
                    }
                }
            }else
            {
                gap.pop_back();
            }
        }else
        {
            deletes.push_back(GapEx(w,u,gap));
        }
    }
}

void PartialOrderGraph::find_delete_at_level(int i,vector<GapEx>& deletes)
{
    PartialOrderGraphNode* u = nodes[i];
    find_delete_from(u,deletes);
    for (auto it=u->sibling.begin(); it!=u->sibling.end(); ++it)
    {
        find_delete_from(*it,deletes);
    }
}

void PartialOrderGraph::canonize_delete_at_level(int i)
{
    vector<GapEx> deletes;
    find_delete_at_level(i,deletes);

    int ul,vl,dl,dd,rid,rcn;
    PartialOrderGraphNode *u, *v, *w;
    PartialOrderGraphNodeState s(del,"=");
    map<PartialOrderGraphNode*,int> nl;

    if (deletes.empty())
    {
        return;
    }

    for (auto it=deletes.begin(); it!=deletes.end(); ++it)
    {
        u = get<0>(*it);
        v = get<1>(*it);
        
        auto uit = nl.find(u);
        if (uit==nl.end())
        {
            ul = node_level_exclude_delete(u);
            nl[u] = ul;
        }
        ul = nl[u];

        auto vit = nl.find(v);
        if (vit==nl.end())
        {
            vl = node_level_exclude_delete(v);
            nl[v] = vl;
        }
        vl = nl[v];

        dl = vl-ul-1;
        dd = get<2>(*it).size();

        if (dl-dd>0)
        {
            v = get<2>(*it)[0];
            rid = get<0>(v->read_pool[0]);
            rcn = get<2>(v->read_pool[0]);
            Gap g;
            for (int t=dl-dd; t>0; --t)
            {
                w = new PartialOrderGraphNode(N,s);
                add_node(w);
                w->read_pool.push_back(ReadBase(rid,get<1>(s),rcn));
                g.push_back(w);
            }
            add_edge(u,v,g);
            delete_edge(u,v);
        }
    }
}

void PartialOrderGraph::canonize_delete()
{
    for (int i=0; i<N; i++)
    {
        if (get<1>(nodes[i]->state)=="$")
        {
            break;
        }
        canonize_delete_at_level(i);
    }
}

void PartialOrderGraph::canonize_graph()
{
    // canonize insert
    canonize_insert();

    // canonize delete
    canonize_delete();

    // forward merging
    forward_merge();

    // backward merging
    backward_merge();
}

void PartialOrderGraph::node_level()
{
    auto it = level_begin();
    for (; it!=level_end(); ++it)
    {
        get<0>(*it)->level = get<1>(*it);
    }
}

typedef tuple<int,int> read_t;
typedef set<read_t> read_pool_t;
void find_common_read_pool(PartialOrderGraphNode* a, PartialOrderGraphNode* b, read_pool_t& c)
{
    read_pool_t a_read, b_read;
    
    // read pool of node a
    for (auto it=a->read_pool.begin(); it!=a->read_pool.end(); ++it)
    {
        a_read.emplace(read_t(get<0>(*it),get<2>(*it)));
    }
    for (auto it=a->out.begin(); it!=a->out.end(); ++it)
    {
        if (get<0>((*it)->state)==ins or get<0>((*it)->state)==del)
        {
            for (auto itt=(*it)->read_pool.begin(); itt!=(*it)->read_pool.end(); ++itt)
            {
                auto dit = a_read.find(read_t(get<0>(*itt),get<2>(*itt)));
                a_read.erase(dit);
            }
        }
    }

    // read pool of node b
    for (auto it=b->read_pool.begin(); it!=b->read_pool.end(); ++it)
    {
        b_read.emplace(read_t(get<0>(*it),get<2>(*it)));
    }
    for (auto it=b->in.begin(); it!=b->in.end(); ++it)
    {
        if (get<0>((*it)->state)==ins or get<0>((*it)->state)==del)
        {
            for (auto itt=(*it)->read_pool.begin(); itt!=(*it)->read_pool.end(); ++itt)
            {
                auto dit = b_read.find(read_t(get<0>(*itt),get<2>(*itt)));
                b_read.erase(dit);
            }
        }
    }

    // set intersection
    c.clear();
    for (auto ait=a_read.begin(); ait!=a_read.end(); ++ait)
    {
        auto bit = b_read.find(*ait);
        if (bit!=b_read.end())
        {
            c.emplace(*ait);
        }
    }

}

void PartialOrderGraph::add_edge(int i, int l)
{
    PartialOrderGraphNodeState s(ins,"-");
    PartialOrderGraphNode *u,*v,*w;
    read_pool_t crp;
    
    u = nodes[i];
    v = nodes[i+1];

    if (linking(u,v))
    {
        find_common_read_pool(u,v,crp);

        Gap gap;
        for (int t=0; t<l; ++t)
        {
            w = new PartialOrderGraphNode(N,s);
            add_node(w);
            for (auto crpit=crp.begin(); crpit!=crp.end(); ++crpit)
            {
                w->read_pool.push_back(ReadBase(get<0>(*crpit),get<1>(s),get<1>(*crpit)));
            }
            gap.push_back(w);
        }
        add_edge(u,v,gap);
    }

    for (auto it=v->sibling.begin(); it!=v->sibling.end(); ++it)
    {
        if (linking(u,*it))
        {
            find_common_read_pool(u,*it,crp);

            Gap gap;
            for (int t=0; t<l; ++t)
            {
                w = new PartialOrderGraphNode(N,s);
                add_node(w);
                for (auto crpit=crp.begin(); crpit!=crp.end(); ++crpit)
                {
                    w->read_pool.push_back(ReadBase(get<0>(*crpit),get<1>(s),get<1>(*crpit)));
                }
                gap.push_back(w);
            }
            add_edge(u,*it,gap);
        }
    }

    for (auto it=u->sibling.begin(); it!=u->sibling.end(); ++it)
    {
        if (linking(*it,v))
        {
            find_common_read_pool(*it,v,crp);

            Gap gap;
            for (int t=0; t<l; ++t)
            {
                w = new PartialOrderGraphNode(N,s);
                add_node(w);
                for (auto crpit=crp.begin(); crpit!=crp.end(); ++crpit)
                {
                    w->read_pool.push_back(ReadBase(get<0>(*crpit),get<1>(s),get<1>(*crpit)));
                }
                gap.push_back(w);
            }
            add_edge(*it,v,gap);
        }
    }

    for (auto uit=u->sibling.begin(); uit!=u->sibling.end(); ++uit)
    {
        for (auto vit=v->sibling.begin(); vit!=v->sibling.end(); ++vit)
        {
            if (linking(*uit,*vit))
            {
                find_common_read_pool(*uit,*vit,crp);

                Gap gap;
                for (int t=0; t<l; ++t)
                {
                    w = new PartialOrderGraphNode(N,s);
                    add_node(w);
                    for (auto crpit=crp.begin(); crpit!=crp.end(); ++crpit)
                    {
                        w->read_pool.push_back(ReadBase(get<0>(*crpit),get<1>(s),get<1>(*crpit)));
                    }
                    gap.push_back(w);
                }
                add_edge(*uit,*vit,gap);
            }
        }
    }
}

void PartialOrderGraph::delete_edge(int i)
{
    PartialOrderGraphNode* u = nodes[i];
    PartialOrderGraphNode* v = nodes[i+1];

    if (linking(u,v))
    {
        delete_edge(u,v);
    }

    for (auto it=v->sibling.begin(); it!=v->sibling.end(); ++it)
    {
        if (linking(u,*it))
        {
            delete_edge(u,*it);
        }
    }

    for (auto it=u->sibling.begin(); it!=u->sibling.end(); ++it)
    {
        if (linking(*it,v))
        {
            delete_edge(*it,v);
        }
    }

    for (auto uit=u->sibling.begin(); uit!=u->sibling.end(); ++uit)
    {
        for (auto vit=v->sibling.begin(); vit!=v->sibling.end(); ++vit)
        {
            if (linking(*uit,*vit))
            {
                delete_edge(*uit,*vit);
            }
        }
    }
}

void merge_read_pool(PartialOrderGraphNode* u, PartialOrderGraphNode* v)
{
    int i=0,j=0,m=u->read_pool.size(),n=v->read_pool.size();
    sort(u->read_pool.begin(),u->read_pool.end());
    sort(v->read_pool.begin(),v->read_pool.end());

    vector<ReadBase> result;
    while(i<m and j<n)
    {
        int u_rid = get<0>(u->read_pool[i]);
        int v_rid = get<0>(v->read_pool[j]);
        string u_rb = get<1>(u->read_pool[i]);
        string v_rb = get<1>(v->read_pool[j]);
        int u_rcp = get<2>(u->read_pool[i]);
        int v_rcp = get<2>(v->read_pool[j]);
        if (u_rid==v_rid)
        {
            result.push_back(ReadBase(u_rid,u_rb+v_rb,u_rcp));
            i++;
            j++;
        }else if (u_rid<v_rid)
        {
            result.push_back(ReadBase(u_rid,u_rb,u_rcp));
            i++;
        }else
        {
            result.push_back(ReadBase(v_rid,v_rb,v_rcp));
            j++;
        }
    }
    while (i<m)
    {
        result.push_back(u->read_pool[i++]);
    }
    while (j<n)
    {
        result.push_back(v->read_pool[j++]);
    }

    // copy result to u
    u->read_pool.clear(); u->read_pool.resize(0);
    copy(result.begin(),result.end(),back_inserter(u->read_pool));
}

void PartialOrderGraph::merge_node(PartialOrderGraphNode* u, PartialOrderGraphNode* v)
{
    // add linking from v's parent to u
    for (auto it=v->in.begin(); it!=v->in.end(); ++it)
    {
        if (!linking(*it,u) and (*it)!=u)
        {
            add_edge(*it,u);
        }
    }

    // add linking from u to v's children
    for (auto it=v->out.begin(); it!=v->out.end(); ++it)
    {
        if (!linking(u,*it) and u!=(*it))
        {
            add_edge(u,*it);
        }
    }

    // merge state
    if (linking(u,v) and get<0>(u->state)==mat and get<0>(v->state)==mat)
    {
        get<1>(u->state) += get<1>(v->state);
    }

    // add v's read to u
    merge_read_pool(u,v);

    // delete node v
    delete_node(v,false);
}

void PartialOrderGraph::forward_merge()
{
    PartialOrderGraphNode *u,*v,*w;
    vector<tuple<PartialOrderGraphNode*,
                 PartialOrderGraphNode*>> to_merge;
    set<PartialOrderGraphNode*> merged_node;
    queue<PartialOrderGraphNode*> to_visit;   
    set<PartialOrderGraphNode*> visited;
    to_visit.push(nodes[0]);
    while(!to_visit.empty())
    {
        w = to_visit.front(); to_visit.pop();
        if (merged_node.count(w)>0) continue;
        // find pairs of nodes to merge
        for (auto uit=w->out.begin(); uit!=w->out.end(); ++uit)
        {
            u = *uit;
            for (auto vit=next(uit); vit!=w->out.end(); ++vit)
            {
                v = *vit;
                if (u==v) continue;
                if (u->state==v->state)
                {
                    if (merged_node.count(u)==0 and merged_node.count(v)==0)
                    {
                        to_merge.push_back(make_tuple(u,v));
                        merged_node.emplace(v);
                    }
                }
            }
        }

        // merge
        for (auto it=to_merge.begin(); it!=to_merge.end(); ++it)
        {
            u = get<0>(*it);
            v = get<1>(*it);
            merge_node(u,v);
        }

        // push children of w
        for (auto it=w->out.begin(); it!=w->out.end(); ++it)
        {
            if (visited.count(*it)==0)
            {
                to_visit.push(*it);
                visited.emplace(*it);
            }
        }

        // clean temporary storage
        to_merge.clear(); to_merge.resize(0);
//        merged_node.clear(); 
    }
}


void PartialOrderGraph::backward_merge()
{
    PartialOrderGraphNode *u,*v,*w;
    queue<PartialOrderGraphNode*> to_visit;
    set<PartialOrderGraphNode*> visited;
    vector<tuple<PartialOrderGraphNode*,
                 PartialOrderGraphNode*>> to_merge;
    set<PartialOrderGraphNode*> merged_node;

    for (auto it=nodes.begin(); it!=nodes.end(); ++it)
    {
        if (get<1>((*it)->state)=="$")
        {
            to_visit.push(*it);
        }
    }

    while (!to_visit.empty())
    {
        w = to_visit.front(); to_visit.pop();
        if (merged_node.count(w)>0) continue;
        // find pairs of nodes to merge
        for (auto uit=w->in.begin(); uit!=w->in.end(); ++uit)
        {
            u = *uit;
            for (auto vit=next(uit); vit!=w->in.end(); ++vit)
            {
                v = *vit;
                if (u==v) continue;
                if (u->state==v->state)
                {
                    if (merged_node.count(u)==0 and merged_node.count(v)==0)
                    {
                        to_merge.push_back(make_tuple(u,v));
                        merged_node.emplace(v);
                    }
                }
            }
        }

        // merge
        for (auto it=to_merge.begin(); it!=to_merge.end(); ++it)
        {
            u = get<0>(*it);
            v = get<1>(*it);
            merge_node(u,v);
        }

        // push parent of w to visit
        for (auto it=w->in.begin(); it!=w->in.end(); ++it)
        {
            if (visited.count(*it)==0)
            {
                to_visit.push(*it);
                visited.emplace(*it);
            }
        }

        // clean
        to_merge.clear(); to_merge.resize(0);
        //merged_node.clear();
    }
}

int PartialOrderGraph::indegree(PartialOrderGraphNode* u)
{
    return u->in.size();
}

int PartialOrderGraph::outdegree(PartialOrderGraphNode* u)
{
    return u->out.size();
}

void PartialOrderGraph::path_collapse()
{
    int level = 0,level_size=0;
    PartialOrderGraphNode *u,*v;
    queue<PartialOrderGraphNode*> level_node;
    queue<PartialOrderGraphNode*> sublevel_node;
    set<PartialOrderGraphNode*> multi_in;

    level_node.push(nodes[0]);
    while (!level_node.empty())
    {
        u = level_node.front(); level_node.pop();
        
        if (level_size==1 and outdegree(u)==1)
        {
            v = u->out[0];
            while (outdegree(v)==1)
            {
                merge_node(u,v);
                v = u->out[0];
            }
        }

        for (auto it=u->out.begin(); it!=u->out.end(); ++it)
        {
            v = *it;
            if (multi_in.find(v)==multi_in.end())
            {
                sublevel_node.push(v);
                multi_in.emplace(v);
            }
        }

        if (level_node.empty())
        {
            while (!sublevel_node.empty())
            {
                v = sublevel_node.front(); sublevel_node.pop();
                level_node.push(v);
            }
            level += 1;
            level_size = level_node.size();
            multi_in.clear();
        }
    }
}

int PartialOrderGraph::number_of_reads_cover_nodes(PartialOrderGraphNode* u, PartialOrderGraphNode* v)
{
    int n = 0;
    if (u==nodes[0])
    {
        for (auto vit=v->read_pool.begin(); vit!=v->read_pool.end(); ++vit)
            n += get<2>(*vit);
    }else if (get<1>(v->state)=="$")
    {
        for (auto uit=u->read_pool.begin(); uit!=u->read_pool.end(); ++uit)
            n += get<2>(*uit);
    }else
    {
        for (auto uit=u->read_pool.begin(); uit!=u->read_pool.end(); ++uit)
        {
            for (auto vit=v->read_pool.begin(); vit!=v->read_pool.end(); ++vit)
            {
                if (get<0>(*uit)==get<0>(*vit))
                {
                    n += get<2>(*vit);
                }
            }
        }
    }

    return n;
}

