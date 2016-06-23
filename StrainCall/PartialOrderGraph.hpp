// Partial Order Graph Class
// -------------------------
// It provides an engine to construct and manipulate
// a partial order graph.
//
// It provides the following features:
//     (a) iterator for topological order traversal
//     (b) iterator for level order traversal
//
//
//
//
// Examples
// --------
// [1] Create a partial order graph
//
//     #include "PartialOrderGraph.hpp"
//     #include <iostream>
//
//     GenomeSeq G = "ACGTACGT";
//
//     vector<AlignRead> R;
//     R.push_back(AlignRead(0,"8M","ACCTACGT"));
//     R.push_back(AlignRead(0,"3M1I5M","ACCCTACGT"));
//     R.push_back(AlignRead(0,"3M1D4M","ACCACGT"));
//
//     PartialOrderGraph* pog = new PartialOrderGraph(G,R);
//
//     pog->output_edge(std::cout);
//
//
//
//
// TODO
// ----
// [1] Partial order graph construction
//     (a) Insert canonizer                   [DONE]
//     (b) Graph canonizer                    [DONE]
//     (c) Vertex-read matrix                 [DONE]
//     (d) path collapse                      [DONE]
// [2] Topological-order visiting             [DONE]
// [3] Level-order visiting                   [DONE]
// [4] Streaming clustering                   [DONE]
// [5] graph based error correction
// [6] sequence to graph alignment
// 
// 
//
//
// Last Changed
// ------------
// July 14, 2015    Feng Zeng    Create it.
// July 15, 2015    Feng Zeng    Add iterator for topological sorting.
// July 23, 2015    Feng Zeng    Add extended gap type.
// July 27, 2015    Feng Zeng    Add deletion canonization method.
//                               Add distance of two nodes.
//                               Add iterator for level visiting.
// July 28, 2015    Feng Zeng    Add forward merging and backward merging.
// July 29, 2015    Feng Zeng    Add node indegree and outdegree.
//                               Add path collapsing method.
// August 1, 2015   Feng Zeng    Add Strain class.                               
//                               Add nonparametric bayesian clustering
//
//

#ifndef _PARTIALORDERGRAPH_HPP
#define _PARTIALORDERGRAPH_HPP

#include <string>
#include <vector>
#include <set>
#include <map>
#include <tuple>
#include <iterator>
#include <iostream>
#include <cstddef>
#include <stack>
#include <queue>
using namespace std;


typedef enum {mat=0,mis,ins,del} AlignState;
typedef tuple<AlignState,string> State;
typedef State PartialOrderGraphEdgeState;
typedef State PartialOrderGraphNodeState;

typedef tuple<int,string,int> ReadBase;
typedef map<int,vector<int>> ReadPairs;
// Partial Order Graph Node Class
class PartialOrderGraphNode
{
public:
    PartialOrderGraphNode(){}
    PartialOrderGraphNode(int id, PartialOrderGraphNodeState state);
    ~PartialOrderGraphNode(){}
    
public:
    typedef vector<PartialOrderGraphNode*>::iterator in_iterator;
    typedef vector<PartialOrderGraphNode*>::iterator out_iterator;
    typedef vector<PartialOrderGraphNode*>::iterator sibling_iterator;

public:
    void insert_in(PartialOrderGraphNode* node);
    void delete_in(PartialOrderGraphNode* node);

    void insert_out(PartialOrderGraphNode* node);
    void delete_out(PartialOrderGraphNode* node);

    void insert_sibling(PartialOrderGraphNode* node);

    out_iterator find_out(PartialOrderGraphNodeState out_state);
    sibling_iterator find_sibling(PartialOrderGraphNodeState sibling_state);
    

public:
    int id;
    int level;
    int pos;
    PartialOrderGraphNodeState state;
    vector<PartialOrderGraphNode*> in;
    vector<PartialOrderGraphNode*> out;
    vector<PartialOrderGraphNode*> sibling;
    vector<ReadBase> read_pool;
};


typedef vector<PartialOrderGraphNode*> PartialOrderGraphEdge;
typedef vector<PartialOrderGraphNode*> Gap;
typedef tuple<PartialOrderGraphNode*, PartialOrderGraphNode*, Gap> GapEx;


// iterator for topological sorting
class TopoSortIterator:public iterator<forward_iterator_tag,
                                       PartialOrderGraphNode,
                                       ptrdiff_t,
                                       PartialOrderGraphNode*,
                                       PartialOrderGraphNode&>
{
public:
    PartialOrderGraphNode* iter;
    stack<PartialOrderGraphNode*> toposort;
    stack<PartialOrderGraphNode*> to_visit;
    stack<int> visited_count;
    int n;

public:
    explicit TopoSortIterator(PartialOrderGraphNode* begin, int n);

public:
    TopoSortIterator& operator++();

    PartialOrderGraphNode& operator*() const
    {
        return *iter;
    }

    PartialOrderGraphNode* operator->() const
    {
        return iter;
    }

    bool operator!= (const TopoSortIterator& rhs) const
    {
        return n != rhs.n;
    }

    bool operator< (const TopoSortIterator& rhs) const
    {
        return n < rhs.n;
    }
};

// iterator for level sorting
typedef tuple<PartialOrderGraphNode*, int> PartialOrderGraphNodeL;
class LevelOrderIterator:public iterator<forward_iterator_tag,
                                         PartialOrderGraphNode,
                                         ptrdiff_t,
                                         PartialOrderGraphNode*,
                                         PartialOrderGraphNode&>
{
public:
    PartialOrderGraphNodeL iter;
    stack<PartialOrderGraphNode*> level_node;
    stack<PartialOrderGraphNode*> sublevel_node;
    set<PartialOrderGraphNode*> visited_node;
    int n;
    int level;
public:
    explicit LevelOrderIterator(PartialOrderGraphNode* begin, int n);
public:
    LevelOrderIterator& operator++();
    PartialOrderGraphNodeL& operator*() 
    {
        return iter;
    }
    PartialOrderGraphNodeL* operator->() 
    {
        return &iter;
    }
    bool operator!=(const LevelOrderIterator& rhs) const
    {
        return n != rhs.n;
    }
    bool operator<(const LevelOrderIterator& rhs) const
    {
        return n < rhs.n;
    }
};

// Partial Order Graph Class
typedef string GenomeSeq;
typedef string ReadSeq;
typedef int GenomePosition;
typedef int RelativePosition;
typedef int CopyNumber;
typedef string AlignCigar;
typedef string ReadQual;
typedef tuple<RelativePosition,AlignCigar,ReadSeq,ReadQual,CopyNumber> AlignRead;
typedef tuple<int,string> ReadAux; // <id,read_name>

typedef string CigarOp;
typedef int CigarOpLen;
typedef tuple<CigarOp,CigarOpLen> CigarRecord;
void parse_cigar(AlignCigar& cigar, vector<CigarRecord>& records);

typedef vector<PartialOrderGraphNode*> Path;

class Strain;
typedef long double DoubleL;

class PartialOrderGraph
{
public:
    PartialOrderGraph(){N=0;}
    // TODO to specify
    PartialOrderGraph(GenomeSeq& G, vector<AlignRead>& R); 
    ~PartialOrderGraph()
    {
        for (auto it=nodes.begin(); it!=nodes.end(); ++it)
        {
            if (*it) delete (*it);
        }
        for (auto it=deleted_nodes.begin(); it!=deleted_nodes.end(); ++it)
        {
            if (*it) delete (*it);
        }
    }

public:
    // builder
    void build(GenomeSeq& G, vector<AlignRead>& R);

    // add a node
    void add_node(PartialOrderGraphNode* w);
    // delete a node
    void delete_node(PartialOrderGraphNode* w, bool bridging=true);
    // add a node w in between u and v
    void add_edge(PartialOrderGraphNode* u, PartialOrderGraphNode* w);
    // add a gap after u
    void add_edge(PartialOrderGraphNode* u, Gap& gap);
    // add a gap, that could be insert or delete, in between u and v
    void add_edge(PartialOrderGraphNode* u, PartialOrderGraphNode* v, Gap& gap);
    // add a gap between level i and i+1
    void add_edge(int i, int l);
    // delete the edge linking u and v
    void delete_edge(PartialOrderGraphNode* u, PartialOrderGraphNode* v);
    // delete edge between level i and i+1
    void delete_edge(int i);
    // check whether there is a linking u and v
    bool linking(PartialOrderGraphNode* u, PartialOrderGraphNode* v);

    // find the insertions starting from u 
    void find_insert_from(PartialOrderGraphNode* u, vector<GapEx>& inserts);
    // find the insertions starting at level i
    void find_insert_at_level(int i, vector<GapEx>& inserts);
    // canonize the insertions starting at level i
    void canonize_insert_at_level(int i);
    // canonize the insertions on graph
    void canonize_insert();

    // find the deletions starting from u
    void find_delete_from(PartialOrderGraphNode* u, vector<GapEx>& deletes);
    // find the deletions starting at level i
    void find_delete_at_level(int i, vector<GapEx>& deletes);
    // canonize the deletions starting at level i
    void canonize_delete_at_level(int i);
    // canonize the deletion on graph
    void canonize_delete();

    // merge node
    void merge_node(PartialOrderGraphNode* u, PartialOrderGraphNode* v);
    // forward merging
    void forward_merge();
    // backward merging
    void backward_merge();
    // canonize the partial order graph
    void canonize_graph();
    
    // get node level for non-delete nodes
    int node_level_exclude_delete(PartialOrderGraphNode* w);
    // get node level for all kinds of nodes
    int node_level(PartialOrderGraphNode* w);
    // set node level
    void node_level();
    // compute the distance between two nodes assuming these two nodes are linked by a path
    int distance(PartialOrderGraphNode* u, PartialOrderGraphNode* v);

    // output edge list
    void output_edge(ostream& f);

    // in degree
    int indegree(PartialOrderGraphNode* u);
    // out degree
    int outdegree(PartialOrderGraphNode* u);
    // collapse simple path
    void path_collapse();

    // reads cover u and v
    int number_of_reads_cover_nodes(PartialOrderGraphNode* u, PartialOrderGraphNode* v);

private:
    // nonparametric bayesian clustering
    void np_bayes_clustering(vector<Strain>& strains, 
                             vector<ReadBase>& reads, ReadPairs& read_pairs, int n, 
                             vector<DoubleL>& abundance,
                             vector<DoubleL>& posterior);
    // hard clustering
    void hard_clustering(vector<Strain>& strains, vector<ReadBase>& reads, 
                         vector<int>& new_reads, ReadPairs& read_pairs);
    // streaming clustering
    void streaming_clustering(vector<Strain>& strains, ReadPairs& read_pairs,
                              int n, DoubleL e, DoubleL tau, DoubleL diff);
public:
    // strain calling
    void infer_strains(vector<Strain>& strains, ReadPairs& read_pairs, 
                       int n, DoubleL e, DoubleL tau, DoubleL diff);

    // read assignment
    void read_assign(vector<Strain>& strains,vector<ReadAux>& reads,int n);
    void read_assign(vector<Strain>& strains,vector<AlignRead>& reads,ReadPairs& read_pairs,int n);

public:
    // iterator
    typedef TopoSortIterator topo_iterator;
    topo_iterator topo_begin() { return topo_iterator(nodes[0],0); }
    topo_iterator topo_end() { return topo_iterator(nodes[0],N); }

    typedef LevelOrderIterator level_iterator;
    level_iterator level_begin() { return level_iterator(nodes[0],0); }
    level_iterator level_end() { return level_iterator(nodes[0],N); }
   
public:
    int N;
    vector<PartialOrderGraphNode*> nodes;
    int M;
    vector<PartialOrderGraphNode*> deleted_nodes;
};

typedef map<string,DoubleL> CompositionCount;
typedef tuple<string,string> Substitution;
typedef map<Substitution,DoubleL> SubstitutionCount;
class Strain
{
public:
    Strain();
    Strain(int N,DoubleL e);
    Strain(const Strain& other);
    ~Strain(){}

public:
    void update_model(DoubleL al,SubstitutionCount& sc);
    void update_read_loglik(int id, DoubleL loglik);
    void update_read_info(int id, SubstitutionCount& c);
    void erase_read(int id);
    void clean_read(vector<ReadBase>& current_reads);
    void path_extend(PartialOrderGraphNode* n);
    string strain_seq();
    string plain_seq();

public:
    DoubleL operator()(int id);
    DoubleL operator()(string a);           // probability
    DoubleL operator()(string a,string b);  // probability
    DoubleL logprob(int id);
    DoubleL logprob(string a);              // log probability
    DoubleL logprob(string a, string b);    // log probability

public:
    Strain& operator=(const Strain& other);

public:
    DoubleL Z;
    CompositionCount comp_count; // composition model
    SubstitutionCount sub_count; // substitution model

    DoubleL abundance;
    map<int,DoubleL> read_loglik;
    map<int,SubstitutionCount> read_info;
    vector<tuple<string,DoubleL>> assign_reads;
    
    Path path;
};
#endif
