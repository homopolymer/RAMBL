#include "PartialOrderGraph.hpp"

PartialOrderGraphNode::PartialOrderGraphNode(int id, PartialOrderGraphNodeState state)
{
    this->id = id;
    this->state = state;
    this->level = -1;
    this->pos = -1;
}

void PartialOrderGraphNode::insert_in(PartialOrderGraphNode* node)
{
   in.push_back(node);
}

void PartialOrderGraphNode::delete_in(PartialOrderGraphNode* node)
{
    in_iterator it = in.begin();
    for (;it!=in.end();++it)
    {
        if (*it==node)
        {
            in.erase(it);
            break;
        }
    }
    
}

void PartialOrderGraphNode::insert_out(PartialOrderGraphNode* node)
{
    out.push_back(node);
}

void PartialOrderGraphNode::delete_out(PartialOrderGraphNode* node)
{
    out_iterator it = out.begin();
    for (;it!=out.end();++it)
    {
        if (*it==node)
        {
            out.erase(it);
            break;
        }
    }
}


void PartialOrderGraphNode::insert_sibling(PartialOrderGraphNode* node)
{
    sibling.push_back(node);
}

PartialOrderGraphNode::out_iterator PartialOrderGraphNode::find_out(PartialOrderGraphNodeState out_state)
{
    auto it = out.begin();
    for (;it!=out.end();++it)
    {
        if ((*it)->state==out_state)
            return it;
    }
    
    return out.end();
}

PartialOrderGraphNode::sibling_iterator PartialOrderGraphNode::find_sibling(PartialOrderGraphNodeState sibling_state)
{
    auto it = sibling.begin();
    for (;it!=sibling.end();++it)
    {
        if ((*it)->state==sibling_state)
            return it;
    }
    return sibling.end();
}
