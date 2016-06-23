#include "PartialOrderGraph.hpp"

LevelOrderIterator::LevelOrderIterator(PartialOrderGraphNode* begin, int n)
{
    this->n = n;
    this->level = 0;
    this->iter = PartialOrderGraphNodeL(begin,0);

    for (auto it=begin->out.begin(); it!=begin->out.end(); ++it)
    {
        level_node.push(*it);
    }
    visited_node.emplace(begin);
}

LevelOrderIterator& LevelOrderIterator::operator++()
{
    if (sublevel_node.empty()) 
    {
        this->level += 1;
        visited_node.clear();
    }

    PartialOrderGraphNode* w;
    if (!level_node.empty()){
        while (true)
        {
            w = level_node.top(); level_node.pop();
            this->iter = PartialOrderGraphNodeL(w,level);
            this->n += 1;

            for (auto it=w->out.begin(); it!=w->out.end(); ++it)
            {
                sublevel_node.push(*it);
            }       
  
            if (level_node.empty())
            {
                while (!sublevel_node.empty())
                {
                    w = sublevel_node.top(); sublevel_node.pop();
                    if (visited_node.find(w)!=visited_node.end())
                    {
                        continue;
                    }
                    level_node.push(w);
                    visited_node.emplace(w);
                }
            }
            break;
        }
    }else
    {
        this->n += 1;
    }
}
