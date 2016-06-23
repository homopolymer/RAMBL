#include "PartialOrderGraph.hpp"

TopoSortIterator::TopoSortIterator(PartialOrderGraphNode* begin, int n)
{
    n = n;

    // postorder traversal
    PartialOrderGraphNode* u;
    int c;
    to_visit.push(begin);
    visited_count.push(0);
    while(!to_visit.empty())
    {
        u = to_visit.top();
        c = visited_count.top();

        if (u->out.empty())
        {
            toposort.push(u);
            to_visit.pop();
            visited_count.pop();
        }else
        {
            if (c==1)
            {
                toposort.push(u);
                to_visit.pop();
                visited_count.pop();
            }else if (c<1)
            {
                visited_count.pop();
                visited_count.push(c+1);
                for (auto it=u->out.begin(); it!=u->out.end(); ++it)
                {
                    to_visit.push(*it);
                    visited_count.push(0);
                }
            }
        }
    }

    iter = toposort.top();
    toposort.pop();
    n += 1;
}

TopoSortIterator& TopoSortIterator::operator++()
{
    iter = toposort.top();
    toposort.pop();
    n += 1;
}
