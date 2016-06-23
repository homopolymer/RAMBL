#include "Score.hpp"
#include <vector>
using namespace std;

SimpleDnaScore::SimpleDnaScore(SimpleScoreModel& model)
{
    set(model);
}

double SimpleDnaScore::operator[] (const Index2D& index)
{
    return score_matrix[index];
}

void SimpleDnaScore::set(SimpleScoreModel& model)
{
    vector<char> alphabet = {'A','a','C','c','G','g','T','t','+','-'};
    for (auto x : alphabet)
    {
        for (auto y : alphabet)
        {
            if (x==y) 
                score_matrix[Index2D(x,y)] = model.match;
            else if ((x=='A' and y=='a') or (x=='a' and y=='A'))
                score_matrix[Index2D(x,y)] = model.match;
            else if ((x=='C' and y=='c') or (x=='c' and y=='C'))
                score_matrix[Index2D(x,y)] = model.match;
            else if ((x=='G' and y=='g') or (x=='g' and y=='G'))
                score_matrix[Index2D(x,y)] = model.match;
            else if ((x=='T' and y=='t') or (x=='t' and y=='T'))
                score_matrix[Index2D(x,y)] = model.match;
            else if ((x=='+' and y=='-') or (x=='-' and y=='+'))
                score_matrix[Index2D(x,y)] = model.match;
            else if (x=='+' or y=='+')
                score_matrix[Index2D(x,y)] = model.gap_open+model.gap_extend;
            else if (x=='-' or y=='-')
                score_matrix[Index2D(x,y)] = model.gap_extend;
            else
                score_matrix[Index2D(x,y)] = model.mismatch;
        }
    }
}
