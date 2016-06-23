// Scoring Function Class
// ----------------------
// It provides an engine to build and manipulate a
// scoring function for sequence alignment.
//
//
// Last Changed
// ------------
// July 25, 2015    Feng Zeng    Create it.
//

#ifndef _SCORE_HPP
#define _SCORE_HPP

#include <map>
#include <tuple>
using namespace std;

// abstract score class
template<class Index, class Model>
class Score
{
public:
    Score(){}
    ~Score(){}
public:
    virtual double operator[] (const Index& index) = 0;
    virtual void set(Model& model) = 0;
};

// simple score class
typedef tuple<char,char> Index2D;
struct SimpleScoreModel
{
    SimpleScoreModel(double m=3,double x=-5,double o=-4,double e=-2)
    {
        match      = m;
        mismatch   = x;
        gap_open   = o;
        gap_extend = e;
    }
    double match;
    double mismatch;
    double gap_open;
    double gap_extend;
};

class SimpleDnaScore: public Score<Index2D, SimpleScoreModel>
{
public:
    SimpleDnaScore()
    {
        SimpleScoreModel model;
        set(model);
    }
    SimpleDnaScore(SimpleScoreModel& model);
    ~SimpleDnaScore(){}

public:
    double operator[] (const Index2D& index);
    void set(SimpleScoreModel& model);

public:
    map<Index2D,double> score_matrix;
};

#endif
