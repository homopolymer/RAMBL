// Multiple Sequence Alignment Class
// ---------------------------------
// It provides an engine to build and manipulate a
// multiple sequence alignment.
//
//
// Last Changed
// ------------
// July 25, 2015    Feng Zeng    Create it.
//
//
//

#ifndef _MULTIPLE_SEQUENCE_ALIGNMENT_HPP
#define _MULTIPLE_SEQUENCE_ALIGNMENT_HPP

#include "Score.hpp"
#include <vector>
#include <array>
#include <list>
#include <iterator>
using namespace std;

//
// generic class for MSA object
// ----------------------------
//

template<template<class,class> class C, class T>
class MSA
{
public:
    MSA(){}
    ~MSA(){}

public:
    class C<C<T,allocator<T>>,allocator<C<T,allocator<T>>>>::iterator begin() { return align_data.begin(); }
    class C<C<T,allocator<T>>,allocator<C<T,allocator<T>>>>::iterator end()   { return align_data.end(); }

    class C<C<T,allocator<T>>,allocator<C<T,allocator<T>>>>::reverse_iterator rbegin() { return align_data.rbegin(); }
    class C<C<T,allocator<T>>,allocator<C<T,allocator<T>>>>::reverse_iterator rend() { return align_data.rend(); }

    void push_back(const C<T,allocator<T>>& a)
    {
        align_data.push_back(a);
    }

    int size()
    {
        return align_data.size();
    }

    void resize(int z)
    {
        align_data.clear();
        align_data.resize(z);
    }

    void get(int i, C<T,allocator<T>>& seq)
    {
        auto it = begin();
        for (; it!=end(); ++it)
        {
            int ii = 0;
            auto itt = (*it).begin();
            for (; ii<i; ++ii,++itt)
            {}
            if (itt!=(*it).end()) seq.push_back(*itt);
        }
    }

public:
    C<C<T,allocator<T>>,allocator<C<T,allocator<T>>>> align_data;
};


// generic class for alignment engine
template<class Index, class Model, template<class,class> class C, class Data, class T>
class MultipleSequenceAlignment
{
public:
    MultipleSequenceAlignment(){}
    ~MultipleSequenceAlignment(){}
public:
    virtual double align(C<Data,allocator<Data>>& data, MSA<C,T>& msa) = 0;
};


// sum of pairs
template<class Index, class Model, template<class,class> class C, class Data, class T>
class MultipleSequenceAlignmentSP: public MultipleSequenceAlignment<Index,Model,C,Data,T>
{
public:
    MultipleSequenceAlignmentSP()
    { this->score = new SimpleDnaScore; use_default = true; }
    MultipleSequenceAlignmentSP(Score<Index,Model>& s)
    { this->score = &s; use_default = false; }
    ~MultipleSequenceAlignmentSP()
    { if (use_default) delete this->score; }

public:
    double align(C<Data,allocator<Data>>& data, MSA<C,T>& msa);
    void forward(Data& seq, MSA<C,T>& msa, double* SC, int* PI, int* SI, int* SJ, int* PP, int m, int n, int s);
    void backward(Data& seq, MSA<C,T>& msa, int* SI, int* SJ, int m, int n, int s);

public:
    bool use_default;
    Score<Index,Model>* score;
};

#endif
