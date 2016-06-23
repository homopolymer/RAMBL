#include "Score.hpp"
#include "MultipleSequenceAlignment.hpp"
#include <array>
#include <vector>
#include <string>
#include <iostream>
using namespace std;


template<class Index, class Model, template<class,class> class C, class Data, class T>
double MultipleSequenceAlignmentSP<Index,Model,C,Data,T>::align(C<Data,allocator<Data>>& data, MSA<C,T>& msa)
{
    int m,n,s;
    
    // insert first element of data to be msa
    auto it = data.begin();   
    for (auto iit=(*it).begin(); iit!=(*it).end(); ++iit)
    {
        msa.push_back(C<T,allocator<T>>(1,*iit));
    }

    // start scan from 2nd data
    ++it;
    ++s;
    for (; it!=data.end(); ++it,++s)
    {
        m = msa.size()+1;
        n = (*it).size()+1;

        double *SC = new double[m*n]; // score of optimal alignment
        int    *PI = new int[m*n];    // state of optimal alignment
        int    *SI = new int[m*n];    // i-shift of optimal alignment
        int    *SJ = new int[m*n];    // j-shift of optimal alignment
        int    *PP = new int[m*n*s];  // separate state of optimal alignment

        // forward computation
        forward(*it,msa,SC,PI,SI,SJ,PP,m,n,s);

        // backward computation
        backward(*it,msa,SI,SJ,m,n,s);

        delete SC;
        delete PI;
        delete SI;
        delete SJ;
        delete PP;
    }
    return 0.0;
}


template<class Index, class Model, template<class,class> class C, class Data, class T>
void MultipleSequenceAlignmentSP<Index,Model,C,Data,T>::forward(Data& seq, MSA<C,T>& msa,
                                                                double *SC,
                                                                int    *PI,
                                                                int    *SI, int *SJ,
                                                                int    *PP,
                                                                int m, int n, int s)
{
    int mat=0, ins=1, del=2;
    int i,j,k;
    double sp;
    Index2D ind;

    // i=0,j=0 
    SC[0] = 0;
    PI[0] = mat;
    SI[0] = 0;
    SJ[0] = 0;
    for (k=0; k<s; k++)
    {
        PP[k] = mat;
    } 

    // i=0,j>0
    for (i=0,j=1; j<n; ++j)
    {
        if (j==1)
        {
            sp = 0;
            for (k=0; k<s; ++k)
            {
                sp += (*score)[Index2D('A','+')];
            }
        }else
        {
            sp = 0;
            for (k=0; k<s; ++k)
            {
                sp += (*score)[Index2D('A','-')];
            }
        }
        sp += SC[i*n+j-1];
        SC[i*n+j] = sp;
        PI[i*n+j] = ins;
        SI[i*n+j] = 0;
        SJ[i*n+j] = -1;

        for (k=0; k<s; ++k)
        {
            PP[i*n*s+j*s+k] = ins;
        }
    }    

    // i>0,j=0
    auto it1 = msa.begin();
    for (i=1,j=0; i<m; ++i,++it1)
    {
        if (i==1)
        {
            sp = 0;
            auto it3 = (*it1).begin();
            for (k=0; k<s; ++k,++it3)
            {
                sp += (*score)[Index2D(*it3,'+')];
            }
        }else
        {
            sp = 0;
            auto it3 = (*it1).begin();
            for (k=0; k<s; ++k,++it3)
            {
                sp += (*score)[Index2D(*it3,'-')];
            }
        }
        sp += SC[(i-1)*n+j];
        SC[i*n+j] = sp;
        PI[i*n+j] = del;
        SI[i*n+j] = -1;
        SJ[i*n+j] = 0;
        auto it3 = (*it1).begin();
        for (k=0; k<s; ++k,++it3)
        {
            if (*it3=='-') PP[i*n*s+j*s+k] = mat;
            else PP[i*n*s+j*s+k] = del;
        }
    }

    // i>0,j>0
    i=1;
    for (auto it1=msa.begin(); it1!=msa.end(); ++it1,++i)
    {
        j=1;
        for (auto it2=seq.begin(); it2!=seq.end(); ++it2,++j)
        {
            double r1,r2,r3;       
            // match
            r1 = 0;
            k = 0;
            for (auto it3=(*it1).begin(); k<s; ++k,++it3)
            {
                if (*it3=='-')
                {
                    if (PP[(i-1)*n*s+(j-1)*s+k]==ins)
                        r1 += (*score)[Index2D('-',*it2)];
                    else
                        r1 += (*score)[Index2D('+',*it2)];
                }else
                {
                    r1 += (*score)[Index2D(*it3,*it2)];
                }
            }
            r1 += SC[(i-1)*n+(j-1)];
           
            // insert
            r2 = 0;
            k = 0;
            for (auto it3=(*it1).begin(); k<s; ++k,++it3)
            {
                if (PP[i*n*s+(j-1)*s+k]==ins)
                {
                    r2 += (*score)[Index2D('-',*it2)];
                }else
                {
                    r2 += (*score)[Index2D('+',*it2)];
                }
            }
            r2 += SC[i*n+(j-1)];

            // delete
            r3 = 0;
            k = 0;
            for (auto it3=(*it1).begin(); k<s; ++k,++it3)
            {
                if (*it3!='-')
                {
                    if (PP[(i-1)*n*s+j*s+k]==del)
                    {
                       r3 += (*score)[Index2D(*it3,'-')];
                    }else
                    {
                       r3 += (*score)[Index2D(*it3,'+')];
                    }
                }else
                {
                    r3 += (*score)[Index2D(*it3,'-')];
                }
            }
            r3 += SC[(i-1)*n+j];

            // maximize
            if (r1>=r2 and r1>=r3)
            {
                SC[i*n+j] = r1;
                PI[i*n+j] = mat;
                SI[i*n+j] = -1;
                SJ[i*n+j] = -1;
                auto it3 = (*it1).begin();
                for (k=0; k<s; ++k)
                {
                    if (*it3=='-')
                    {
                        PP[i*n*s+j*s+k] = ins;
                    }else
                    {
                        PP[i*n*s+j*s+k] = mat;
                    }
                }
            }else if (r2>=r1 and r2>=r3)
            {
                SC[i*n+j] = r2;
                PI[i*n+j] = ins;
                SI[i*n+j] = 0;
                SJ[i*n+j] = -1;
                for (k=0; k<s; ++k)
                {
                    PP[i*n*s+j*s+k] = ins;
                }
            }else
            {
                SC[i*n+j] = r3;
                PI[i*n+j] = del;
                SI[i*n+j] = -1;
                SJ[i*n+j] = 0;
                auto it3 = (*it1).begin();
                for (k=0; k<s; ++k)
                {
                    if (*it3=='-')
                    {
                        PP[i*n*s+j*s+k] = mat;
                    }else
                    {
                        PP[i*n*s+j*s+k] = del;
                    }
                }
            }
        }
    }
}


template<class Index, class Model, template<class,class> class C, class Data, class T>
void MultipleSequenceAlignmentSP<Index,Model,C,Data,T>::backward(Data& seq, MSA<C,T>& msa,
                                                                 int *SI, int *SJ, int m, int n, int s)
{
    MSA<C,T> msa2;
    auto rit1 = msa.rbegin();
    auto rit2 = seq.rbegin();
    int i,j,x,y;
    for (x=m-1; x>=0; )
    {
        for (y=n-1; y>=0; )
        {
            if (x==0 and y==0)
            {
                x -= 1;
                y -= 1;
                continue;
            }

            i = SI[x*n+y];
            j = SJ[x*n+y];

            if (i==-1 and j==-1)
            {
                C<T,allocator<T>> tmp(*rit1);
                tmp.push_back(T(*rit2));
                msa2.push_back(tmp);
                ++rit1;
                ++rit2;
            }else if (i==0 and j==-1)
            {
                C<T,allocator<T>> tmp(s,T('-'));
                tmp.push_back(T(*rit2));
                msa2.push_back(tmp);
                ++rit2;
            }else if (i==-1 and j==0)
            {
                C<T,allocator<T>> tmp(*rit1);
                tmp.push_back(T('-'));
                msa2.push_back(tmp);
                ++rit1;
            }
            x += i;
            y += j;
        }
    }

    msa.resize(msa2.size());
    copy(msa2.rbegin(),msa2.rend(),msa.begin());
}
// explicit template instantiation
template class MultipleSequenceAlignmentSP<Index2D,SimpleScoreModel,vector,string,char>;
template class MultipleSequenceAlignmentSP<Index2D,SimpleScoreModel,list,string,char>;

