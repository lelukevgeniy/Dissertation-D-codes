#ifndef VECTORSPACEF2_H
#define VECTORSPACEF2_H

#include "../Codes/SBitsArray.h"
#include "../Codes/combinatorics.h"

typedef struct{
    int list_size;
    int current;
    int * list;
} TLIST;

class VectorSpaceF2
{
public:
    VectorSpaceF2(VectorSpaceF2& vs_for_init)
    {
        space = vs_for_init.space;
        dimension = vs_for_init.dimension;
        weights = vs_for_init.weights;
    };

    VectorSpaceF2(quint64 dimension_v)
        :space(nullptr),
          dimension(0)
    {
        if( (dimension_v > 20) || (dimension_v == 0) )
            cout << "Too big dimesion" << endl;
        else{
            dimension = dimension_v;
            weights = new TLIST[dimension + 1];

            for(int i = 0; i < dimension + 1; i++)
            {
                weights[i].list_size = comb_choose(dimension, i);
                weights[i].list = new int[weights[i].list_size];
                weights[i].current = 0;
            }

            space = new SBitsArray[(int)pow(2,dimension)];
            if(space)
            {
                for(int i = 0; i < (int)pow(2,dimension); i++)
                {
                    space[i].setFromByteArray((unsigned char*)(&i), dimension);
                    int w = space[i].getHammingWeight();
                    weights[w].list[weights[w].current] = i;
                    weights[w].current++;
                }
            }else{
                cout << "error: there is not enough memory" << endl;
            }
        }
    }

    ~VectorSpaceF2()
    {
        if(space)
        {
            for(int i = 0; i < dimension + 1; i++)
                delete [] weights[i].list;
            delete [] weights;
            delete [] space;
        }
    }

    void print_vectors_for_weight(int w)
    {
        if( (w >= 0) && (w <= dimension) )
        {
            for(int i = 0; i < weights[w].list_size; i++)
                cout << space[weights[w].list[i]] << endl;
        }
    }

    SBitsArray* get_first_vector_for_weight(int w)
    {
        if( (w >= 0) && (w <= dimension))
        {
            weights[w].current = 0;
            return &(space[weights[w].list[0]]);
        }else{
            return nullptr;
        }
    };

    SBitsArray* get_next_vector_for_weight(int w)
    {
        if( (w >= 0) && (w <= dimension))
        {
            if(weights[w].current + 1 < weights[w].list_size)
            {
                weights[w].current++;
                return &(space[weights[w].list[weights[w].current]]);
            }else{
                return nullptr;
            }
        }else{
            return nullptr;
        }
    };

    SBitsArray& operator[] (const int index)
    {
        return space[index];
    }

private:
    SBitsArray* space;
    quint64 dimension;
    TLIST * weights;
};

#endif // VECTORSPACEF2_H
