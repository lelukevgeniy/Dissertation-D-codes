#ifndef FACTORCLASSES_H
#define FACTORCLASSES_H

#include "SBitsArray.h"
#include "./vectorspacef2.h"

class FactorClasses
{
public:
    FactorClasses();
    FactorClasses(uint8_t dimension);
    ~FactorClasses();
    void init(uint8_t dimension);
    SBitsArray Proj(SBitsArray &y, SBitsArray& Bi);
    SBitsArray Proj_fast(SBitsArray &y, SBitsArray& Bi);
    SBitsArray* get_factor_classes_for_Bi(SBitsArray& Bi);
    uint16_t* get_factor_classes_idx_for_Bi(SBitsArray& Bi);

    void make_factor_classes();
    void print_factor_classes();
    void print_factor_classes_idx();
    void print_factor_class_for_Bi(SBitsArray& Bi);

private:
    VectorSpaceF2* space;
    SBitsArray** factor_classes;
    uint16_t** factor_classes_idx;
    uint8_t dim;
};

#endif // FACTORCLASSES_H
