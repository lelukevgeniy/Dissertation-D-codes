#include "factorclasses.h"

FactorClasses::FactorClasses()
    :space(nullptr),
      factor_classes(nullptr),
      factor_classes_idx(nullptr),
      dim(0)
{
}

FactorClasses::FactorClasses(uint8_t dimension)
    :FactorClasses()
{
    init(dimension);
}

void FactorClasses::init(uint8_t dimension)
{
    dim = dimension;
    space = new VectorSpaceF2(dimension);
}

FactorClasses::~FactorClasses()
{
    if(space)
        delete space;

    if(factor_classes)
    {
        for(uint16_t i = 0; i < ((1 << dim) - 1); i++)
        {
            if(factor_classes[i])
                delete [] factor_classes[i];
        }
        delete [] factor_classes;
    }

    if(factor_classes_idx)
    {
        for(uint16_t i = 0; i < ((1 << dim) - 1); i++)
        {
            if(factor_classes_idx[i])
                delete [] factor_classes_idx[i];
        }
        delete [] factor_classes_idx;
    }
}

void FactorClasses::make_factor_classes()
{
    SBitsArray space_mask(1<<dim);

    // память под хранение фактор-классов по всем одномерным подпространствам
    factor_classes = new SBitsArray*[(1 << dim) - 1];
    factor_classes_idx = new uint16_t*[(1 << dim) - 1];

    for(uint16_t j = 0; j < (1 << dim) - 1; j++)
    {
        SBitsArray base_space = (*space)[j + 1];
        //cout << "base_space " << base_space << endl;
        // память под хранение фактор-классов по дному одномерному подпространству
        factor_classes[j] = new SBitsArray[1 << dim];
        factor_classes_idx[j] = new uint16_t[1 << dim];
        uint16_t first_vector_of_factor_class = 0;
        space_mask.toZeroAllBits();

        for(uint16_t current_factor_class_idx = 0; current_factor_class_idx < (1 << (dim -1)); current_factor_class_idx++)
        {
            uint16_t second_vector_of_factor_class = 0;
            // находим первый вектор из смежного класса
            while(space_mask.getBit(first_vector_of_factor_class))
                first_vector_of_factor_class++;

            space_mask.setBitValue(first_vector_of_factor_class, true); // помечаем его занятым

            // находоим второй вектор из смежного класса
            for(second_vector_of_factor_class = 0; second_vector_of_factor_class < (1 << dim); second_vector_of_factor_class++)
            {
                SBitsArray temp = (*space)[second_vector_of_factor_class];
                temp.xorWithSBitArray_fast((*space)[first_vector_of_factor_class]);
                temp.xorWithSBitArray_fast(base_space);
                if(temp.getHammingWeight_fast() == 0)
                    break;
            }
            space_mask.setBitValue(second_vector_of_factor_class, true); // помечаем его занятым
            factor_classes[j][2*current_factor_class_idx].allocateBufferForBits(dim);
            factor_classes[j][2*current_factor_class_idx] = (*space)[first_vector_of_factor_class];
            factor_classes_idx[j][2*current_factor_class_idx] = (*space)[first_vector_of_factor_class].getLowestBitsAsInteger(dim);
            factor_classes[j][2*current_factor_class_idx + 1].allocateBufferForBits(dim);
            factor_classes[j][2*current_factor_class_idx + 1] = (*space)[second_vector_of_factor_class];
            factor_classes_idx[j][2*current_factor_class_idx + 1] = (*space)[second_vector_of_factor_class].getLowestBitsAsInteger(dim);
            //cout << factor_classes[j][2*current_factor_class_idx] << " " << factor_classes[j][2*current_factor_class_idx + 1] <<  endl;
        }
    }
}

SBitsArray FactorClasses::Proj(SBitsArray &y, SBitsArray& Bi)
{
    SBitsArray result(y.getSizeInBits()/2);
    uint16_t Bi_idx = Bi.getLowestBitsAsInteger(dim) - 1;
    uint16_t len = result.getSizeInBits();

    for(uint16_t i = 0; i < len; i++)
    {
        result.setBitValue_fast(i, y.getBit_fast(factor_classes_idx[Bi_idx][2*i]) ^ y.getBit_fast(factor_classes_idx[Bi_idx][2*i + 1]));
    }
    return result;
}

SBitsArray FactorClasses::Proj_fast(SBitsArray &y, SBitsArray& Bi)
{
    SBitsArray result(y.getSizeInBits()/2);
    uint16_t Bi_idx = Bi.getLowestBitsAsInteger_fast(dim) - 1;
    uint16_t len = result.getSizeInBits();

    for(uint16_t i = 0; i < len; i++)
    {
        result.setBitValue_fast(i, y.getBit_fast(factor_classes_idx[Bi_idx][2*i]) ^ y.getBit_fast(factor_classes_idx[Bi_idx][2*i + 1]));
    }
    return result;
}

void FactorClasses::print_factor_classes()
{
    if(factor_classes)
    {
        for(uint16_t i = 0; i < ((1 << dim) - 1); i++)
        {
            cout << endl << "factor class " << i << endl;
            for(uint16_t j = 0; j < (1 << dim); j++)
                cout << factor_classes[i][j] << endl;
        };
    }
}

void FactorClasses::print_factor_classes_idx()
{
    if(factor_classes)
    {
        for(uint16_t i = 0; i < ((1 << dim) - 1); i++)
        {
            cout << endl << "factor class " << i << endl;
            for(uint16_t j = 0; j < (1 << dim); j++)
                cout << factor_classes_idx[i][j] << endl;
        };
    }
}

void FactorClasses::print_factor_class_for_Bi(SBitsArray& Bi)
{
    if(factor_classes)
    {
        cout << endl << "factor class " << Bi.getLowestBitsAsInteger(dim) - 1 << endl;
        for(uint16_t j = 0; j < (1 << dim); j++)
            cout << factor_classes[Bi.getLowestBitsAsInteger(dim) - 1][j] << "\t" << factor_classes_idx[Bi.getLowestBitsAsInteger(dim) - 1][j] << endl;
    }
}

SBitsArray* FactorClasses::get_factor_classes_for_Bi(SBitsArray& Bi)
{
    return factor_classes[Bi.getLowestBitsAsInteger(dim) - 1];
}

uint16_t* FactorClasses::get_factor_classes_idx_for_Bi(SBitsArray& Bi)
{
    return factor_classes_idx[Bi.getLowestBitsAsInteger_fast(dim) - 1];
}
