#include "hadamar_matrix_fast.h"
#include <cstring>

hadamar_matrix_fast::hadamar_matrix_fast(uint16_t m, MatrixInt16*** &MMatrices)
    :MatrixInt16(1<<m, 1<<m),
    number_of_M_matrices(m)
{
    M_matrices = MMatrices[m];
}

hadamar_matrix_fast::~hadamar_matrix_fast()
{
}

void hadamar_matrix_fast::print_M_matrices()
{
    for(uint16_t j = 0; j < number_of_M_matrices; j++)
    {
        cout << (*M_matrices[j]) << endl;
    }
}

int16_t* hadamar_matrix_fast::vector_left_product(int16_t *vec)
{
    int16_t *result = new int16_t[cols], *temp;
    memcpy(result, vec, sizeof(int16_t)*cols);
    for(uint16_t j = 0; j < number_of_M_matrices; j++)
    {
        //temp = M_matrices[j]->vector_left_product(result);
        temp = M_matrices[j]->vector_left_product_compact(result);
        delete [] result;
        result = temp;
    }
    return result;
}
