#ifndef HADAMAR_MATRIX_FAST_H
#define HADAMAR_MATRIX_FAST_H

#include "matrix_tensor_product.h"

class hadamar_matrix_fast : public MatrixInt16
{
public:
    hadamar_matrix_fast(uint16_t m, MatrixInt16*** &MMatrices);
    ~hadamar_matrix_fast();
    void print_M_matrices();
    int16_t* vector_left_product(int16_t *vec) override;
private:
    MatrixInt16 **M_matrices;
    uint16_t number_of_M_matrices;
};

#endif // HADAMAR_MATRIX_FAST_H
