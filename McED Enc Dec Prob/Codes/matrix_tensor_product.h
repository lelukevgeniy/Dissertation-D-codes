#ifndef MATRIX_TENSOR_PRODUCT_H
#define MATRIX_TENSOR_PRODUCT_H

#include <cstdint>
#include <iostream>

using namespace std;

typedef struct _SparseColumn{
    uint16_t size;
    uint16_t* row_idx;
    int16_t* row_val;
} SparseColumn;

class MatrixInt16{
public:
    MatrixInt16():rows(0), cols(0), matrix(nullptr), sparse_cols(nullptr) {};
    MatrixInt16(MatrixInt16& matrix_v);
    MatrixInt16(uint16_t rows_v, uint16_t cols_v);
    void init(uint16_t rows_v, uint16_t cols_v);
    void destroy();
    ~MatrixInt16();
    void set(uint16_t row_idx, uint16_t col_idx, int16_t value);
    int16_t get(uint16_t row_idx, uint16_t col_idx);
    MatrixInt16* tensor_product(MatrixInt16& B);
    virtual int16_t* vector_left_product(int16_t *vec);
    int16_t* vector_left_product_compact(int16_t *vec);
    void make_identity_matrix();
    void make_compact();

    friend ostream& operator<<(ostream& os, const MatrixInt16& matrix_v);

protected:
    uint16_t rows;
    uint16_t cols;
    int16_t** matrix;
    SparseColumn *sparse_cols;
};

#endif // MATRIX_TENSOR_PRODUCT_H
