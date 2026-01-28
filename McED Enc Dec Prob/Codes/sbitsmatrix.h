#ifndef SBITSMATRIX_H
#define SBITSMATRIX_H

#include "SBitsArray.h"
#include "permutation.h"

#define SYMBOL_DELIMS ","

class SBitsMatrix
{
public:
    SBitsMatrix(quint64 number_of_rows_v, quint64 number_of_cols_v);
    SBitsMatrix(quint64 number_of_rows_v, quint64 number_of_cols_v, quint64 potential_number_of_cols_v);
    SBitsMatrix(SBitsMatrix& copy);
    void fast_copy(SBitsMatrix& copy);
    ~SBitsMatrix()
    {
        if(rows)
        {
            delete [] rows;
            rows = nullptr;
        }
    }

    SBitsArray& operator[](quint64 row_index);

    void print_nonzero_rows();
    friend ostream& operator<<(ostream& os, const SBitsMatrix& bit_array);
    friend SBitsArray operator*(SBitsMatrix& matrix, const SBitsArray& array);

    // transpose
    SBitsMatrix& T();

    quint64 get_number_of_rows()
    {
        return number_of_rows;
    }

    quint64 get_number_of_cols()
    {
        return number_of_cols;
    }

    void swap_rows(quint64 row_index1, quint64 row_index2);
    void swap_rows_fast(quint64 row_index1, quint64 row_index2);
    int64_t get_row_number_of_nonzero_in_col(quint64 col_index, quint64 start_row_index);
    int64_t get_row_number_of_nonzero_in_col_fast(quint64 col_index, quint64 start_row_index);
    quint64 low_triangular();
    quint64 low_triangular_fast();
    quint64 low_triangular_fast2();
    quint64 pseudo_systematic_form();
    void append_column(SBitsArray &col);
    void append_column_fast(SBitsArray &col);
    void random(quint64 seed);
    SBitsArray random_with_bound_column_support(quint64 seed, quint64 bound);
    void concat_columns(SBitsMatrix& matrix);
    void concat_rows(SBitsMatrix& matrix);
    SBitsMatrix* get_last_n_columns(quint64 n);
    SBitsMatrix* get_projection(SBitsArray& projection);
    SBitsMatrix* get_inverse();
    SBitsMatrix* permute_columns(Permutation &perm);
    void delete_last_rows(quint64 n_to_delete);
    void transpose();

    SBitsArray mul_left(SBitsArray& vector);
    SBitsMatrix* mul_left_matrix(SBitsMatrix& matrix);

    static SBitsMatrix* load_from_file(const char* file_name);
    static SBitsMatrix* identity_matrix(quint64 k);
    static SBitsMatrix* tensor_product(SBitsMatrix& matrix1, SBitsMatrix& matrix2);

private:
    SBitsArray *rows;
    quint64 number_of_rows;
    quint64 number_of_cols;
};

#endif // SBITSMATRIX_H
