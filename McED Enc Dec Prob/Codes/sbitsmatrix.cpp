#include "sbitsmatrix.h"
#include <sstream>
#include <fstream>
#include <iostream>
using namespace std;

SBitsMatrix::SBitsMatrix(quint64 number_of_rows_v, quint64 number_of_cols_v)
    :number_of_rows(0),
          number_of_cols(0),
          rows(nullptr)
{
    if((number_of_cols_v > 0) && (number_of_rows_v > 0))
    {
        number_of_cols = number_of_cols_v;
        number_of_rows = number_of_rows_v;
        rows = new SBitsArray[number_of_rows_v];
        for(int i = 0; i < number_of_rows; i++)
            rows[i].allocateBufferForBits(number_of_cols);
    }
}

SBitsMatrix::SBitsMatrix(quint64 number_of_rows_v, quint64 number_of_cols_v, quint64 potential_number_of_cols_v)
    :number_of_rows(0),
          number_of_cols(0),
          rows(nullptr)
{
    if((number_of_cols_v > 0) && (number_of_rows_v > 0))
    {
        number_of_cols = number_of_cols_v;
        number_of_rows = number_of_rows_v;
        rows = new SBitsArray[number_of_rows_v];
        for(int i = 0; i < number_of_rows; i++)
        {
            rows[i].allocateBufferForPotentialBits(number_of_cols_v, potential_number_of_cols_v);
        }
    }
}

SBitsMatrix::SBitsMatrix(SBitsMatrix& copy)
    :SBitsMatrix(copy.number_of_rows, copy.number_of_cols)
{
    for(uint64_t i = 0; i < copy.number_of_rows; i++)
    {
        rows[i] = copy[i];
    };
}

void SBitsMatrix::fast_copy(SBitsMatrix& copy)
{
    for(uint64_t i = 0; i < copy.number_of_rows; i++)
    {
        rows[i].fast_copy_without_check(copy[i]);
    };
}

SBitsArray& SBitsMatrix::operator[](quint64 row_index)
{
    return rows[row_index];
}

ostream& operator<<(ostream& os, const SBitsMatrix& bit_array)
{
    for(int i = 0; i < bit_array.number_of_rows; i++)
        os << bit_array.rows[i] << endl;

    return os;
}

void SBitsMatrix::print_nonzero_rows()
{
    for(int i = 0; i < number_of_rows; i++)
    {
        if(rows[i].getHammingWeight() != 0)
            cout << rows[i] << endl;
    }
}

 SBitsMatrix& SBitsMatrix::T()
{
    SBitsArray *cols;
    if(rows != nullptr)
    {
        cols = new SBitsArray[number_of_cols];

        for(int i = 0; i < number_of_cols; i++)
        {
            cols[i].allocateBufferForBits(number_of_rows);

            for(int j = 0; j < number_of_rows; j++)
                cols[i].setBitValue(j, rows[j].getBit(i));
        }

        delete [] rows;
        rows = cols;

        quint64 temp = number_of_cols;
        number_of_cols = number_of_rows;
        number_of_rows = temp;
    }

    return *this;
}

SBitsArray operator*(SBitsMatrix& matrix, const SBitsArray& array)
{
	SBitsArray result(0);

    if(matrix.number_of_cols != array.getSizeInBits())
    {
        cout << "matrix * vector: wrong dimensions" << endl;
        return result;
    }

    result.allocateBufferForBits(matrix.number_of_rows);

    for(int i = 0; i < matrix.number_of_rows; i++)
    {
        SBitsArray temp(matrix.number_of_cols);
        temp = matrix[i];
        temp.andWithSBitArray(array);
        result.setBitValue(i, temp.getHammingWeightMod2());
    }

    return result;
}

void SBitsMatrix::swap_rows(quint64 row_index1, quint64 row_index2)
{
    bool temp_val;
    for(uint64_t i = 0; i < number_of_cols; i++)
    {
        temp_val = (*this)[row_index1].getBit(i);
        (*this)[row_index1].setBitValue(i, (*this)[row_index2].getBit(i));
        (*this)[row_index2].setBitValue(i, temp_val);
    }
}

void SBitsMatrix::swap_rows_fast(quint64 row_index1, quint64 row_index2)
{
    SBitsArray temp_val(rows[row_index1]);
    (*this)[row_index1].fast_copy(rows[row_index2]);
    (*this)[row_index2].fast_copy(temp_val);
}

int64_t SBitsMatrix::get_row_number_of_nonzero_in_col(quint64 col_index, quint64 start_row_index)
{
    int64_t row_number = -1;

    if(col_index < number_of_cols)
    {
        for(quint64 i = start_row_index; i < number_of_rows; i++)
        {
            if((*this)[i].getBit(col_index))
            {
                row_number = i;
                break;
            }
        }
    }
    return row_number;
}

int64_t SBitsMatrix::get_row_number_of_nonzero_in_col_fast(quint64 col_index, quint64 start_row_index)
{
    int64_t row_number = -1;

    if(col_index < number_of_cols)
    {
        for(quint64 i = start_row_index; i < number_of_rows; i++)
        {
            if((*this)[i].getBit_fast(col_index))
            {
                row_number = i;
                break;
            }
        }
    }
    return row_number;
}

quint64 SBitsMatrix::low_triangular()
{
    quint64 current_row = 0;

    for(quint64 i = 0; i < number_of_rows; i++)
    {
        int64_t first_nonzero = get_row_number_of_nonzero_in_col(i, current_row);

        if(first_nonzero != -1)
        {
            if(first_nonzero != current_row)
            {
                swap_rows(current_row, first_nonzero);
            };

            for(quint64 j = current_row + 1; j < number_of_rows; j++)
            {
                if((*this)[j].getBit(i))
                {
                    (*this)[j].xorWithSBitArray((*this)[current_row]);
                }
            }

            current_row++;
        }else{
            continue;
        }
    }

    return current_row;
}

quint64 SBitsMatrix::low_triangular_fast()
{
    quint64 current_row = 0;

    for(quint64 i = 0; i < number_of_rows; i++)
    {
        int64_t first_nonzero = get_row_number_of_nonzero_in_col_fast(i, current_row);

        if(first_nonzero != -1)
        {
            if(first_nonzero != current_row)
            {
                swap_rows_fast(current_row, first_nonzero);
            };

            for(quint64 j = current_row + 1; j < number_of_rows; j++)
            {
                if((*this)[j].getBit_fast(i))
                {
                    (*this)[j].xorWithSBitArray_fast((*this)[current_row]);
                }
            }
            current_row++;
        }
    }

    return current_row;
}

quint64 SBitsMatrix::low_triangular_fast2()
{
    quint64 current_row = 0;

    for(quint64 i = 0; i < number_of_cols; i++)
    {
        int64_t first_nonzero = get_row_number_of_nonzero_in_col_fast(i, current_row);

        if(first_nonzero != -1)
        {
            if(first_nonzero != current_row)
            {
                swap_rows_fast(current_row, first_nonzero);
            };

            for(quint64 j = current_row + 1; j < number_of_rows; j++)
            {
                if((*this)[j].getBit_fast(i))
                {
                    (*this)[j].xorWithSBitArray_fast((*this)[current_row]);
                }
            }
            current_row++;
        }
    }

    return current_row;
}

quint64 SBitsMatrix::pseudo_systematic_form()
{
    quint64 current_row = low_triangular_fast2();
    if(current_row < number_of_rows)
        return current_row;

    //cout << *this << endl;

    for(quint64 i = number_of_rows - 1; i >= 1; i--)
    {
        quint64 idx_fnz = rows[i].getIndexOfFirstNonzero();
        for(qint64 j = i - 1; j >= 0; j--)
        {
            if(rows[j].getBit_fast(idx_fnz) == true)
            {
                rows[j].xorWithSBitArray_fast(rows[i]);
            }
        }
    }

    return current_row;
}

void SBitsMatrix::append_column(SBitsArray &col)
{
    for(uint64_t i = 0; i < col.getSizeInBits(); i++)
    {
        rows[i].concat_with_bit(col.getBit(i));
    }
    number_of_cols++;
}

void SBitsMatrix::append_column_fast(SBitsArray &col)
{
    for(uint64_t i = 0; i < col.getSizeInBits(); i++)
    {
        rows[i].concat_with_bit_fast(col.getBit_fast(i));
    }
    number_of_cols++;
}

void SBitsMatrix::random(quint64 seed)
{
    for(quint32 i = 0; i < number_of_rows; i++)
        rows[i].randomBitArray(seed + i);
};

SBitsArray SBitsMatrix::random_with_bound_column_support(quint64 seed, quint64 bound)
{
    SBitsArray support(number_of_cols);
    support.toZeroAllBits();
    srand(seed);
    transpose();
    for(quint32 i = 0; i < bound; i++)
        rows[rand()%number_of_rows].randomBitArray(seed + i);

    for(quint32 i = 0; i < number_of_rows; i++)
    {
        if(rows[i].getHammingWeight_fast() != 0)
            support.setBit_fast(i);
    };
    transpose();
    return support;
}

void SBitsMatrix::concat_columns(SBitsMatrix& matrix)
{
    if(number_of_rows != matrix.get_number_of_rows())
        return;

    for(quint64 i = 0; i < number_of_rows; i++)
    {
        rows[i].concat(matrix[i]);
    };
    number_of_cols = number_of_cols + matrix.get_number_of_cols();
}

void SBitsMatrix::concat_rows(SBitsMatrix& matrix)
{
    if(number_of_cols != matrix.get_number_of_cols())
        return;

    SBitsArray *new_rows = new SBitsArray[number_of_rows + matrix.get_number_of_rows()];

    for(quint64 i = 0; i < number_of_rows; i++)
    {
        new_rows[i].allocateBufferForBits(number_of_cols);
        new_rows[i].fast_copy_without_check(rows[i]);
    };

    for(quint64 i = 0; i <matrix.get_number_of_rows(); i++)
    {
        new_rows[number_of_rows + i].allocateBufferForBits(number_of_cols);
        new_rows[number_of_rows + i].fast_copy_without_check(matrix.rows[i]);
    };

    delete [] rows;
    rows = new_rows;
    number_of_rows += matrix.get_number_of_rows();
}

SBitsMatrix* SBitsMatrix::load_from_file(const char* file_name)
{
    SBitsMatrix *matrix = nullptr;

    {
        std::ifstream file;
        file.open(file_name);
        std::string delims = SYMBOL_DELIMS;
        auto number_of_rows = 0;
        auto number_of_cols = 0;

        // get number of rows and columns
        for(std::string line; std::getline(file, line);)
        {
            number_of_rows++;
            auto start = 0U;
            auto end = line.find(delims);
            auto number_of_columns_temp = 0;
            while (end != std::string::npos)
            {
                number_of_columns_temp++;
                start = end + delims.length();
                end = line.find(delims, start);
            }
            number_of_cols = max(number_of_cols, number_of_columns_temp + 1);
        }
        file.close();

        matrix = new SBitsMatrix(number_of_rows, number_of_cols);
        number_of_rows = 0;
        number_of_cols = 0;
        // read matrix
        file.open(file_name);
        for(std::string line; std::getline(file, line);)
        {
            auto start = 0U;
            auto end = line.find(delims);
            while (end != std::string::npos)
            {
                (*matrix)[number_of_rows].setBitValue(number_of_cols, atoi(line.substr(start, end - start).c_str()));
                start = end + delims.length();
                end = line.find(delims, start);
                number_of_cols++;
            }
            (*matrix)[number_of_rows].setBitValue(number_of_cols, atoi(line.substr(start, line.length() - start).c_str()));
            number_of_cols = 0;
            number_of_rows++;
        }
    }

    return matrix;
}

SBitsMatrix* SBitsMatrix::identity_matrix(quint64 k)
{
    SBitsMatrix* matrix = new SBitsMatrix(k, k);
    for(quint64 i = 0; i < k; i++)
        (*matrix)[i].setBit_fast(i);

    return matrix;
}

SBitsMatrix* SBitsMatrix::tensor_product(SBitsMatrix& matrix1, SBitsMatrix& matrix2)
{
    quint64 number_of_rows_in_matrix2 = matrix2.get_number_of_rows();
    quint64 number_of_cols_in_matrix2 = matrix2.get_number_of_cols();
    SBitsMatrix* matrix = new SBitsMatrix(matrix1.get_number_of_rows()*number_of_rows_in_matrix2,
                                          matrix1.get_number_of_cols()*number_of_cols_in_matrix2);

    for(quint64 i = 0; i < matrix1.get_number_of_rows(); i++)
        for(quint64 j = 0; j < matrix1.get_number_of_cols(); j++)
        {
            if(matrix1[i].getBit_fast(j))
            {
                for(quint64 k = 0; k < number_of_rows_in_matrix2; k++)
                {
                    (*matrix).rows[i*number_of_rows_in_matrix2 + k].setSubvectorInPosition_fast(j*number_of_cols_in_matrix2, matrix2.rows[k]);
                }
            }
        }

    return matrix;
}

void SBitsMatrix::delete_last_rows(quint64 n_to_delete)
{
    if(n_to_delete >= number_of_rows)
        return;

    SBitsArray *new_rows = new SBitsArray[number_of_rows - n_to_delete];

    for(quint64 i = 0; i < number_of_rows - n_to_delete; i++)
    {
        new_rows[i].allocateBufferForBits(number_of_cols);
        new_rows[i].fast_copy_without_check(rows[i]);
    };

    delete [] rows;
    rows = new_rows;
    number_of_rows -= n_to_delete;
}

SBitsMatrix* SBitsMatrix::permute_columns(Permutation& perm)
{
    if (perm.get_length() != get_number_of_cols())
        return NULL;

    SBitsMatrix* matrix = new SBitsMatrix(*this);

    for (uint16_t i = 0; i < get_number_of_rows(); i++)
        (*matrix).rows[i].permute(perm);

    return matrix;
}

void SBitsMatrix::transpose()
{
    SBitsArray *new_rows = new SBitsArray[number_of_cols];

    for(quint64 i = 0; i < number_of_cols; i++)
    {
        new_rows[i].allocateBufferForBits(number_of_rows);
    };

    for(quint64 i = 0; i < number_of_rows; i++)
    {
        for(quint64 j = 0; j < number_of_cols; j++)
            new_rows[j].setBitValue_fast(i, rows[i].getBit_fast(j));
    };

    delete [] rows;
    rows = new_rows;
    qint64 temp = number_of_rows;
    number_of_rows = number_of_cols;
    number_of_cols = temp;
}

SBitsArray SBitsMatrix::mul_left(SBitsArray& vector)
{
    SBitsArray result(number_of_cols);

    for(quint64 i = 0; i < number_of_cols; i++)
    {
        quint64 temp = 0;
        for(quint64 j = 0; j < number_of_rows; j++)
        {
            temp = temp ^ ((vector.getBit_fast(j) && rows[j].getBit_fast(i)) == true ? 1 : 0);
            //temp = temp + ((vector.getBit_fast(j) && rows[j].getBit_fast(i)) == true ? 1 : 0);
        };

        //cout << temp << "#";
        result.setBitValue_fast(i, temp == 0 ? false : true);
        //result.setBitValue_fast(i, (temp%2) == 0 ? false : true);
    };

    return result;
}

SBitsMatrix* SBitsMatrix::mul_left_matrix(SBitsMatrix& matrix)
{
    SBitsMatrix *result = new SBitsMatrix(number_of_rows, number_of_cols);

    for(uint64_t i = 0; i < number_of_rows; i++)
        result->rows[i] = mul_left(matrix[i]);

    return result;
}

SBitsMatrix* SBitsMatrix::get_last_n_columns(quint64 n)
{
    SBitsMatrix* matrix = new SBitsMatrix(number_of_rows, n);

    for(quint64 i = 0; i < number_of_rows; i++)
    {
        for(quint64 j = 0; j < n; j++)
        {
            (*matrix)[i].setBitValue_fast(j, rows[i].getBit_fast(number_of_cols - n + j));
        };
    };

    return matrix;
}

SBitsMatrix* SBitsMatrix::get_projection(SBitsArray& projection)
{
    quint64 cols_in_projection = projection.getHammingWeight_fast();
    SBitsMatrix* matrix = new SBitsMatrix(number_of_rows, cols_in_projection);
    quint64 columns_counter = 0;
    for(quint64 i = 0; i < number_of_rows; i++)
    {
        columns_counter = 0;
        for(quint64 j = 0; j < number_of_cols; j++)
        {
            if(projection.getBit_fast(j))
            {
                (*matrix)[i].setBitValue_fast(columns_counter, rows[i].getBit_fast(j));
                columns_counter++;
            }
        };
    };

    //cout << *matrix << endl;

    return matrix;
}

SBitsMatrix* SBitsMatrix::get_inverse()
{
    if(number_of_cols != number_of_rows)
        return nullptr;

    SBitsMatrix copy(number_of_cols, number_of_cols);
    SBitsMatrix* id = SBitsMatrix::identity_matrix(number_of_cols);
    copy.fast_copy(*this);
    copy.concat_columns(*id);
    delete id;
    copy.pseudo_systematic_form();

    if(copy[number_of_rows - 1].getIndexOfFirstNonzero() == (number_of_rows - 1))
        return copy.get_last_n_columns(number_of_rows);
    else
        return nullptr;
}
