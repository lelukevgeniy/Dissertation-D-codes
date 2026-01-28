#include "matrix_tensor_product.h""
#include <cstring>

void MatrixInt16::init(uint16_t rows_v, uint16_t cols_v)
{
    destroy();
    rows = rows_v;
    cols = cols_v;

    matrix = (int16_t**) new int16_t*[rows];
    for(int16_t i = 0; i < rows; i++)
    {
        matrix[i] = (int16_t*) new int16_t[cols];
    };

    sparse_cols = new SparseColumn[cols];
    for(int16_t i = 0; i < cols; i++)
    {
        sparse_cols[i].row_idx = nullptr;
        sparse_cols[i].row_val = nullptr;
        sparse_cols[i].size = 0;
    };
}

void MatrixInt16::destroy()
{
    if(matrix)
    {
        for(int16_t i = 0; i < rows; i++)
            if(matrix[i]) delete [] matrix[i];

        delete [] matrix;
        matrix = nullptr;
    }

    if(sparse_cols)
    {
        for(int16_t i = 0; i < cols; i++)
        {
            if(sparse_cols[i].row_idx) delete [] sparse_cols[i].row_idx;
            if(sparse_cols[i].row_val) delete [] sparse_cols[i].row_val;
        }
        delete [] sparse_cols;
        sparse_cols = nullptr;
    }
}

MatrixInt16::MatrixInt16(uint16_t rows_v, uint16_t cols_v)
    :MatrixInt16()
{
    init(rows_v, cols_v);
}

MatrixInt16::MatrixInt16(MatrixInt16& matrix_v)
{
    MatrixInt16(matrix_v.rows, matrix_v.cols);
    for(int16_t i = 0; i < rows; i++)
        memcpy(matrix[i], matrix_v.matrix[i], sizeof(int16_t)*cols);
}

MatrixInt16::~MatrixInt16()
{
    destroy();
};

void MatrixInt16::set(uint16_t row_idx, uint16_t col_idx, int16_t value)
{
    matrix[row_idx][col_idx]= value;
}

int16_t MatrixInt16::get(uint16_t row_idx, uint16_t col_idx)
{
    return matrix[row_idx][col_idx];
}

ostream& operator<<(ostream& os, const MatrixInt16& matrix_v)
{
    if(matrix_v.matrix)
    {
        for(int16_t i = 0; i < matrix_v.rows; i++)
        {
            for(int16_t j = 0; j < matrix_v.cols; j++)
                os << (int16_t) matrix_v.matrix[i][j] << " ";
            cout << endl;
        }
    }else{
        os << "null";
    }
    return os;
}

MatrixInt16* MatrixInt16::tensor_product(MatrixInt16& B)
{
    MatrixInt16* tens_prod = new MatrixInt16(rows*B.rows, cols*B.cols);

    for(int16_t i = 0; i < rows; i++)
    {
        for(int16_t j = 0; j < cols; j++)
        {
            for(int16_t Bi = 0; Bi < B.rows; Bi++)
            {
                for(int16_t Bj = 0; Bj < B.cols; Bj++)
                {
                    tens_prod->set(i*B.rows + Bi, j*B.cols + Bj, matrix[i][j]*B.matrix[Bi][Bj]);
                }
            }
        }
    }

    return tens_prod;
}

int16_t* MatrixInt16::vector_left_product(int16_t *vec)
{
    int16_t* result = new int16_t[cols];
    memset(result, 0, sizeof(int16_t)*cols);

    for(int16_t j = 0; j < cols; j++)
    {
        for(int16_t i = 0; i < rows; i++)
        {
            result[j] += (vec[i] * matrix[i][j]);
        };
    };

    return result;
};

int16_t* MatrixInt16::vector_left_product_compact(int16_t *vec)
{
    int16_t* result = new int16_t[cols];

    for(int16_t j = 0; j < cols; j++)
    {
        result[j] = vec[sparse_cols[j].row_idx[0]] * matrix[sparse_cols[j].row_idx[0]][j] + vec[sparse_cols[j].row_idx[1]] * matrix[sparse_cols[j].row_idx[1]][j];
    };

    return result;
}

void MatrixInt16::make_identity_matrix()
{
    if((cols != rows) || (cols == 0) || (rows == 0))
    {
        cout << "Matrix is not square or number of rows of cols is equal to zero " << endl;
        return;
    };

    for(int16_t i = 0; i < rows; i++)
    {
        memset(matrix[i], 0, sizeof(int16_t)*cols);
        matrix[i][i] = 1;
    }
}

void MatrixInt16::make_compact()
{
    uint16_t non_zero_count = 0;

    for(int16_t i = 0; i < cols; i++)
    {
        non_zero_count = 0;
        for(int16_t j = 0; j < rows; j++)
        {
            if(matrix[j][i] != 0) non_zero_count++;
        };

        if(non_zero_count)
        {
            sparse_cols[i].row_idx = new uint16_t[non_zero_count];
            sparse_cols[i].row_val = new int16_t[non_zero_count];
            sparse_cols[i].size = non_zero_count;

            non_zero_count = 0;
            for(int16_t j = 0; j < rows; j++)
            {
                if(matrix[j][i] != 0)
                {
                    sparse_cols[i].row_idx[non_zero_count] = j;
                    sparse_cols[i].row_val[non_zero_count] = matrix[j][i];
                    non_zero_count++;
                };
            };
        }
    }
}
