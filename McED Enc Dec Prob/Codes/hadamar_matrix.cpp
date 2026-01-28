#include "hadamar_matrix.h"
#include <cstring>

hadamar_matrix::hadamar_matrix(uint16_t m)
    :MatrixInt16(1<<m, 1<<m),
    number_of_M_matrices(m)
{
    MatrixInt16 H1(2,2);
    H1.set(0,0,1);
    H1.set(0,1,1);
    H1.set(1,0,1);
    H1.set(1,1,-1);
    M_matrices = new MatrixInt16*[m];
    memset(M_matrices, 0, sizeof(MatrixInt16*)*m);

    for(uint16_t j = 0; j < m; j++)
    {
        MatrixInt16 Il(1<<(m-(j+1)), 1<<(m-(j+1))), Ir(1<<j, 1<<j), *temp;
        Il.make_identity_matrix();
        Ir.make_identity_matrix();
        temp = Il.tensor_product(H1);
        M_matrices[j] = temp->tensor_product(Ir);
        delete temp;
    }
}

hadamar_matrix::~hadamar_matrix()
{
    for(uint16_t j = 0; j < number_of_M_matrices; j++)
    {
        delete M_matrices[j];
    }
    delete [] M_matrices;
}

void hadamar_matrix::print_M_matrices()
{
    for(uint16_t j = 0; j < number_of_M_matrices; j++)
    {
        cout << (*M_matrices[j]) << endl;
    }
}

int16_t* hadamar_matrix::vector_left_product(int16_t *vec)
{
    int16_t *result = new int16_t[cols], *temp;
    memcpy(result, vec, sizeof(int16_t)*cols);
    for(uint16_t j = 0; j < number_of_M_matrices; j++)
    {
        temp = M_matrices[j]->vector_left_product(result);
        delete [] result;
        result = temp;
    }
    return result;
}
