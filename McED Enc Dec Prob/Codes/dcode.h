#ifndef DCODE_H
#define DCODE_H

#include "rmcode.h"
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <iterator>

class DCode
{
public:
    DCode():gen_matrix(nullptr){};
    ~DCode()
    {
        if(gen_matrix)
            delete gen_matrix;
        gen_matrix = nullptr;
    };

    void add_tensor_product(quint16 r1, quint16 m1, quint16 r2, quint16 m2);
    SBitsMatrix* get_K_random_subblocks(quint64 K);
    SBitsMatrix* get_subblocks(std::vector<int> &myvector);
    void print_params()
    {
        if(gen_matrix)
            cout << "D-code: dim=" <<  gen_matrix->get_number_of_rows() << ", len=" << gen_matrix->get_number_of_cols() << endl;
    };

    void print_gen_matrix()
    {
        cout << (*gen_matrix) << endl;
    };

    quint64 get_length()
    {
        if(gen_matrix)
            return gen_matrix->get_number_of_cols();
        else
            return 0;
    };

    quint64 get_dimension()
    {
        if(gen_matrix)
            return gen_matrix->get_number_of_rows();
        else
            return 0;
    };

    /* Одна итерация алгоритма ISDDecoder */
    SBitsArray isd_one_iteration(quint64 K, SBitsArray &z, bool *full_rank);

    /*  Варианты MatrixDecoder */
    static void test_MatrixDecoder(char* f_params, uint32_t N);
    static void test_MatrixDecoder_with_common_information(char* f_params, uint32_t N);

    /*  Кодирование */
    SBitsArray encode(SBitsArray &m)
    {
        return gen_matrix->mul_left(m);
    };

    /* Тестирование шифрования и расшифрования для McE(D), на основе ISDDecoder */
    static void test_EncryptionDecryption_ISD(const char* f_params, uint16_t ISD_iter, uint32_t N);
    /* Тестирование шифрования и расшифрования для McE(D), на основе MatrixDecoder */
    static void test_EncryptionDecryption_MatrixDecoder(const char* f_params, uint32_t N);


public:
    SBitsMatrix *gen_matrix;
    quint64 rm_m1;
    quint64 rm_m2;

};

#endif // DCODE_H
