#ifndef RMCODE_H
#define RMCODE_H

#include "SBitsArray.h"
#include "sbitsmatrix.h"
#include "vectorspacef2.h"
#include "factorclasses.h"
#include "matrix_tensor_product.h"
#include "combinatorics.h"
#include "hadamar_matrix_fast.h"
#include "hadamar_matrix.h"

class RMCode
{
public:
    RMCode(uint8_t r_val, uint8_t m_val);
    RMCode(uint8_t r_val, uint8_t m_val, bool pc);
    RMCode(uint8_t r_val, uint8_t m_val, RMCode& rm_for_init);
    ~RMCode();
    void print_monoms();
    void print_selected_monoms(char* shift, uint64_t* selectes_idx, uint64_t number_of_idx);
    void print_gen_matrix();
    SBitsMatrix* get_gen_matrix();

    SBitsArray encode(SBitsArray &input_vector);
    SBitsArray encode_fast(SBitsArray &input_vector);
    uint64_t dim() {return dimension; };
    SBitsArray linear_automotphism(SBitsArray &vector, SBitsMatrix& la);
    SBitsArray FHT_decoder(SBitsArray &noisy_vector); // only for the first order
    SBitsArray FHT_decoder_fast(SBitsArray &noisy_vector); // only for the first order
    SBitsArray FHT_decoder_fast2(SBitsArray &noisy_vector); // only for the first order
    SBitsArray ns_FHT_decoder_fast2(SBitsArray &noisy_vector); // only for the first order
    SBitsArray FHT_decoder_ret_codeword(SBitsArray &noisy_vector); // only for the first order
    SBitsArray FHT_decoder_ret_codeword_fast(SBitsArray &noisy_vector); // only for the first order
    SBitsArray FHT_decoder_ret_codeword_fast2(SBitsArray &noisy_vector); // only for the first order
    SBitsArray ns_FHT_decoder_ret_codeword_fast2(SBitsArray &noisy_vector); // only for the first order
    SBitsArray SSV_decoder(SBitsArray &noisy_vector);
    SBitsArray SSV_decoder_fast(SBitsArray &noisy_vector);
    SBitsArray RPA_decoder(SBitsArray &y, uint16_t N_max);
    SBitsArray RPA_decoder_fast(SBitsArray &y, uint16_t N_max);
    SBitsArray RPA_decoder_fast2(SBitsArray &y, uint16_t N_max);
    SBitsArray ns_RPA_decoder_fast2(SBitsArray &y, uint16_t N_max);
    SBitsArray ns_RPA_decoder_fast2_with_common_data(SBitsArray &y, uint16_t N_max, uint16_t val_r, uint16_t max_r);
    bool SSV_decoder_fast_in_one_position(uint64_t v);
    static void test_on_one_position();
    static void test_on_all_positions(float bsc_p, uint8_t proc);
    static void test_FHT(float bsc_p, uint8_t proc);
    static void test_RPA(uint8_t r, uint8_t m, uint32_t N, float bsc_p);
    static void test_RPA_by_blocks(uint8_t r, uint8_t m, uint16_t num_blocks, uint16_t weight,  uint32_t N);
    static void ns_test_RPA_by_blocks(uint8_t r1, uint8_t m1, uint8_t r2, uint8_t m2, uint16_t weight,  uint32_t N);
    SBitsArray r_syndrom_fast(SBitsArray &input_vector);
    void make_base_lin_system_fast(SBitsArray &r_s);
    static void init_factor_classes_for_RPA(uint16_t max_m);
    static void deinit_factor_classes_for_RPA();

    static void init_M_matrices_for_Hadamar(uint16_t max_m);
    static void deinit_M_matrices_for_Hadamar();

    void init_codes_for_RPA_decoder(uint16_t max_r, uint16_t max_m);
    RMCode* get_code_for_RPA_decoder(uint16_t val_r, uint16_t max_r);

private:
    SBitsArray r_syndrom(SBitsArray &input_vector);
    SBitsArray* get_monom(uint64_t i) const {return monoms[i];};
    uint64_t* get_array_monoms_after_shift(SBitsArray &shift_monom);
    uint64_t* get_array_monoms_after_shift_fast(SBitsArray &shift_monom);
    bool eval_monom_on_vector(uint64_t monom_idx, SBitsArray &vector);
    bool eval_monom_on_vector_fast(uint64_t monom_idx, SBitsArray &vector);
    void generate_all_monoms();
    void generate_gen_matrix();
    void generate_par_matrix();
    uint64_t E_r() {return ((m-r-2)/2);};
    uint64_t E_dim() {return comb_rm_dimension(E_r(),m);}
    void make_base_lin_system(SBitsArray &r_s);

public:
    static FactorClasses *fc;
    static MatrixInt16*** MMatrices;
    static uint16_t m_limit;

    FactorClasses *ns_fc;
    MatrixInt16*** ns_MMatrices;
    uint16_t ns_m_limit;

private:
    VectorSpaceF2 vs;
    uint8_t r;
    uint8_t m;
    uint64_t dimension;
    uint64_t codimension;
    SBitsArray **monoms;
    SBitsArray **gen_matrix; // columns of matrix
    SBitsArray **par_matrix; // columns of matrix
    SBitsMatrix base_lin_system;
    SBitsArray b;
    bool isCommonData;
    RMCode** codes_for_RPA_decoding;
};

#endif // RMCODE_H
