
#include "rmcode.h"


FactorClasses * RMCode::fc = nullptr;
MatrixInt16 *** RMCode::MMatrices = nullptr;
uint16_t RMCode::m_limit = 0;

/*
SBitsArray parallel_proc(RMCode rm, SBitsArray cw, uint16_t iter)
{
    return rm.RPA_decoder(cw, iter);
}
*/

RMCode::RMCode(uint8_t r_val, uint8_t m_val, bool pc)
    :r(r_val),
      m(m_val),
      monoms(nullptr),
      gen_matrix(nullptr),
      par_matrix(nullptr),
      vs(m_val),
      base_lin_system(1 + E_dim() + m_val*E_dim(),E_dim()),
      b(1 + E_dim() + m_val*E_dim()),
      isCommonData(false)
 {
    if(r_val < m_val)
    {
        r = r_val;
        m = m_val;
        dimension = comb_rm_dimension(r,m);
        codimension = (1<<m) - dimension;
        generate_all_monoms();
        generate_gen_matrix();
        if(pc)
        {
            generate_par_matrix();
        }
    }

    /*
    cout << "dimension = " << dimension << endl;
    cout << "codimension = " << codimension << endl;
    cout << "E_dim = " << E_dim() << endl;
    /**/
}

RMCode::RMCode(uint8_t r_val, uint8_t m_val, RMCode& rm_for_init)
    :r(r_val),
      m(m_val),
      monoms(nullptr),
      gen_matrix(nullptr),
      par_matrix(nullptr),
      vs(rm_for_init.vs), //
      base_lin_system(1 + E_dim() + m_val*E_dim(),E_dim()),
      b(1 + E_dim() + m_val*E_dim()),
      isCommonData(true)
{
    if(r_val < m_val)
    {
        r = r_val;
        m = m_val;
        dimension = comb_rm_dimension(r,m);
        codimension = (1<<m) - dimension;
        monoms = rm_for_init.monoms;
        gen_matrix = rm_for_init.gen_matrix;
        par_matrix = rm_for_init.par_matrix;
        //generate_all_monoms();
        //generate_gen_matrix();
        //generate_par_matrix();
    }

    /*
    cout << "dimension = " << dimension << endl;
    cout << "codimension = " << codimension << endl;
    cout << "E_dim = " << E_dim() << endl;
    /**/
}

RMCode::RMCode(uint8_t r_val, uint8_t m_val)
    :r(r_val),
      m(m_val),
      monoms(nullptr),
      gen_matrix(nullptr),
      par_matrix(nullptr),
      vs(m_val),
      base_lin_system(1 + E_dim() + m_val*E_dim(),E_dim()),
      b(1 + E_dim() + m_val*E_dim()),
      isCommonData(false)
{
    if(r_val < m_val)
    {
        r = r_val;
        m = m_val;
        dimension = comb_rm_dimension(r,m);
        codimension = (1<<m) - dimension;
        generate_all_monoms();
        generate_gen_matrix();
        generate_par_matrix();
    }

    /*
    cout << "dimension = " << dimension << endl;
    cout << "codimension = " << codimension << endl;
    cout << "E_dim = " << E_dim() << endl;
    /**/
}

RMCode::~RMCode()
{
    if(isCommonData)
    {
        if(monoms)
        {
            delete [] monoms;
        }

        if(gen_matrix)
        {
            for(uint64_t arg_counter = 0; arg_counter < (1<<m); arg_counter++) // цикл по аргументам мономов
            {
                if(gen_matrix[arg_counter])
                    delete gen_matrix[arg_counter];
                if(par_matrix)
                    if(par_matrix[arg_counter])
                        delete par_matrix[arg_counter];
            };

            delete [] gen_matrix;
            delete [] par_matrix;
        };
    }
}

SBitsArray RMCode::encode(SBitsArray &input_vector)
{
    SBitsArray result(1<<m);

    for(uint64_t i = 0; i < (1<<m); i++)
    {
        result.setBitValue(i, input_vector.andWithSBitArray(*(gen_matrix[i]), input_vector).getHammingWeightMod2());
    }

    return result;
}

SBitsArray RMCode::encode_fast(SBitsArray &input_vector)
{
    SBitsArray result(1<<m);

    for(uint64_t i = 0; i < (1<<m); i++)
    {
        result.setBitValue_fast(i, input_vector.andWithSBitArray_fast(*(gen_matrix[i]), input_vector).getHammingWeightMod2_fast());
    }

    return result;
}


void RMCode::generate_all_monoms()
{
    uint64_t current_index = 0;

    if(dimension)
    {
        if(monoms) delete[] monoms;
        monoms = new SBitsArray*[1<<m];

        for(uint8_t w = 0; w <= m; w++) // цикл по весам
        {
            monoms[current_index++] = vs.get_first_vector_for_weight(w);

            uint64_t one_block_size = comb_choose(m, w);

            for(uint64_t i = 1; i < one_block_size; i++)
            {
                monoms[current_index++] = vs.get_next_vector_for_weight(w);
            }
        }
    }
}

void RMCode::generate_gen_matrix()
{
    gen_matrix = new SBitsArray*[1<<m];

    for(uint64_t arg_counter = 0; arg_counter < (1<<m); arg_counter++) // цикл по аргументам мономов
    {
        uint8_t v = 0;
        gen_matrix[arg_counter] = new SBitsArray(dimension);
        for(uint64_t monom_counter = 0; monom_counter < dimension; monom_counter++) // цикл по мономам
        {
            SBitsArray supp_inclusion = SBitsArray::andWithSBitArray(*(monoms[monom_counter]), vs[arg_counter]);
            gen_matrix[arg_counter]->setBitValue(monom_counter, supp_inclusion == *(monoms[monom_counter]));
        }
    }
}

void RMCode::generate_par_matrix()
{
    par_matrix = new SBitsArray*[1<<m];

    for(uint64_t arg_counter = 0; arg_counter < (1<<m); arg_counter++) // цикл по аргументам мономов
    {
        uint8_t v = 0;
        par_matrix[arg_counter] = new SBitsArray(((1<<m)-dimension));
        for(uint64_t monom_counter = 0; monom_counter < ((1<<m)-dimension); monom_counter++) // цикл по мономам
        {
            SBitsArray supp_inclusion = SBitsArray::andWithSBitArray(*(monoms[monom_counter]), vs[arg_counter]);
            par_matrix[arg_counter]->setBitValue(monom_counter, supp_inclusion == *(monoms[monom_counter]));
        }
    }
}

void RMCode::print_monoms()
{
    for(uint64_t i = 0; i < dimension; i++)
        cout << (*monoms[i]) << endl;
}

void RMCode::print_selected_monoms(char* shift, uint64_t* selectes_idx, uint64_t number_of_idx)
{
    for(uint64_t i = 0; i < codimension; i++)
    {
        cout << shift << (*monoms[i]);
        for(uint64_t j = 0; j < number_of_idx; j++)
        {
            if(i == selectes_idx[j])
            {
                cout << "*";
                break;
            }
        }
        cout << endl;
    }
}

void RMCode::print_gen_matrix()
{
    for(uint64_t monom_counter = 0; monom_counter < dimension; monom_counter++) // цикл по мономам
    {
        for(uint64_t arg_counter = 0; arg_counter < (1<<m); arg_counter++) // цикл по аргументам мономов
        {
            cout << (gen_matrix[arg_counter]->getBit(monom_counter) ? "1" : "0" );
        }
        cout << endl;
    }
}

SBitsMatrix* RMCode::get_gen_matrix()
{
    SBitsMatrix* matrix = new SBitsMatrix(dimension, 1<<m);

    for(quint64 i = 0; i < dimension; i++)
        for(quint64 j = 0; j < 1<<m; j++)
        {
            if(gen_matrix[j]->getBit_fast(i))
                (*matrix)[i].setBit_fast(j);
        };

    return matrix;
}

SBitsArray RMCode::r_syndrom(SBitsArray &input_vector)
{
    SBitsArray r_s(codimension);

    //cout << "monom\tvalue"<< endl;
    for(uint64_t monom_counter = 0; monom_counter < codimension; monom_counter++) // цикл по мономам
    {
        uint8_t v = 0;
        for(uint64_t arg_counter = 0; arg_counter < (1<<m); arg_counter++) // цикл по аргументам мономов
        {
            SBitsArray supp_inclusion = SBitsArray::andWithSBitArray(*(monoms[monom_counter]), vs[arg_counter]);
            if((supp_inclusion == *(monoms[monom_counter])) && (input_vector.getBit(arg_counter)))
            {
                //cout << *(monoms[monom_counter]) << ", " << vs[arg_counter] << endl;
                v = v ^ 1;
            }
        }
        r_s.setBitValue(monom_counter, bool(v));
        //cout << *(monoms[monom_counter]) << "\t" << uint32_t(v) << endl;
    }

    return r_s;
}

SBitsArray RMCode::r_syndrom_fast(SBitsArray &input_vector)
{
    SBitsArray r_s(codimension);

    //cout << "monom\tvalue"<< endl;
    for(uint64_t monom_counter = 0; monom_counter < codimension; monom_counter++) // цикл по мономам
    {
        uint8_t v = 0;
        for(uint64_t arg_counter = 0; arg_counter < (1<<m); arg_counter++) // цикл по аргументам мономов
        {
            SBitsArray supp_inclusion = SBitsArray::andWithSBitArray_fast(*(monoms[monom_counter]), vs[arg_counter]);
            if((supp_inclusion == *(monoms[monom_counter])) && (input_vector.getBit_fast(arg_counter)))
            {
                //cout << *(monoms[monom_counter]) << ", " << vs[arg_counter] << endl;
                v = v ^ 1;
            }
        }
        r_s.setBitValue_fast(monom_counter, bool(v));
        //cout << *(monoms[monom_counter]) << "\t" << uint32_t(v) << endl;
    }

    return r_s;
}

uint64_t* RMCode::get_array_monoms_after_shift(SBitsArray &shift_monom)
{
    uint64_t *array = new uint64_t[E_dim()];

    for(uint64_t i = 0; i <E_dim(); i++)
    {
        SBitsArray shift = SBitsArray::orWithSBitArray(*(monoms[i]), shift_monom);

        for(uint64_t mon = 0; mon < codimension; mon++)
        {
            if(SBitsArray::xorWithSBitArray(*(monoms[mon]), shift).getHammingWeight() == 0)
            {
                array[i] = mon;
                break;
            }
        }
    };

    return array;
}

uint64_t* RMCode::get_array_monoms_after_shift_fast(SBitsArray &shift_monom)
{
    uint64_t *array = new uint64_t[E_dim()];

    for(uint64_t i = 0; i <E_dim(); i++)
    {
        SBitsArray shift(shift_monom.getSizeInBits());
        shift.fast_copy(SBitsArray::orWithSBitArray_fast(*(monoms[i]), shift_monom));

        for(uint64_t mon = 0; mon < codimension; mon++)
        {
            if(SBitsArray::xorWithSBitArray_fast(*(monoms[mon]), shift).getHammingWeight_fast() == 0)
            {
                array[i] = mon;
                break;
            }
        }
    };

    return array;
}

bool RMCode::eval_monom_on_vector(uint64_t monom_idx, SBitsArray &vector)
{
    bool result = SBitsArray::andWithSBitArray(*(monoms[monom_idx]), vector) == (*(monoms[monom_idx]));
    //cout << "evaluation " << *(monoms[monom_idx]) << " on " << vector << " is " << uint32_t(result) << endl;
    return result;
}

bool RMCode::eval_monom_on_vector_fast(uint64_t monom_idx, SBitsArray &vector)
{
    bool result = SBitsArray::andWithSBitArray_fast(*(monoms[monom_idx]), vector) == (*(monoms[monom_idx]));
    //cout << "evaluation " << *(monoms[monom_idx]) << " on " << vector << " is " << uint32_t(result) << endl;
    return result;
}

void RMCode::make_base_lin_system(SBitsArray &r_s)
{
    uint64_t *arr0 = get_array_monoms_after_shift(*get_monom(0));
    base_lin_system[0] = r_s.getProjection(arr0, E_dim());
    delete [] arr0;
    base_lin_system[0].setBitValue(E_dim(), true);

    for(uint64_t M = 0; M < E_dim(); M++) // цикл по мономам  степени не выше E_r
    {
        uint64_t *arr1 = get_array_monoms_after_shift(*get_monom(M));
        base_lin_system[1 + (m+1)*M] = r_s.getProjection(arr1, E_dim());
        delete [] arr1;
        for(uint64_t l = 0; l < m; l++) // цикл по позициям
        {
            SBitsArray Mx_l(m);
            Mx_l = SBitsArray::orWithSBitArray(*get_monom(M), *get_monom(1 + l));
            uint64_t *arr2 = get_array_monoms_after_shift(Mx_l);
            base_lin_system[1 + (m+1)*M + 1 + l] = r_s.getProjection(arr2, E_dim());
            delete [] arr2;
        }
    };
}

void RMCode::make_base_lin_system_fast(SBitsArray &r_s)
{
    uint64_t *arr0 = get_array_monoms_after_shift_fast(*get_monom(0));
    base_lin_system[0].fast_copy(r_s.getProjection_fast(arr0, E_dim()));
    delete [] arr0;
    base_lin_system[0].setBitValue(E_dim(), true);

    for(uint64_t M = 0; M < E_dim(); M++) // цикл по мономам  степени не выше E_r
    {
        uint64_t *arr1 = get_array_monoms_after_shift_fast(*get_monom(M));
        base_lin_system[1 + (m+1)*M].fast_copy(r_s.getProjection_fast(arr1, E_dim()));
        delete [] arr1;
        for(uint64_t l = 0; l < m; l++) // цикл по позициям
        {
            SBitsArray Mx_l(m);
            //Mx_l.fast_copy(SBitsArray::orWithSBitArray_fast(*get_monom(M), *get_monom(1 + l)));
            Mx_l.fast_copy(*(monoms[M]));
            Mx_l.orWithSBitArray_fast(*(monoms[1+l]));
            uint64_t *arr2 = get_array_monoms_after_shift_fast(Mx_l);
            base_lin_system[1 + (m+1)*M + 1 + l].fast_copy(r_s.getProjection_fast(arr2, E_dim()));
            delete [] arr2;
        }
    };
}

SBitsArray RMCode::FHT_decoder(SBitsArray &noisy_vector)
{
    SBitsArray message(m+1);
    if(r == 1)
    {
        hadamar_matrix H(m);
        int16_t* Q = new int16_t[1<<m], *WH;
        for(uint16_t i = 0; i < (1 << m); i++)
        {
            if(noisy_vector.getBit(i))
                Q[i] = -1;
            else
                Q[i] = 1;
        };

        WH = H.vector_left_product(Q);
        delete [] Q;

        int16_t max = 0, idx_max = -1;
        for(uint16_t i = 0; i < (1 << m); i++)
        {
            if(std::abs(WH[i]) > max)
            {
                max = std::abs(WH[i]);
                idx_max = i;
            };
        };
        message.setBitsValueFromQint64(0, m+1, (WH[idx_max] >= 0 ? (idx_max << 1) : ((idx_max << 1) | 1)));
        delete [] WH;
    };
    return message;
}

SBitsArray RMCode::FHT_decoder_fast(SBitsArray &noisy_vector)
{
    SBitsArray message(m+1);
    if(r == 1)
    {
        hadamar_matrix H(m);
        int16_t* Q = new int16_t[1<<m], *WH;
        for(uint16_t i = 0; i < (1 << m); i++)
        {
            if(noisy_vector.getBit_fast(i))
                Q[i] = -1;
            else
                Q[i] = 1;
        };

        WH = H.vector_left_product(Q);
        delete [] Q;

        int16_t max = 0, idx_max = -1;
        for(uint16_t i = 0; i < (1 << m); i++)
        {
            if(std::abs(WH[i]) > max)
            {
                max = std::abs(WH[i]);
                idx_max = i;
            };
        };
        message.setBitsValueFromQint64_fast(WH[idx_max] >= 0 ? (idx_max << 1) : ((idx_max << 1) | 1));
        delete [] WH;
    };
    return message;
}

SBitsArray RMCode::FHT_decoder_fast2(SBitsArray &noisy_vector)
{
    SBitsArray message(m+1);
    if(r == 1)
    {
        hadamar_matrix_fast H(m, MMatrices);
        int16_t* Q = new int16_t[1<<m], *WH;
        for(uint16_t i = 0; i < (1 << m); i++)
        {
            if(noisy_vector.getBit_fast(i))
                Q[i] = -1;
            else
                Q[i] = 1;
        };

        WH = H.vector_left_product(Q);
        delete [] Q;

        int16_t max = 0, idx_max = -1;
        for(uint16_t i = 0; i < (1 << m); i++)
        {
            if(std::abs(WH[i]) > max)
            {
                max = std::abs(WH[i]);
                idx_max = i;
            };
        };
        message.setBitsValueFromQint64_fast(WH[idx_max] >= 0 ? (idx_max << 1) : ((idx_max << 1) | 1));
        delete [] WH;
    };
    return message;
}

SBitsArray RMCode::ns_FHT_decoder_fast2(SBitsArray &noisy_vector)
{
    SBitsArray message(m+1);
    if(r == 1)
    {
        hadamar_matrix_fast H(m, ns_MMatrices);
        int16_t* Q = new int16_t[1<<m], *WH;
        for(uint16_t i = 0; i < (1 << m); i++)
        {
            if(noisy_vector.getBit_fast(i))
                Q[i] = -1;
            else
                Q[i] = 1;
        };

        WH = H.vector_left_product(Q);
        delete [] Q;

        int16_t max = 0, idx_max = -1;
        for(uint16_t i = 0; i < (1 << m); i++)
        {
            if(std::abs(WH[i]) > max)
            {
                max = std::abs(WH[i]);
                idx_max = i;
            };
        };
        message.setBitsValueFromQint64_fast(WH[idx_max] >= 0 ? (idx_max << 1) : ((idx_max << 1) | 1));
        delete [] WH;
    };
    return message;
}

SBitsArray RMCode::FHT_decoder_ret_codeword(SBitsArray &noisy_vector)
{
    SBitsArray message = FHT_decoder(noisy_vector);
    return encode(message);
}

SBitsArray RMCode::FHT_decoder_ret_codeword_fast(SBitsArray &noisy_vector)
{
    SBitsArray message = FHT_decoder_fast(noisy_vector);
    return encode_fast(message);
}

SBitsArray RMCode::FHT_decoder_ret_codeword_fast2(SBitsArray &noisy_vector)
{
    SBitsArray message = FHT_decoder_fast2(noisy_vector);
    return encode_fast(message);
}

SBitsArray RMCode::ns_FHT_decoder_ret_codeword_fast2(SBitsArray &noisy_vector)
{
    SBitsArray message = ns_FHT_decoder_fast2(noisy_vector);
    return encode_fast(message);
}

SBitsArray RMCode::SSV_decoder(SBitsArray &noisy_vector)
{
    SBitsArray error(1<<m);
    bool skip_position;

    // вычислить r-syndrom
    SBitsArray r_s = r_syndrom(noisy_vector);
    make_base_lin_system(r_s);

    for(uint64_t v = 0; v < (1 << m); v++) // цикл по возможным позициям ошибок
    {
        SBitsMatrix lin_system(base_lin_system);
        b.setBitValue(0, true);
        skip_position = false;
        for(uint64_t M = 0; (M < E_dim()) && (!skip_position); M++) // цикл по мономам  степени не выше E_r
        {
            bool temp_b = eval_monom_on_vector(M, vs[v]);
            b.setBitValue(1 + (m+1)*M, temp_b);

            for(uint64_t l = 0; l < m; l++) // цикл по позициям
            {
                if(!(vs[v].getBit(l)))
                {
                    lin_system[1 + (m+1)*M + 1 + l].xorWithSBitArray(lin_system[1 + (m+1)*M]);
                }

                if((lin_system[1 + (m+1)*M + 1 + l].getHammingWeight() == 0) && (temp_b))
                {
                    skip_position = true;
                    break;
                }

                b.setBitValue(1 + (m+1)*M + 1 + l, temp_b);
            }
        };

        if(skip_position)
        {
            //cout << "skip" << endl;
            continue;
        }

        SBitsMatrix lin_system_ext(lin_system);
        lin_system_ext.append_column(b);

        uint64_t r1 = lin_system.low_triangular();
        uint64_t r2 = lin_system_ext.low_triangular();


        error.setBitValue(v, r1 == r2);
    }
    return error;
}

SBitsArray RMCode::SSV_decoder_fast(SBitsArray &noisy_vector)
{
    SBitsArray error(1<<m);
    bool skip_position;

    // вычислить r-syndrom
    SBitsArray r_s = r_syndrom_fast(noisy_vector);
    make_base_lin_system_fast(r_s);
    SBitsMatrix lin_system(base_lin_system.get_number_of_rows(), base_lin_system.get_number_of_cols());

    for(uint64_t v = 0; v < (1 << m); v++) // цикл по возможным позициям ошибок
    {
        lin_system.fast_copy(base_lin_system);
        b.setBitValue_fast(0, true);
        skip_position = false;
        for(uint64_t M = 0; (M < E_dim()) && (!skip_position); M++) // цикл по мономам степени не выше E_r
        {
            bool temp_b = eval_monom_on_vector_fast(M, vs[v]);
            b.setBitValue_fast(1 + (m+1)*M, temp_b);

            for(uint64_t l = 0; l < m; l++) // цикл по позициям
            {
                if(!(vs[v].getBit_fast(l)))
                {
                    lin_system[1 + (m+1)*M + 1 + l].xorWithSBitArray_fast(lin_system[1 + (m+1)*M]);
                }

                if((lin_system[1 + (m+1)*M + 1 + l].getHammingWeight_fast() == 0) && (temp_b))
                {
                    skip_position = true;
                    break;
                }

                b.setBitValue_fast(1 + (m+1)*M + 1 + l, temp_b);
            }
        };

        if(!skip_position)
        {
            SBitsMatrix lin_system_ext(lin_system.get_number_of_rows(), lin_system.get_number_of_cols());
            lin_system_ext.fast_copy(lin_system);
            lin_system_ext.append_column(b);
            uint64_t r2 = lin_system_ext.low_triangular_fast();
            error.setBitValue_fast(v, !((lin_system_ext[r2-1].getHammingWeight() == 1) && (lin_system_ext[r2-1].getBit_fast(lin_system.get_number_of_cols()) == true)));
        }
    }
    return error;
}

bool RMCode::SSV_decoder_fast_in_one_position(uint64_t v)
{
    bool skip_position;
    bool result = false;
    SBitsMatrix lin_system(base_lin_system.get_number_of_rows(), base_lin_system.get_number_of_cols());

    lin_system.fast_copy(base_lin_system);
    b.setBitValue_fast(0, true);
    skip_position = false;
    for(uint64_t M = 0; (M < E_dim()) && (!skip_position); M++) // цикл по мономам степени не выше E_r
    {
        bool temp_b = eval_monom_on_vector_fast(M, vs[v]);
        b.setBitValue_fast(1 + (m+1)*M, temp_b);

        for(uint64_t l = 0; l < m; l++) // цикл по позициям
        {
            if(!(vs[v].getBit_fast(l)))
            {
                lin_system[1 + (m+1)*M + 1 + l].xorWithSBitArray_fast(lin_system[1 + (m+1)*M]);
            }

            if((lin_system[1 + (m+1)*M + 1 + l].getHammingWeight_fast() == 0) && (temp_b))
            {
                skip_position = true;
                break;
            }

            b.setBitValue_fast(1 + (m+1)*M + 1 + l, temp_b);
        }
    };

    if(!skip_position)
    {
        SBitsMatrix lin_system_ext(lin_system.get_number_of_rows(), lin_system.get_number_of_cols());
        lin_system_ext.fast_copy(lin_system);
        lin_system_ext.append_column(b);
        uint64_t r2 = lin_system_ext.low_triangular_fast();
        result = !((lin_system_ext[r2-1].getHammingWeight() == 1) && (lin_system_ext[r2-1].getBit_fast(lin_system.get_number_of_cols()) == true));
    }
    return result;
}

void RMCode::init_factor_classes_for_RPA(uint16_t max_m)
{
    RMCode::fc = new FactorClasses[max_m + 1];
    for(uint16_t i = 1; i <= max_m; i++)
    {
        //cout << "i=" << i << endl;
        (RMCode::fc)[i].init(i);
        (RMCode::fc)[i].make_factor_classes();
        //(RMCode::fc)[i].print_factor_classes();
    }
}

void RMCode::deinit_factor_classes_for_RPA()
{
    delete [] RMCode::fc;
    RMCode::fc = nullptr;
}

void RMCode::init_M_matrices_for_Hadamar(uint16_t max_m)
{
    MatrixInt16 H1(2,2);
    H1.set(0,0,1);
    H1.set(0,1,1);
    H1.set(1,0,1);
    H1.set(1,1,-1);

    m_limit = max_m;

    MMatrices = new MatrixInt16**[max_m + 1];

    for(uint16_t m_ = 1; m_ <= max_m; m_++)
    {
        MMatrices[m_] = new MatrixInt16*[m_];
        memset(MMatrices[m_], 0, sizeof(MatrixInt16*)*m_);

        for(uint16_t j = 0; j < m_; j++)
        {
            MatrixInt16 Il(1<<(m_-(j+1)), 1<<(m_-(j+1))), Ir(1<<j, 1<<j), *temp;
            Il.make_identity_matrix();
            Ir.make_identity_matrix();
            temp = Il.tensor_product(H1);
            MMatrices[m_][j] = temp->tensor_product(Ir);
            MMatrices[m_][j]->make_compact();
            delete temp;
        }
    }
}

void RMCode::deinit_M_matrices_for_Hadamar()
{

    for(uint16_t m_ = 1; m_ <= m_limit; m_++)
    {
        for(uint16_t j = 0; j < m_; j++)
        {
            delete MMatrices[m_][j];
        };
        delete [] MMatrices[m_];
    };
    delete [] MMatrices;
}

void RMCode::init_codes_for_RPA_decoder(uint16_t max_r, uint16_t max_m)
{
    if(max_r > 1)
    {
        codes_for_RPA_decoding = (RMCode**)malloc(sizeof(RMCode**) * max_r-1);
        for(int i = 1; i < max_r - 1; i++)
            codes_for_RPA_decoding[i - 1] = new RMCode(max_r - i, max_m - i);
    }
}

RMCode* RMCode::get_code_for_RPA_decoder(uint16_t val_r, uint16_t max_r)
{
    return codes_for_RPA_decoding[(max_r - 1) - val_r];
}

SBitsArray RMCode::RPA_decoder(SBitsArray &y, uint16_t N_max)
{
    SBitsArray c(y);
    if(r == 1)
    {
        // если код первого порядка, то используем декодер с помощью быстрого преобразования Адамара
        c = FHT_decoder_ret_codeword_fast(y);
    }else{
        // иначе используем рекурсию
        RMCode rm_small(r-1, m-1);
        for(uint16_t i = 0; i < N_max; i++)
        {
            uint16_t* changevote = new uint16_t[1<<m];
            memset(changevote, 0, sizeof(uint16_t) * (1<<m));
            for(uint64_t z0 = 1; z0 < (1<<m); z0++)
            {
                SBitsArray B((unsigned char*)(&z0), 1<<m);
                SBitsArray y_B = RMCode::fc[m].Proj(c, B);
                SBitsArray hat_y_B = rm_small.RPA_decoder(y_B, N_max);
                uint16_t *factor_classes = RMCode::fc[m].get_factor_classes_idx_for_Bi(B);

                // цикл по классам смежности
                for(uint64_t z = 0; z < (1<<(m-1)); z++)
                {
                    if(y_B.getBit_fast(z) != hat_y_B.getBit_fast(z))
                    {
                        // для всех векторов текущего класса смежности увеличить счетчик
                        changevote[factor_classes[2*z]]++;
                        changevote[factor_classes[2*z+1]]++;
                    };
                }
            }

            uint16_t numberofchange = 0;
            uint64_t n = 1<<m;

            for(uint64_t z = 0; z < 1<<m; z++)
            {
                if(changevote[z] > ((n-1)/2.))
                {
                    c.invertBit_fast(z);
                    numberofchange++;
                }
            };

            delete[]changevote;

            if(numberofchange == 0)
            {
                break;
            }
        };
    };
    return c;
}

SBitsArray RMCode::RPA_decoder_fast(SBitsArray &y, uint16_t N_max)
{
    SBitsArray c(y);
    if(r == 1)
    {
        // если код первого порядка, то используем декодер с помощью быстрого преобразования Адамара
        c = FHT_decoder_ret_codeword_fast(y);
    }else{
        // иначе используем рекурсию
        RMCode rm_small(r-1, m-1);
        for(uint16_t i = 0; i < N_max; i++)
        {
            uint16_t *changevote = new uint16_t[1<<m];
            memset(changevote, 0, sizeof(uint16_t) * (1<<m));
            for(uint64_t z0 = 1; z0 < (1<<m); z0++)
            {
                SBitsArray B((unsigned char*)(&z0), 1<<m);
                SBitsArray y_B = RMCode::fc[m].Proj_fast(c, B);
                SBitsArray hat_y_B = rm_small.RPA_decoder_fast(y_B, N_max);
                uint16_t *factor_classes = RMCode::fc[m].get_factor_classes_idx_for_Bi(B);

                // цикл по классам смежности
                for(uint64_t z = 0; z < (1<<(m-1)); z++)
                {
                    if(y_B.getBit_fast(z) != hat_y_B.getBit_fast(z))
                    {
                        // для всех векторов текущего класса смежности увеличить счетчик
                        changevote[factor_classes[2*z]]++;
                        changevote[factor_classes[2*z+1]]++;
                    };
                }
            }

            uint16_t numberofchange = 0;
            uint64_t n = 1<<m;

            for(uint64_t z = 0; z < 1<<m; z++)
            {
                if(changevote[z] > ((n-1)/2.))
                {
                    c.invertBit_fast(z);
                    numberofchange++;
                }
            };

            delete[] changevote;
            if(numberofchange == 0)
            {
                break;
            }
        };
    };
    return c;
}

SBitsArray RMCode::ns_RPA_decoder_fast2(SBitsArray &y, uint16_t N_max)
{
    SBitsArray c(y);
    if(r == 1)
    {
        // если код первого порядка, то используем декодер с помощью быстрого преобразования Адамара
        c = ns_FHT_decoder_ret_codeword_fast2(y);
    }else{
        // иначе используем рекурсию
        RMCode rm_small(r-1, m-1);
        rm_small.ns_fc = ns_fc;
        rm_small.ns_MMatrices = ns_MMatrices;
        rm_small.ns_m_limit = ns_m_limit;

        for(uint16_t i = 0; i < N_max; i++)
        {
            uint16_t* changevote = new uint16_t[1<<m];
            memset(changevote, 0, sizeof(uint16_t) * (1<<m));
            for(uint64_t z0 = 1; z0 < (1<<m); z0++)
            {
                SBitsArray B((unsigned char*)(&z0), 1<<m);
                SBitsArray y_B = ns_fc[m].Proj_fast(c, B);
                SBitsArray hat_y_B = rm_small.ns_RPA_decoder_fast2(y_B, N_max);
                uint16_t *factor_classes = ns_fc[m].get_factor_classes_idx_for_Bi(B);

                // цикл по классам смежности
                for(uint64_t z = 0; z < (1<<(m-1)); z++)
                {
                    if(y_B.getBit_fast(z) != hat_y_B.getBit_fast(z))
                    {
                        // для всех векторов текущего класса смежности увеличить счетчик
                        changevote[factor_classes[2*z]]++;
                        changevote[factor_classes[2*z+1]]++;
                    };
                }
            }

            uint16_t numberofchange = 0;
            uint64_t n = 1<<m;

            for(uint64_t z = 0; z < 1<<m; z++)
            {
                if(changevote[z] > ((n-1)/2.))
                {
                    c.invertBit_fast(z);
                    numberofchange++;
                }
            };
            delete[] changevote;
            if(numberofchange == 0)
            {
                break;
            }
        };
    };
    return c;
}

SBitsArray RMCode::ns_RPA_decoder_fast2_with_common_data(SBitsArray &y, uint16_t N_max, uint16_t val_r, uint16_t max_r)
{
    SBitsArray c(y);
    if(r == 1)
    {
        // если код первого порядка, то используем декодер с помощью быстрого преобразования Адамара
        c = ns_FHT_decoder_ret_codeword_fast2(y);
    }else{
        // иначе используем рекурсию
        RMCode rm_small(r-1, m-1, get_code_for_RPA_decoder(val_r-1, max_r));
        rm_small.ns_fc = ns_fc;
        rm_small.ns_MMatrices = ns_MMatrices;
        rm_small.ns_m_limit = ns_m_limit;

        for(uint16_t i = 0; i < N_max; i++)
        {
            uint16_t* changevote = new uint16_t[1<<m];
            memset(changevote, 0, sizeof(uint16_t) * (1<<m));
            for(uint64_t z0 = 1; z0 < (1<<m); z0++)
            {
                SBitsArray B((unsigned char*)(&z0), 1<<m);
                SBitsArray y_B = ns_fc[m].Proj_fast(c, B);
                SBitsArray hat_y_B = rm_small.ns_RPA_decoder_fast2_with_common_data(y_B, N_max, val_r-1, max_r);
                uint16_t *factor_classes = ns_fc[m].get_factor_classes_idx_for_Bi(B);

                // цикл по классам смежности
                for(uint64_t z = 0; z < (1<<(m-1)); z++)
                {
                    if(y_B.getBit_fast(z) != hat_y_B.getBit_fast(z))
                    {
                        // для всех векторов текущего класса смежности увеличить счетчик
                        changevote[factor_classes[2*z]]++;
                        changevote[factor_classes[2*z+1]]++;
                    };
                }
            }

            uint16_t numberofchange = 0;
            uint64_t n = 1<<m;

            for(uint64_t z = 0; z < 1<<m; z++)
            {
                if(changevote[z] > ((n-1)/2.))
                {
                    c.invertBit_fast(z);
                    numberofchange++;
                }
            };

            delete[] changevote;

            if(numberofchange == 0)
            {
                break;
            }
        };
    };
    return c;
}

SBitsArray RMCode::RPA_decoder_fast2(SBitsArray &y, uint16_t N_max)
{
    SBitsArray c(y);
    if(r == 1)
    {
        // если код первого порядка, то используем декодер с помощью быстрого преобразования Адамара
        c = FHT_decoder_ret_codeword_fast2(y);
    }else{
        // иначе используем рекурсию
        RMCode rm_small(r-1, m-1);
        for(uint16_t i = 0; i < N_max; i++)
        {
            uint16_t* changevote = new uint16_t[1<<m];
            memset(changevote, 0, sizeof(uint16_t) * (1<<m));
            for(uint64_t z0 = 1; z0 < (1<<m); z0++)
            {
                SBitsArray B((unsigned char*)(&z0), 1<<m);
                SBitsArray y_B = RMCode::fc[m].Proj_fast(c, B);
                SBitsArray hat_y_B = rm_small.RPA_decoder_fast2(y_B, N_max);
                uint16_t *factor_classes = RMCode::fc[m].get_factor_classes_idx_for_Bi(B);

                // цикл по классам смежности
                for(uint64_t z = 0; z < (1<<(m-1)); z++)
                {
                    if(y_B.getBit_fast(z) != hat_y_B.getBit_fast(z))
                    {
                        // для всех векторов текущего класса смежности увеличить счетчик
                        changevote[factor_classes[2*z]]++;
                        changevote[factor_classes[2*z+1]]++;
                    };
                }
            }

            uint16_t numberofchange = 0;
            uint64_t n = 1<<m;

            for(uint64_t z = 0; z < 1<<m; z++)
            {
                if(changevote[z] > ((n-1)/2.))
                {
                    c.invertBit_fast(z);
                    numberofchange++;
                }
            };

            delete[] changevote;

            if(numberofchange == 0)
            {
                break;
            }
        };
    };
    return c;
}

void RMCode::test_on_one_position()
{
    uint8_t r = 2, m = 14, t = 30;
    RMCode RM(m-(2*r+2),m);
    SBitsArray info_vector(RM.dim());
    SBitsArray code_vector(1<<m);
    SBitsArray error_vector(1<<m);
    SBitsArray corrupted_vector(1<<m);
    SBitsArray r_syndrom((1<<m) - RM.dim());

    info_vector.randomBitArray();
    cout << "info generated" << endl;
    code_vector = RM.encode(info_vector);
    cout << "code generated" << endl;
    error_vector.randomBitArrayWixedWeight(t);
    cout << "error generated" << endl;
    corrupted_vector = SBitsArray::xorWithSBitArray_fast(code_vector, error_vector);
    cout << "corrupted vector computed" << endl;
    r_syndrom = RM.r_syndrom_fast(corrupted_vector);
    cout << "r-syndrom computed" << endl;
    RM.make_base_lin_system_fast(r_syndrom);
    cout << "base linear system generated" << endl;

    clock_t t1 = clock();
    RM.SSV_decoder_fast_in_one_position(0);
    cout << "decoding time in one position is " << (clock()-t1) << endl;
}

void RMCode::test_on_all_positions(float bsc_p, uint8_t proc)
{
    uint8_t r = 2, m = 8, t = 31;
    RMCode RM(m-(2*r+2),m);
    SBitsArray info_vector(RM.dim());
    SBitsArray code_vector(1<<m);
    SBitsArray error_vector(1<<m);
    SBitsArray corrupted_vector(1<<m);
    SBitsArray r_syndrom((1<<m) - RM.dim());
    uint32_t counter = 0, N = 10;

    info_vector.randomBitArray();
    //cout << "info " << info_vector << endl;
    code_vector = RM.encode(info_vector);
    //cout << "code " << code_vector << endl;

    for(uint32_t i = 0; i < N; i++)
    {
        error_vector.pass_through_bsc(bsc_p, time(NULL) + i + proc);
                //randomBitArrayWixedWeight(t);
        corrupted_vector = SBitsArray::xorWithSBitArray_fast(code_vector, error_vector);
        //clock_t t1 = clock();
        SBitsArray ret = RM.SSV_decoder_fast(corrupted_vector);
        //cout << "time=" << (clock()-t1) << endl;
        ret.xorWithSBitArray_fast(error_vector);
        if(ret.getHammingWeight_fast() != 0)
            counter++;
        if( (i%1000) == 0)
            cout << i << "\r";
    };
    cout << endl << "p=" << bsc_p << ", DFR=" << counter/float(N) << endl;
    getchar();
}

void RMCode::test_FHT(float bsc_p, uint8_t proc)
{
    uint8_t r = 1, m = 8, t = 85;
    RMCode RM(r,m);
    SBitsArray info_vector(RM.dim());
    SBitsArray code_vector(1<<m);
    SBitsArray error_vector(1<<m);
    SBitsArray corrupted_vector(1<<m);
    uint32_t counter = 0, N = 100000;

    //RM.print_gen_matrix();

    for(uint32_t i = 0; i < N; i++)
    {
        info_vector.randomBitArray();
        code_vector = RM.encode(info_vector);
        error_vector.randomBitArrayWixedWeight(t);
        corrupted_vector = SBitsArray::xorWithSBitArray_fast(code_vector, error_vector);
        SBitsArray ret = RM.FHT_decoder_ret_codeword(corrupted_vector);

        ret.xorWithSBitArray(code_vector);
        if(ret.getHammingWeight() != 0)
            counter++;
        //cout << info_vector << " " << ret << endl;

    };
    cout << "DFR=" << counter/(1.*N) << endl;
}

void RMCode::test_RPA(uint8_t r, uint8_t m, uint32_t N, float bsc_p)
{
    RMCode RM(r,m);
    SBitsArray info_vector(RM.dim());
    SBitsArray code_vector(1<<m);
    SBitsArray error_vector(1<<m);
    SBitsArray corrupted_vector(1<<m);
    uint32_t counter = 0;
    uint32_t avg_weight = 0;

    //RM.print_gen_matrix();

    cout << "start making factor classes" << endl;
    RMCode::init_factor_classes_for_RPA(m);
    cout << "stop making factor classes" << endl;

    cout << "average error weight " << uint16_t(bsc_p*(1<<m)) << endl;
    cout << "start decoding" << endl;
    for(uint32_t i = 0; i < N; i++)
    {
        info_vector.randomBitArray();
        code_vector = RM.encode_fast(info_vector);
        //error_vector.toZeroAllBits();
        //error_vector.pass_through_bsc_secure_rand_isaak(bsc_p, clock());
        error_vector.randomBitArrayWixedWeight_fast(uint16_t(bsc_p*(1<<m)));
        corrupted_vector = SBitsArray::xorWithSBitArray_fast(code_vector, error_vector);
        SBitsArray ret = RM.RPA_decoder(corrupted_vector, m/2);
        ret.xorWithSBitArray(code_vector);
        if(ret.getHammingWeight_fast() != 0)
            counter++;
        cout << i+1 << "/" << N << ", current error weight=" << error_vector.getHammingWeight_fast() << ", real average weight " << uint32_t(avg_weight*1./(i+1)) << ", current DFR=" << (counter*1.0)/(i+1) << "\r";
        avg_weight += error_vector.getHammingWeight_fast();
    };
    cout << "\nnumber of decoding failures = " << counter << ", real average weight " << avg_weight*1./N << endl;
    RMCode::deinit_factor_classes_for_RPA();
}

void RMCode::test_RPA_by_blocks(uint8_t r, uint8_t m, uint16_t num_blocks, uint16_t weight, uint32_t N)
{
    RMCode RM(r,m);
    SBitsArray info_vector(RM.dim());
    SBitsArray* code_vectors = new SBitsArray[num_blocks];
    SBitsArray corrupted_code_vector((1<<m));
    SBitsArray full_code_vector((1<<m)*num_blocks);
    SBitsArray full_error_vector((1<<m)*num_blocks);
    SBitsArray full_corrupted_vector((1<<m)*num_blocks);
    SBitsArray full_decoded_vector((1<<m)*num_blocks);
    uint32_t* counters = new uint32_t[num_blocks+1];
    uint32_t number_of_wrong_codewords = 0;
    //uint32_t number_of_right_decryptions = 0;

    // initialize counters for decoding failures
    for(uint16_t i = 0; i < num_blocks+1; i++)
    {
        counters[i] = 0;
    };

    for(uint16_t i = 0; i < num_blocks; i++)
    {
        code_vectors[i].allocateBufferForBits(1<<m);
    };

    cout << "r=" << uint16_t(r) << ", m=" << uint16_t(m) << ", num_blocks=" << num_blocks << ", weight=" << weight << endl;
    cout << "start making factor classes" << endl;
    RMCode::init_factor_classes_for_RPA(m);
    cout << "stop making factor classes" << endl;

    cout << "start making Hadamar matrices" << endl;
    RMCode::init_M_matrices_for_Hadamar(m);
    cout << "stop making Hadamar matrices" << endl;

    //cout << "start making permutation" << endl;
    //Permutation perm((1<<m)*num_blocks, clock());
    //cout << "stop making permutation" << endl;

    cout << "start decoding at " << time(NULL) << " s." << endl;
    for(uint32_t it = 0; it < N; it++)
    {
        Permutation perm((1<<m)*num_blocks, clock());

        uint64_t t1 = time(NULL);
        uint64_t t1_clock = clock();
        // concatenate blocks
        for(uint16_t j = 0; j < num_blocks; j++)
        {
            info_vector.randomBitArray();
            code_vectors[j] = RM.encode_fast(info_vector);
            //cout << "code_vectors[j]\t\t\t\t" << code_vectors[j] << endl;
            full_code_vector.setSubvectorInPosition_fast(j*(1<<m), code_vectors[j]);
        };

        //cout << "full_code_vector\t\t\t" << full_code_vector << endl;
        // permute code vector
        full_code_vector.permute(perm);
        //cout << "permuted full_code_vector\t\t" << full_code_vector << endl;

        // add random error
        full_error_vector.randomBitArrayWixedWeight_fast(weight);
        full_corrupted_vector = SBitsArray::xorWithSBitArray_fast(full_code_vector, full_error_vector);
        //cout << "full_corrupted_vector\t\t\t" << full_corrupted_vector << endl;

        // inverse permute corrupted vector
        full_corrupted_vector.inverse_permute(perm);
        //cout << "inverse permuted full_corrupted_vector\t" << full_corrupted_vector << endl;

        // decode each codeword separately and count number of decoding failures
        uint16_t local_counter = 0;


        //QFuture<void> future[num_blocks];

        uint64_t t2_clock = clock();
        uint64_t distance_to_decoded = 0;
        for(uint16_t j = 0; j < num_blocks; j++)
        {
            corrupted_code_vector = full_corrupted_vector.getSubvectorInPosition_fast(j*(1<<m), (1<<m));
            //cout << "corrupted_code_vector\t\t\t" << corrupted_code_vector << endl;
            //SBitsArray ret = RM.RPA_decoder_fast2(corrupted_code_vector, m/2);
            SBitsArray ret = RM.RPA_decoder_fast2(corrupted_code_vector, (m%2 == 0) ? m/2 : m/2+1);
            //SBitsArray ret_copy(ret);
            //ret_copy.xorWithSBitArray_fast(corrupted_code_vector);
            //distance_to_decoded += ret_copy.getHammingWeight_fast();

            SBitsArray synd = RM.r_syndrom_fast(ret);
            ret.xorWithSBitArray_fast(code_vectors[j]);
            if(ret.getHammingWeight_fast() != 0)
            {
                local_counter++;
                if(synd.getHammingWeight_fast() == 0)
                    number_of_wrong_codewords++;
            }
        };
        uint64_t t3_clock = clock();

        /*
        if(distance_to_decoded == weight)
            number_of_right_decryptions++;
        */

        // count number of decoding failures
        counters[local_counter]++;
        cout << it+1 << "/" << N << " DFR counters ";
        for(uint16_t i = 0; i < num_blocks+1; i++)
        {
            cout << i << ":" << counters[i] << ",";
        };
        cout << " iter_time=" << time(NULL)-t1 << "s., dec_clocks=" << t3_clock-t2_clock << ", total_clocks=" << clock()-t1_clock << " \r";
    };

    cout << "\nstop decoding at " << time(NULL) << " s., number_of_wrong_codewords=" << number_of_wrong_codewords << /*", number_of_right_decryptions_by_weight=" << number_of_right_decryptions <<*/ endl;
    RMCode::deinit_factor_classes_for_RPA();
    RMCode::deinit_M_matrices_for_Hadamar();
}

SBitsArray RMCode::linear_automotphism(SBitsArray &vector, SBitsMatrix& la)
{
    SBitsArray permuted_vec(vector);

    for(uint64_t z = 0; z < 1<<m; z++)
    {
        SBitsArray z_vec((unsigned char*) (&z), m);
        SBitsArray res = la.mul_left(z_vec);
        uint64_t idx = res.getLowestBitsAsInteger_fast(m);
        permuted_vec.setBitValue_fast(z, vector.getBit_fast(idx));
    };

    return permuted_vec;
}
