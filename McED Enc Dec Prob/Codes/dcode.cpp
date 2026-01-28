
#include "dcode.h"

void DCode::add_tensor_product(quint16 r1, quint16 m1, quint16 r2, quint16 m2)
{
    RMCode rm1(r1, m1, false);
    RMCode rm2(r2, m2, false);
    quint64 rank = 0;

    SBitsMatrix *gen1= rm1.get_gen_matrix();
    SBitsMatrix *gen2= rm2.get_gen_matrix();
    SBitsMatrix *tensor = SBitsMatrix::tensor_product(*gen1, *gen2);

    if(gen_matrix == nullptr)
    {
        gen_matrix = new SBitsMatrix(*tensor);
    }else{
        gen_matrix->concat_rows(*tensor);
    };

    rank = gen_matrix->low_triangular_fast2();
    if(rank != gen_matrix->get_number_of_rows())
    {
        gen_matrix->delete_last_rows(gen_matrix->get_number_of_rows() - rank);
    };

    rm_m1 = m1;
    rm_m2 = m2;

    delete gen1;
    delete gen2;
    delete tensor;
}

SBitsMatrix* DCode::get_K_random_subblocks(quint64 K)
{
    quint64 nr = gen_matrix->get_number_of_rows();
    SBitsMatrix *matrix = new SBitsMatrix(nr, (1<<rm_m2)*K);

    std::srand (clock());
    std::vector<int> myvector;

    for(uint16_t i = 0; i < (1<<rm_m1); i++)
    {
        myvector.push_back(i);
    };
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(myvector.begin(), myvector.end(), g);
    sort(myvector.begin(), myvector.begin() + K);
    //sort(&(myvector[0]), &(myvector[0]) + K);

    // get fist K

    for(uint16_t r = 0; r < nr; r++)
        for(uint16_t i = 0; i < K; i++)
        {
            SBitsArray temp = (*gen_matrix)[r].getSubvectorInPosition_fast(myvector[i]*(1<<rm_m2), (1<<rm_m2));
            (*matrix)[r].setSubvectorInPosition_fast(i*(1<<rm_m2), temp);
        };

    return matrix;
}

SBitsMatrix* DCode::get_subblocks(std::vector<int> &myvector)
{
    quint64 nr = gen_matrix->get_number_of_rows();
    SBitsMatrix *matrix = new SBitsMatrix(nr, (1<<rm_m2)*myvector.size());

    for(uint16_t r = 0; r < nr; r++)
    {
        for(uint16_t i = 0; i < myvector.size(); i++)
        {
            SBitsArray temp = (*gen_matrix)[r].getSubvectorInPosition_fast(myvector[i]*(1<<rm_m2), (1<<rm_m2));
            (*matrix)[r].setSubvectorInPosition_fast(i*(1<<rm_m2), temp);
        };
    }

    return matrix;
}

SBitsArray DCode::isd_one_iteration(quint64 K, SBitsArray &z, bool *full_rank)
{
    quint64 dimension = get_dimension();
    SBitsMatrix *id_matrix = SBitsMatrix::identity_matrix(dimension);
    quint64 *column_indices = new quint64[dimension];

    std::srand (clock());
    std::vector<int> block_indices;

    for(uint16_t i = 0; i < (1<<rm_m1); i++)
    {
        block_indices.push_back(i);
    };

    //std::random_shuffle(block_indices.begin(), block_indices.end());
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(block_indices.begin(), block_indices.end(), g);
    std::vector<int> sub_indices(block_indices.begin(), block_indices.begin() + K);
    //std::vector<int> sub_indices(&(block_indices[0]), &(block_indices[0]) + K);
    sort(sub_indices.begin(), sub_indices.end());
    cout << "random choosen blocks: ";
    for (int k_rnd = 0; k_rnd < sub_indices.size(); k_rnd++)
        cout << sub_indices[k_rnd] << " ";
    SBitsMatrix *proj = get_subblocks(sub_indices);
    proj->concat_columns(*id_matrix);
    delete id_matrix;

    proj->pseudo_systematic_form();

    SBitsArray decoded(dimension);
    if((*proj)[dimension-1].getIndexOfFirstNonzero() < (proj->get_number_of_cols() - dimension))
    {
        SBitsMatrix *inverse = proj->get_last_n_columns(dimension);
        for(uint16_t k = 0; k < dimension; k++)
        {
            qint64 nz = (*proj)[k].getIndexOfFirstNonzero();
            column_indices[k] = sub_indices[nz/(1<<rm_m2)]*(1<<rm_m2) + nz%(1<<rm_m2);
        };
        SBitsArray z_proj(z.getProjection_fast(column_indices, dimension));
        decoded = inverse->mul_left(z_proj);
        delete inverse;
        *full_rank = true;
    }else{
        *full_rank = false;
    };

    delete proj;
    delete [] column_indices;

    return decoded;
}

void DCode::test_EncryptionDecryption_ISD(const char* f_params, uint16_t ISD_iter, uint32_t N)
{
    std::ifstream dcode_file(f_params);
    if(dcode_file.is_open())
    {
        DCode dcode;
        uint8_t r, m;
        uint16_t weight;
        dcode_file >> weight;
        uint16_t K;
        dcode_file >> K;
        uint8_t iter = 0;
        uint64_t total_block_decoding_clocks = 0;
        uint64_t total_block_decoding_iters = 0;
        uint64_t total_isd_decoding_clocks = 0;
        uint64_t total_isd_decoding_iters = 0;

        int r1, m1, r2, m2;
        while (dcode_file >> r1 >> m1 >> r2 >> m2)
        {
            dcode.add_tensor_product(r1, m1, r2, m2);
            if(iter == 0)
            {
                r = r2;
                m = m2;
            };

            iter++;
        }
        dcode.print_params();

        RMCode RM(r,m,false);
        SBitsArray corrupted_code_vector((1<<m));
        SBitsArray info_vector(dcode.get_dimension());
        SBitsArray info_vector_dec(dcode.get_dimension());
        SBitsArray full_code_vector((1<<m)*(1<<dcode.rm_m1));
        SBitsArray full_error_vector((1<<m)*(1<<dcode.rm_m1));
        SBitsArray full_corrupted_vector((1<<m)*(1<<dcode.rm_m1));
        SBitsArray full_decoded_vector((1<<m)*(1<<dcode.rm_m1));
        uint32_t DFR = 0;
        uint32_t wrong_check = 0;
        SBitsMatrix *S_inv = nullptr;
        SBitsMatrix S_copy(dcode.get_dimension(), dcode.get_dimension());
        SBitsMatrix *identity = SBitsMatrix::identity_matrix(dcode.get_dimension());

        cout << "r=" << uint16_t(r) << ", m=" << uint16_t(m) << ", num_blocks=" << (1<<dcode.rm_m1) << ", weight=" << weight << endl;

        cout << "Generator matrix" << endl;
        dcode.gen_matrix->print_nonzero_rows();

        cout << "start making factor classes" << endl;
        RMCode::init_factor_classes_for_RPA(m);
        cout << "stop making factor classes" << endl;

        cout << "start making Hadamar matrices" << endl;
        RMCode::init_M_matrices_for_Hadamar(m);
        cout << "stop making Hadamar matrices" << endl;


        cout << "start making information set for nonzero vector" << endl;
        SBitsMatrix *id_matrix = SBitsMatrix::identity_matrix(dcode.get_dimension());
        SBitsMatrix gen_matr_copy(dcode.gen_matrix->get_number_of_rows(), dcode.gen_matrix->get_number_of_cols());
        gen_matr_copy.fast_copy(*(dcode.gen_matrix));
        gen_matr_copy.pseudo_systematic_form();
        quint64* information_set_base = new quint64[dcode.get_dimension()];
        for(uint64_t is = 0; is < dcode.get_dimension(); is++)
            information_set_base[is] = gen_matr_copy[is].getIndexOfFirstNonzero();
        cout << "stop making information set for nonzero vector" << endl;

        uint64_t t_total_keygen = 0;
        uint64_t t_total_encoding = 0;
        uint64_t t_total_decoding = 0;
        uint64_t t_total_ISD_one_attempt_end = 0;
        uint64_t rnd_seed = 0;
        for(uint32_t it = 0; it < N; it++)
        {
            SBitsMatrix gen_matr_copy2(dcode.gen_matrix->get_number_of_rows(), dcode.gen_matrix->get_number_of_cols());
            cout << "iter#" << it+1 << "(start):" << endl;
            // Key generation
            uint64_t t_keygen_begin = clock();
            Permutation perm((1<<dcode.rm_m1)*((1<<dcode.rm_m2)), clock());
            cout << "Column permutation (aka matrix P):" << endl;
            for (int p = 0; p < (1 << dcode.rm_m1) * (1 << dcode.rm_m2); p++)
                cout << perm.get(p) << ' ';
            cout << endl;
            do{
                rnd_seed++;
                //cout << "1" << endl;
                SBitsMatrix S(dcode.get_dimension(), dcode.get_dimension());
                //cout << "2" << endl;
                S.random(it + time(NULL) + rnd_seed);
                S_copy.fast_copy(S);
                //cout << "3" << endl;
                if(S.low_triangular_fast2() == dcode.get_dimension())
                {
                    //cout << "3.1" << endl;
                    S.fast_copy(S_copy);
                    //cout << "3.2" << endl;
                    S.concat_columns(*identity);
                    //cout << "3.3" << endl;
                    S.pseudo_systematic_form();
                    //cout << "3.4" << endl;
                    S_inv = S.get_last_n_columns(dcode.get_dimension());
                    //cout << "3.5" << endl;
                    break;
                };
                //cout << "4" << endl;
            }while(1);
            cout << "matrix S:" << endl;
            S_copy.print_nonzero_rows();
            SBitsMatrix* semi_pub_key = nullptr;
            SBitsMatrix *pub_key = nullptr;
            semi_pub_key = dcode.gen_matrix->mul_left_matrix(S_copy);
            pub_key = semi_pub_key->permute_columns(perm);
            

            uint64_t t_keygen_end = clock() - t_keygen_begin;
            cout << "Public key (G_pub):" << endl;
            pub_key->print_nonzero_rows();

            // Encryption
            info_vector.randomBitArray();
            cout << "Message (m):" << endl;
            cout << info_vector << endl;

            uint64_t t_encoding_begin = clock();
            full_code_vector = pub_key->mul_left(info_vector);
            cout << "Data encoded using public matrix (m*G_pub):" << endl;
            cout << full_code_vector << endl;

            full_error_vector.randomBitArrayWixedWeight_fast(weight);
            cout << "Error vector (e):" << endl;
            cout << full_error_vector << endl;
            full_corrupted_vector = SBitsArray::xorWithSBitArray_fast(full_code_vector, full_error_vector);
            uint64_t t_encoding_end = clock() - t_encoding_begin;
            cout << "Ciphertext (c=mG_pub+e):" << endl;
            cout << full_corrupted_vector << endl;

            // Decryption
            uint64_t t_total_decoding_begin = clock();
            full_corrupted_vector.inverse_permute(perm);
            cout << "Start decryption:" << endl << "Apply inverse permutation (c*P^{-1}):" << endl;
            cout << full_corrupted_vector << endl;

            uint64_t t_block_decoding_begin = clock();
            for(uint16_t j = 0; j < (1<<dcode.rm_m1); j++)
            {
                corrupted_code_vector = full_corrupted_vector.getSubvectorInPosition_fast(j*(1<<m), (1<<m));
                SBitsArray ret = RM.RPA_decoder_fast2(corrupted_code_vector, (m%2 == 0) ? m/2 : m/2+1);
                full_decoded_vector.setSubvectorInPosition_fast(j*(1<<m), ret);
            };
            uint64_t t_block_decoding_end = clock() - t_block_decoding_begin;

            total_block_decoding_clocks += t_block_decoding_end;
            total_block_decoding_iters++;

            cout << "Noisy codeword decoded by blocks (block_decoder(c*P^{-1})):" << endl;
            cout << full_decoded_vector << endl;

            // First attemp

            uint64_t t_ISD_one_attempt_begin= clock();
            SBitsArray z_proj(full_decoded_vector.getProjection_fast(information_set_base, dcode.get_dimension()));
            gen_matr_copy2.fast_copy(*(semi_pub_key));
            gen_matr_copy2.concat_columns(*id_matrix);
            gen_matr_copy2.pseudo_systematic_form();
            SBitsMatrix *inverse_base = gen_matr_copy2.get_last_n_columns(dcode.get_dimension());
            info_vector_dec = inverse_base->mul_left(z_proj);
            delete inverse_base;

            // Checking
            SBitsArray c(dcode.get_length());
            c = semi_pub_key->mul_left(info_vector_dec);
            c.xorWithSBitArray_fast(full_corrupted_vector);
            uint64_t t_ISD_one_attempt_end = clock() - t_ISD_one_attempt_begin;

            total_isd_decoding_clocks += t_ISD_one_attempt_end;
            total_isd_decoding_iters++;

            if(c.getHammingWeight_fast() != weight)
            {
                // Second attempt
                uint16_t isd = 0;
                bool full_rank = true;

                ///////////////////////////////////////////////////////////
                // ISDDecoder (begin)
                while(isd < ISD_iter)
                {
                    uint64_t t_IDS_decoding_begin = clock();
                    cout << "\t\tISD: " << isd << " ";
                    info_vector_dec = dcode.isd_one_iteration(K, full_decoded_vector, &full_rank);
                    if(!full_rank)
                    {
                        cout << ", not full rank " << endl;
                        continue;
                    }
                    info_vector_dec = S_inv->mul_left(info_vector_dec);
                    cout << " ISD-decoded message: " << info_vector_dec << endl;
                    // Checking
                    SBitsArray c(dcode.get_length());
                    c = semi_pub_key->mul_left(info_vector_dec);
                    c.xorWithSBitArray_fast(full_corrupted_vector);

                    __int64 fix_clocks = clock() - t_IDS_decoding_begin;

                    //cout << "\t\tISD: " << isd << ", clocks=" << fix_clocks << endl;

                    total_isd_decoding_clocks += fix_clocks;
                    total_isd_decoding_iters++;

                    if(c.getHammingWeight_fast() == weight)
                    {
                        info_vector_dec.xorWithSBitArray_fast(info_vector);
                        if (info_vector_dec.getHammingWeight_fast() != 0)
                        {
                            cout << "Decryption failure." << endl;
                            wrong_check++;
                        }
                        else {
                            cout << "Decryption success." << endl;
                        };
                        break;
                    };
                    isd++;
                };
                // ISDDecoder (end)
                ///////////////////////////////////////////////////////////

                if(isd == ISD_iter)
                {
                    cout << "Decryption failure (due limit for ISD iteration)." << endl;
                    DFR++;
                };
            }else{
                info_vector_dec.xorWithSBitArray_fast(info_vector);
                if (info_vector_dec.getHammingWeight_fast() != 0)
                {
                    cout << "Decryption failure." << endl;
                    wrong_check++;
                }
                else {
                    cout << "Decryption success." << endl;
                }
            };
            uint64_t t_total_decoding_end = clock() - t_total_decoding_begin;

            t_total_keygen += t_keygen_end;
            t_total_encoding += t_encoding_end;
            t_total_decoding += t_total_decoding_end;
            t_total_ISD_one_attempt_end += t_ISD_one_attempt_end;

            cout << "iter#" << it+1 << "(stop): keygen clocks=" << t_keygen_end
                 <<  ", enc. clocks=" << t_encoding_end
                  << ", dec. clocks=" << t_total_decoding_end
                  << "(ISD one iter. clocks=" << t_ISD_one_attempt_end
                  << "), DFR=" << DFR << "/" << (it+1)
                  << ", wrong_checks=" << wrong_check
                  << endl;
            if(S_inv)
            {
                delete S_inv;
                S_inv = nullptr;
            };

            if (semi_pub_key)
            {
                delete semi_pub_key;
                semi_pub_key = nullptr;
            };

            if(pub_key)
            {
                delete pub_key;
                pub_key = nullptr;
            };
        };

        cout << "Average block decoding clocks=" << total_block_decoding_clocks/total_block_decoding_iters << ", average isd decoding clocks="  << total_isd_decoding_clocks/total_isd_decoding_iters << endl;

        delete id_matrix;

        RMCode::deinit_factor_classes_for_RPA();
        RMCode::deinit_M_matrices_for_Hadamar();
    }
}

void DCode::test_EncryptionDecryption_MatrixDecoder(const char* f_params, uint32_t N)
{
    std::ifstream dcode_file(f_params);
    if (dcode_file.is_open())
    {
        DCode dcode;
        uint8_t r, m;
        uint16_t weight;
        dcode_file >> weight;
        uint16_t K;
        dcode_file >> K;
        uint8_t iter = 0;
        uint64_t total_block_decoding_clocks = 0;
        uint64_t total_block_decoding_iters = 0;

        int r1, m1, r2, m2;
        while (dcode_file >> r1 >> m1 >> r2 >> m2)
        {
            dcode.add_tensor_product(r1, m1, r2, m2);
            if (iter == 0)
            {
                r = r2;
                m = m2;
            };

            iter++;
        }

        r2 = r;
        m2 = m;

        uint16_t num_blocks = 1 << m1;
        RMCode RM1(r1, m1);
        RMCode RM2(r2, m2);
        SBitsArray* code_vectors_of_small_code = new SBitsArray[1 << m2];
        for (uint16_t i = 0; i < (1 << m2); i++)
        {
            code_vectors_of_small_code[i].allocateBufferForBits(1 << m1);
        };

        dcode.print_params();
        cout << "Generator matrix" << endl;
        dcode.gen_matrix->print_nonzero_rows();

        RMCode RM(r, m, false);
        SBitsArray corrupted_code_vector((1 << m));
        SBitsArray info_vector(dcode.get_dimension());
        SBitsArray info_vector_dec(dcode.get_dimension());
        SBitsArray full_code_vector((1 << m) * (1 << dcode.rm_m1));
        SBitsArray full_error_vector((1 << m) * (1 << dcode.rm_m1));
        SBitsArray full_corrupted_vector((1 << m) * (1 << dcode.rm_m1));
        SBitsArray full_decoded_vector((1 << m) * (1 << dcode.rm_m1));
        uint32_t DFR = 0;
        uint32_t wrong_check = 0;
        SBitsMatrix* S_inv = nullptr;
        SBitsMatrix S_copy(dcode.get_dimension(), dcode.get_dimension());
        SBitsMatrix* identity = SBitsMatrix::identity_matrix(dcode.get_dimension());

        cout << "r1=" << uint16_t(r1) << ", m1=" << uint16_t(m1) << ", r2=" << uint16_t(r2) << ", m2=" << uint16_t(m2) << ", weight=" << weight << endl;
        cout << "start making factor classes" << endl;
        RMCode::init_factor_classes_for_RPA(m1);
        RM1.ns_fc = RMCode::fc;
        RMCode::init_factor_classes_for_RPA(m2);
        RM2.ns_fc = RMCode::fc;
        cout << "stop making factor classes" << endl;

        cout << "start making Hadamar matrices" << endl;
        RMCode::init_M_matrices_for_Hadamar(m1);
        RM1.ns_MMatrices = RMCode::MMatrices;
        RMCode::init_M_matrices_for_Hadamar(m2);
        RM2.ns_MMatrices = RMCode::MMatrices;
        cout << "stop making Hadamar matrices" << endl;

        /*
        cout << "start making factor classes" << endl;
        RMCode::init_factor_classes_for_RPA(m);
        cout << "stop making factor classes" << endl;

        cout << "start making Hadamar matrices" << endl;
        RMCode::init_M_matrices_for_Hadamar(m);
        cout << "stop making Hadamar matrices" << endl;
        */

        cout << "start making information set for nonzero vector" << endl;
        SBitsMatrix* id_matrix = SBitsMatrix::identity_matrix(dcode.get_dimension());
        SBitsMatrix gen_matr_copy(dcode.gen_matrix->get_number_of_rows(), dcode.gen_matrix->get_number_of_cols());
        gen_matr_copy.fast_copy(*(dcode.gen_matrix));
        gen_matr_copy.pseudo_systematic_form();
        quint64* information_set_base = new quint64[dcode.get_dimension()];
        for (uint64_t is = 0; is < dcode.get_dimension(); is++)
            information_set_base[is] = gen_matr_copy[is].getIndexOfFirstNonzero();
        cout << "stop making information set for nonzero vector" << endl;

        uint64_t t_total_keygen = 0;
        uint64_t t_total_encoding = 0;
        uint64_t t_total_decoding = 0;
        uint64_t t_total_ISD_one_attempt_end = 0;
        uint64_t rnd_seed = 0;
        for (uint32_t it = 0; it < N; it++)
        {
            SBitsMatrix gen_matr_copy2(dcode.gen_matrix->get_number_of_rows(), dcode.gen_matrix->get_number_of_cols());
            cout << "iter#" << it + 1 << "(start):" << endl;
            // Key generation
            uint64_t t_keygen_begin = clock();
            Permutation perm((1 << dcode.rm_m1) * ((1 << dcode.rm_m2)), clock());
            cout << "Column permutation (aka matrix P):" << endl;
            for (int p = 0; p < (1 << dcode.rm_m1) * (1 << dcode.rm_m2); p++)
                cout << perm.get(p) << ' ';
            cout << endl;
            do {
                rnd_seed++;
                //cout << "1" << endl;
                SBitsMatrix S(dcode.get_dimension(), dcode.get_dimension());
                //cout << "2" << endl;
                S.random(it + time(NULL) + rnd_seed);
                S_copy.fast_copy(S);
                //cout << "3" << endl;
                if (S.low_triangular_fast2() == dcode.get_dimension())
                {
                    //cout << "3.1" << endl;
                    S.fast_copy(S_copy);
                    //cout << "3.2" << endl;
                    S.concat_columns(*identity);
                    //cout << "3.3" << endl;
                    S.pseudo_systematic_form();
                    //cout << "3.4" << endl;
                    S_inv = S.get_last_n_columns(dcode.get_dimension());
                    //cout << "3.5" << endl;
                    break;
                };
                //cout << "4" << endl;
            } while (1);
            cout << "matrix S:" << endl;
            S_copy.print_nonzero_rows();
            SBitsMatrix* semi_pub_key = nullptr;
            SBitsMatrix* pub_key = nullptr;
            semi_pub_key = dcode.gen_matrix->mul_left_matrix(S_copy);
            pub_key = semi_pub_key->permute_columns(perm);


            uint64_t t_keygen_end = clock() - t_keygen_begin;
            cout << "Public key (G_pub):" << endl;
            pub_key->print_nonzero_rows();

            // Encryption
            info_vector.randomBitArray();
            cout << "Message (m):" << endl;
            cout << info_vector << endl;

            uint64_t t_encoding_begin = clock();
            full_code_vector = pub_key->mul_left(info_vector);
            cout << "Data encoded using public matrix (m*G_pub):" << endl;
            cout << full_code_vector << endl;

            full_error_vector.randomBitArrayWixedWeight_fast(weight);
            cout << "Error vector (e):" << endl;
            cout << full_error_vector << endl;
            full_corrupted_vector = SBitsArray::xorWithSBitArray_fast(full_code_vector, full_error_vector);
            uint64_t t_encoding_end = clock() - t_encoding_begin;
            cout << "Ciphertext (c=mG_pub+e):" << endl;
            cout << full_corrupted_vector << endl;

            // Decryption
            uint64_t t_total_decoding_begin = clock();
            full_corrupted_vector.inverse_permute(perm);
            cout << "Start decryption:" << endl << "Apply inverse permutation (c*P^{-1}):" << endl;
            cout << full_corrupted_vector << endl;

            /// apply MatrixDecoder

            // decode each codeword separately and count number of decoding failures
            uint16_t local_counter = 0;

            uint64_t t2_clock = clock();
            uint64_t distance_to_decoded = 0;
            uint64_t t_block_decoding_begin = clock();

            ///////////////////////////////////////////////////////////
            // MatrixDecoder (begin)
            /* Декодирование по строкам */
            for (uint16_t j = 0; j < num_blocks; j++)
            {
                corrupted_code_vector = full_corrupted_vector.getSubvectorInPosition_fast(j * (1 << m2), (1 << m2));
                SBitsArray ret = RM2.ns_RPA_decoder_fast2(corrupted_code_vector, (m2 % 2 == 0) ? m2 / 2 : m2 / 2 + 1);
                for (uint16_t sv = 0; sv < (1 << m2); sv++)
                {
                    code_vectors_of_small_code[sv].setBitValue_fast(j, ret.getBit_fast(sv));
                };
                full_decoded_vector.setSubvectorInPosition_fast(j * (1 << m2), ret);
            };

            cout << "Noisy codeword decoded by rows:" << endl;
            cout << full_decoded_vector << endl;

            /* Декодирование по столбцам */
            for (uint16_t sv = 0; sv < (1 << m2); sv++)
            {
                code_vectors_of_small_code[sv] = RM1.ns_RPA_decoder_fast2(code_vectors_of_small_code[sv], (m1 % 2 == 0) ? m1 / 2 : m1 / 2 + 1);

                for (uint16_t bv = 0; bv < (1 << m1); bv++)
                    full_decoded_vector.setBitValue_fast(bv * (1 << m2) + sv, code_vectors_of_small_code[sv].getBit_fast(bv));
            };
            // MatrixDecoder (end)
            ///////////////////////////////////////////////////////////

            cout << "Noisy codeword decoded by columns:" << endl;
            cout << full_decoded_vector << endl;

            uint64_t t_block_decoding_end = clock() - t_block_decoding_begin;
            total_block_decoding_clocks += t_block_decoding_end;
            total_block_decoding_iters++;

            // First attemp

            uint64_t t_ISD_one_attempt_begin = clock();
            SBitsArray z_proj(full_decoded_vector.getProjection_fast(information_set_base, dcode.get_dimension()));
            gen_matr_copy2.fast_copy(*(semi_pub_key));
            gen_matr_copy2.concat_columns(*id_matrix);
            gen_matr_copy2.pseudo_systematic_form();
            SBitsMatrix* inverse_base = gen_matr_copy2.get_last_n_columns(dcode.get_dimension());
            info_vector_dec = inverse_base->mul_left(z_proj);
            delete inverse_base;

            // Checking
            SBitsArray c(dcode.get_length());
            c = semi_pub_key->mul_left(info_vector_dec);
            c.xorWithSBitArray_fast(full_corrupted_vector);
            uint64_t t_ISD_one_attempt_end = clock() - t_ISD_one_attempt_begin;

            if (c.getHammingWeight_fast() != weight)
            {
                cout << "Decryption failure." << endl;
                DFR++;
            }
            else {
                info_vector_dec.xorWithSBitArray_fast(info_vector);
                if (info_vector_dec.getHammingWeight_fast() != 0)
                {
                    cout << "Decryption failure." << endl;
                    wrong_check++;
                }
                else {
                    cout << "Decryption success." << endl;
                }
            };
            uint64_t t_total_decoding_end = clock() - t_total_decoding_begin;

            t_total_keygen += t_keygen_end;
            t_total_encoding += t_encoding_end;
            t_total_decoding += t_total_decoding_end;
            t_total_ISD_one_attempt_end += t_ISD_one_attempt_end;

            cout << "iter#" << it + 1 << "(stop): keygen clocks=" << t_keygen_end
                << ", enc. clocks=" << t_encoding_end
                << ", dec. clocks=" << t_total_decoding_end
                << ", DFR=" << DFR << "/" << (it + 1)
                << ", wrong_checks=" << wrong_check
                << endl;
            if (S_inv)
            {
                delete S_inv;
                S_inv = nullptr;
            };

            if (semi_pub_key)
            {
                delete semi_pub_key;
                semi_pub_key = nullptr;
            };

            if (pub_key)
            {
                delete pub_key;
                pub_key = nullptr;
            };
        };

        cout << "Average block decoding clocks=" << total_block_decoding_clocks / total_block_decoding_iters << endl;

        delete id_matrix;

        RMCode::deinit_factor_classes_for_RPA();
        RMCode::deinit_M_matrices_for_Hadamar();
    }
}

void DCode::test_MatrixDecoder(char* f_params, uint32_t N)
{
    std::ifstream dcode_file(f_params);
    DCode dcode;
    uint16_t weight;
    dcode_file >> weight;
    uint16_t K;
    dcode_file >> K;
    uint8_t iter = 0;

    int r1, m1, r2, m2, r, m;
    while (dcode_file >> r1 >> m1 >> r2 >> m2)
    {
        dcode.add_tensor_product(r1, m1, r2, m2);
        if(iter == 0)
        {
            r = r2;
            m = m2;
        };
        iter++;
    };

    r2 = r;
    m2 = m;

    uint16_t num_blocks = 1 << m1;
    RMCode RM1(r1,m1);
    RMCode RM2(r2,m2);
    SBitsArray info_vector(dcode.get_dimension());
    SBitsArray *code_vectors = new SBitsArray[num_blocks];
    SBitsArray *code_vectors_of_small_code = new SBitsArray[1<<m2];
    SBitsArray corrupted_code_vector((1<<m2));
    SBitsArray full_code_vector((1<<m2)*num_blocks);
    SBitsArray full_error_vector((1<<m2)*num_blocks);
    SBitsArray full_corrupted_vector((1<<m2)*num_blocks);
    SBitsArray full_decoded_vector((1<<m2)*num_blocks);
    uint32_t *counters = new uint32_t[num_blocks+1];
    uint32_t number_of_wrong_codewords = 0;
    //uint32_t number_of_right_decryptions = 0;

    // initialize counters for decoding failures
    for(uint16_t i = 0; i < num_blocks+1; i++)
    {
        counters[i] = 0;
    };

    for(uint16_t i = 0; i < num_blocks; i++)
    {
        code_vectors[i].allocateBufferForBits(1<<m2);
    };

    for(uint16_t i = 0; i < (1<<m2); i++)
    {
        code_vectors_of_small_code[i].allocateBufferForBits(1<<m1);
    };

    cout << "r1=" << uint16_t(r1) << ", m1=" << uint16_t(m1) << ", r2=" << uint16_t(r2) << ", m2=" << uint16_t(m2) << ", weight=" << weight << endl;
    cout << "start making factor classes" << endl;
    RMCode::init_factor_classes_for_RPA(m1);
    RM1.ns_fc = RMCode::fc;
    RMCode::init_factor_classes_for_RPA(m2);
    RM2.ns_fc = RMCode::fc;
    cout << "stop making factor classes" << endl;

    cout << "start making Hadamar matrices" << endl;
    RMCode::init_M_matrices_for_Hadamar(m1);
    RM1.ns_MMatrices = RMCode::MMatrices;
    RMCode::init_M_matrices_for_Hadamar(m2);
    RM2.ns_MMatrices = RMCode::MMatrices;
    cout << "stop making Hadamar matrices" << endl;

    //cout << "start making permutation" << endl;
    //Permutation perm((1<<m)*num_blocks, clock());
    //cout << "stop making permutation" << endl;

    cout << "start decoding at " << time(NULL) << " s." << endl;
    for(uint32_t it = 0; it < N; it++)
    {
        Permutation perm((1<<m2)*num_blocks, clock());

        uint64_t t1 = time(NULL);
        uint64_t t1_clock = clock();
        // concatenate blocks
        info_vector.randomBitArray();
        full_code_vector = dcode.encode(info_vector);
        for(uint16_t j = 0; j < num_blocks; j++)
        {
            code_vectors[j] = full_code_vector.getSubvectorInPosition_fast(j*(1<<m2), (1<<m2));
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

        uint64_t t2_clock = clock();
        uint64_t distance_to_decoded = 0;
        for(uint16_t j = 0; j < num_blocks; j++)
        {
            corrupted_code_vector = full_corrupted_vector.getSubvectorInPosition_fast(j*(1<<m2), (1<<m2));
            SBitsArray ret = RM2.ns_RPA_decoder_fast2(corrupted_code_vector, (m2%2 == 0) ? m2/2 : m2/2+1);
            for(uint16_t sv = 0; sv < (1<<m2); sv++)
            {
                code_vectors_of_small_code[sv].setBitValue_fast(j, ret.getBit_fast(sv));
            };
            full_decoded_vector.setSubvectorInPosition_fast(j*(1<<m2), ret);
        };

        /* Декодирование по столбцам */
        for(uint16_t sv = 0; sv < (1<<m2); sv++)
        {
            code_vectors_of_small_code[sv] = RM1.ns_RPA_decoder_fast2(code_vectors_of_small_code[sv], (m1%2 == 0) ? m1/2 : m1/2+1);

            for(uint16_t bv = 0; bv < (1<<m1); bv++)
                full_decoded_vector.setBitValue_fast(bv*(1<<m2)+sv, code_vectors_of_small_code[sv].getBit_fast(bv));
        };
        /**/

        for(uint16_t bv = 0; bv < (1<<m1); bv++)
        {
            SBitsArray ret = full_decoded_vector.getSubvectorInPosition_fast(bv*(1<<m2), (1<<m2));
            SBitsArray synd = RM2.r_syndrom_fast(ret);
            ret.xorWithSBitArray_fast(code_vectors[bv]);
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

void DCode::test_MatrixDecoder_with_common_information(char* f_params, uint32_t N)
{
    std::ifstream dcode_file(f_params);
    DCode dcode;
    uint16_t weight;
    dcode_file >> weight;
    uint16_t K;
    dcode_file >> K;
    uint8_t iter = 0;

    int r1, m1, r2, m2, r, m;
    while (dcode_file >> r1 >> m1 >> r2 >> m2)
    {
        dcode.add_tensor_product(r1, m1, r2, m2);
        if(iter == 0)
        {
            r = r2;
            m = m2;
        };
        iter++;
    };

    r2 = r;
    m2 = m;

    uint16_t num_blocks = 1 << m1;
    RMCode RM1(r1,m1);
    RMCode RM2(r2,m2);
    SBitsArray info_vector(dcode.get_dimension());
    SBitsArray *code_vectors = new SBitsArray[num_blocks];
    SBitsArray *code_vectors_of_small_code = new SBitsArray[1<<m2];
    SBitsArray corrupted_code_vector((1<<m2));
    SBitsArray full_code_vector((1<<m2)*num_blocks);
    SBitsArray full_error_vector((1<<m2)*num_blocks);
    SBitsArray full_corrupted_vector((1<<m2)*num_blocks);
    SBitsArray full_decoded_vector((1<<m2)*num_blocks);
    uint32_t *counters = new uint32_t[num_blocks+1];
    uint32_t number_of_wrong_codewords = 0;
    //uint32_t number_of_right_decryptions = 0;

    // initialize counters for decoding failures
    for(uint16_t i = 0; i < num_blocks+1; i++)
    {
        counters[i] = 0;
    };

    for(uint16_t i = 0; i < num_blocks; i++)
    {
        code_vectors[i].allocateBufferForBits(1<<m2);
    };

    for(uint16_t i = 0; i < (1<<m2); i++)
    {
        code_vectors_of_small_code[i].allocateBufferForBits(1<<m1);
    };

    cout << "r1=" << uint16_t(r1) << ", m1=" << uint16_t(m1) << ", r2=" << uint16_t(r2) << ", m2=" << uint16_t(m2) << ", weight=" << weight << endl;
    cout << "start making factor classes" << endl;
    RMCode::init_factor_classes_for_RPA(m1);
    RM1.ns_fc = RMCode::fc;
    RMCode::init_factor_classes_for_RPA(m2);
    RM2.ns_fc = RMCode::fc;
    cout << "stop making factor classes" << endl;

    cout << "start making Hadamar matrices" << endl;
    RMCode::init_M_matrices_for_Hadamar(m1);
    RM1.ns_MMatrices = RMCode::MMatrices;
    RMCode::init_M_matrices_for_Hadamar(m2);
    RM2.ns_MMatrices = RMCode::MMatrices;
    cout << "stop making Hadamar matrices" << endl;

    cout << "start making codes for RPA decoder" << endl;
    RM1.init_codes_for_RPA_decoder(r1, m1);
    RM2.init_codes_for_RPA_decoder(r2, m2);

    //cout << "start making permutation" << endl;
    //Permutation perm((1<<m)*num_blocks, clock());
    //cout << "stop making permutation" << endl;

    cout << "start decoding at " << time(NULL) << " s." << endl;
    for(uint32_t it = 0; it < N; it++)
    {
        Permutation perm((1<<m2)*num_blocks, clock());

        uint64_t t1 = time(NULL);
        uint64_t t1_clock = clock();
        // concatenate blocks
        info_vector.randomBitArray();
        full_code_vector = dcode.encode(info_vector);
        for(uint16_t j = 0; j < num_blocks; j++)
        {
            code_vectors[j] = full_code_vector.getSubvectorInPosition_fast(j*(1<<m2), (1<<m2));
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

        uint64_t t2_clock = clock();
        uint64_t distance_to_decoded = 0;
        for(uint16_t j = 0; j < num_blocks; j++)
        {
            corrupted_code_vector = full_corrupted_vector.getSubvectorInPosition_fast(j*(1<<m2), (1<<m2));
            SBitsArray ret = RM2.ns_RPA_decoder_fast2_with_common_data(corrupted_code_vector, (m2%2 == 0) ? m2/2 : m2/2+1, r2, r2);
            for(uint16_t sv = 0; sv < (1<<m2); sv++)
            {
                code_vectors_of_small_code[sv].setBitValue_fast(j, ret.getBit_fast(sv));
            };
            full_decoded_vector.setSubvectorInPosition_fast(j*(1<<m2), ret);
        };

        /* Декодирование по столбцам */
        for(uint16_t sv = 0; sv < (1<<m2); sv++)
        {
            code_vectors_of_small_code[sv] = RM1.ns_RPA_decoder_fast2_with_common_data(code_vectors_of_small_code[sv], (m1%2 == 0) ? m1/2 : m1/2+1, r1, r1);

            for(uint16_t bv = 0; bv < (1<<m1); bv++)
                full_decoded_vector.setBitValue_fast(bv*(1<<m2)+sv, code_vectors_of_small_code[sv].getBit_fast(bv));
        };
        /**/

        for(uint16_t bv = 0; bv < (1<<m1); bv++)
        {
            SBitsArray ret = full_decoded_vector.getSubvectorInPosition_fast(bv*(1<<m2), (1<<m2));
            SBitsArray synd = RM2.r_syndrom_fast(ret);
            ret.xorWithSBitArray_fast(code_vectors[bv]);
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

