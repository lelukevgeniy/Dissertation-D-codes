#ifndef SRANDOMGENERATOR_OFBGAMMA_H
#define SRANDOMGENERATOR_OFBGAMMA_H

#include "../Crypto/RandomGenerators/srandomgenerator.h"
#include "../Crypto/BlockCiphers/mycrypto.hpp"

template <typename CipherType>
class SRandomGenerator_OFBGamma : public SRandomGenerator
{
public:
    SRandomGenerator_OFBGamma(const CipherType &alg, unsigned char *seed_v, int seed_length_v)
        :SRandomGenerator(seed_v, seed_length_v)
    {
        if(alg.block_lenght > seed_length_v)
        {
            unsigned char* temp_seed = new unsigned char[alg.block_lenght];
            memset(temp_seed, 0, alg.block_lenght);
            memcpy(temp_seed, seed, seed_length);
            ByteBlock byte_block(temp_seed, alg.block_lenght);
            generator = new CFB_Mode<CipherType>(alg, byte_block);
            delete [] temp_seed;
        }else{
            ByteBlock byte_block(seed, alg.block_lenght);
            generator = new CFB_Mode<CipherType>(alg, byte_block);
        }
    }

    uint32_t generate32() override
    {
        return generator->next_gamma_value();
    }

    uint64_t generate64() override
    {
        uint64_t temp64;

        ((uint32_t*)(&temp64))[0] = generator->next_gamma_value();
        ((uint32_t*)(&temp64))[1] = generator->next_gamma_value();

        return temp64;
    }

private:
    CFB_Mode<CipherType> *generator;
};

#endif // SRANDOMGENERATOR_OFBGAMMA_H
