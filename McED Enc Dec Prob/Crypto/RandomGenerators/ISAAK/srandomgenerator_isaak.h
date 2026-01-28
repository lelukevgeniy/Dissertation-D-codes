#ifndef SRANDOMGENERATOR_ISAAK_H
#define SRANDOMGENERATOR_ISAAK_H

#include "../Crypto/RandomGenerators/srandomgenerator.h"
#include "sec_rand.h"

class SRandomGenerator_ISAAK : public SRandomGenerator
{
public:
    SRandomGenerator_ISAAK(unsigned char * seed_v, int seed_length_v)
        :SRandomGenerator(seed_v, seed_length_v){

        long sum = 0;
        for(int i = 0; i < seed_length; i++)
            sum += seed[i];
        ctx = sec_rand_init(sum);
    }

    uint32_t generate32() override
    {
        return sec_rand(&ctx);
    }

    uint64_t generate64() override
    {
        uint64_t temp64;

        ((uint32_t*)(&temp64))[0] = sec_rand(&ctx);
        ((uint32_t*)(&temp64))[1] = sec_rand(&ctx);

        return temp64;
    }

private:
    randctx ctx;
};

#endif // SRANDOMGENERATOR_ISAAK_H
