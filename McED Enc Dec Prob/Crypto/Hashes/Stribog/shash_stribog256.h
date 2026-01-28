#ifndef SHASH_STRIBOG256_H
#define SHASH_STRIBOG256_H

#include "../Crypto/Hashes/Stribog/types.h"
#include "../Crypto/Hashes/Stribog/stribog.h"
#include "../Crypto/Hashes/shash.h"

class SHash_Stribog256 : public SHash
{
public:
    SHash_Stribog256()
    {
    }

    int32_t get_hash_byte_length() override {return 256/8;}

    virtual void calculate_hash(unsigned char * buffer, int32_t buffer_length, unsigned char* hash_output) {
        init(&ctx, HASH256);
        stribog(&ctx, buffer, buffer_length);

        for (int i = 0; i < OUTPUT_SIZE_256; i++)
            hash_output[i] = ctx.h[i];
    }

private:
    struct stribog_ctx_t ctx;
};

#endif // SHASH_STRIBOG256_H
