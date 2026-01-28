#ifndef SHASH_STRIBOG512_H
#define SHASH_STRIBOG512_H

#include "../Crypto/Hashes/Stribog/types.h"
#include "../Crypto/Hashes/Stribog/stribog.h"
#include "../Crypto/Hashes/shash.h"

class SHash_Stribog512 : public SHash
{
public:
    SHash_Stribog512()
    {
    }

    int32_t get_hash_byte_length() override {return 512/8;}

    virtual void calculate_hash(unsigned char * buffer, int32_t buffer_length, unsigned char* hash_output) {
        init(&ctx, HASH512);
        stribog(&ctx, buffer, buffer_length);

        for (int i = 0; i < OUTPUT_SIZE_512; i++)
            hash_output[i] = ctx.h[i];
    }

private:
    struct stribog_ctx_t ctx;
};

#endif // SHASH_STRIBOG512_H
