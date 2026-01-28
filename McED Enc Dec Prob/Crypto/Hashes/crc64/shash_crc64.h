#ifndef SHASH_CRC64_H
#define SHASH_CRC64_H

#include "../Crypto/Hashes/shash.h"
#include "../utils/crc.h"

class SHash_CRC64 : public SHash
{
public:
    SHash_CRC64()
        :SHash(nullptr,0)
    {}

    int32_t get_hash_byte_length() override
    {
        return 8;
    }

    void calculate_hash(unsigned char * buffer, int32_t buffer_length, unsigned char* hash_output) override
    {
        ((uint64_t*)hash_output)[0] = crc64(0, buffer, buffer_length);
    }
};

#endif // SHASH_CRC64_H
