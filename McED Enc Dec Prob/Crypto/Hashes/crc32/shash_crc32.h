#ifndef SHASH_CRC32_H
#define SHASH_CRC32_H

#include "../Crypto/Hashes/shash.h"
#include "../utils/crc.h"

class SHash_CRC32 : public SHash
{
public:
    SHash_CRC32()
        :SHash(nullptr,0)
    {}

    int32_t get_hash_byte_length() override
    {
        return 4;
    }

    void calculate_hash(unsigned char * buffer, int32_t buffer_length, unsigned char* hash_output) override
    {
        ((uint32_t*)hash_output)[0] = crc32(0, buffer, buffer_length);
    }
};

#endif // SHASH_CRC32_H
