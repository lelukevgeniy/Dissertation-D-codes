#ifndef SHASH_CRC16_H
#define SHASH_CRC16_H

#include "../Crypto/Hashes/shash.h"
#include "../utils/crc.h"

class SHash_CRC16 : public SHash
{
public:
    SHash_CRC16()
        :SHash(nullptr,0)
    {}

    int32_t get_hash_byte_length() override
    {
        return 2;
    }

    void calculate_hash(unsigned char * buffer, int32_t buffer_length, unsigned char* hash_output)
    {
        ((uint16_t *)hash_output)[0] = (uint16_t) crc32(0, buffer, buffer_length);
    }
};

#endif // SHASH_CRC16_H
