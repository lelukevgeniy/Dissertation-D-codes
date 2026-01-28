#ifndef SHASH_CRC8_H
#define SHASH_CRC8_H

#include "../Crypto/Hashes/shash.h"
#include "../utils/crc.h"

class SHash_CRC8 : public SHash
{
public:
    SHash_CRC8()
        :SHash(nullptr,0)
    {}

    int32_t get_hash_byte_length() override
    {
        return 1;
    }

    void calculate_hash(unsigned char * buffer, int32_t buffer_length, unsigned char* hash_output)
    {
        ((unsigned char *)hash_output)[0] = (unsigned char) crc32(0, buffer, buffer_length);
    }
};

#endif // SHASH_CRC8_H
