#ifndef SHASH_H
#define SHASH_H

#include <stdint.h>

class SHash
{
public:
    SHash(){}
    SHash(unsigned char* seed, int32_t seed_length){}
    virtual int32_t get_hash_byte_length() {return 0;}
    virtual void calculate_hash(unsigned char * buffer, int32_t buffer_length, unsigned char* hash_output) {}
};

typedef SHash SHashNone;

#endif // SHASH_H
