#include "crc.h"
#include <iostream>
using namespace std;

uint64_t crc64(uint64_t crc, const unsigned char *s, uint64_t l) {
    uint64_t j;

    for (j = 0; j < l; j++) {
        uint8_t byte = s[j];
        crc = crc64_tab[(uint8_t)crc ^ byte] ^ (crc >> 8);
    }
    return crc;
}

uint32_t crc32(unsigned long   crc, const unsigned char * buf, unsigned int len )
{
    while (len--)
    {
        crc = (crc >> 8) ^ autodintable[(crc ^ *buf++) & 0xff];
    }

    return crc;
}
