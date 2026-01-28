#ifndef PERMUTATION_H
#define PERMUTATION_H


#include <cstdint>
class Permutation
{
public:
    Permutation(uint16_t n, uint64_t seed);
    ~Permutation();
    uint16_t get(uint16_t idx) {return perm[idx];};
    uint16_t get_inverse(uint16_t idx) {return inverse_perm[idx];};
    uint16_t get_length() { return perm_len; };

private:
    uint16_t *perm;
    uint16_t *inverse_perm;
    uint16_t perm_len;
};

#endif // PERMUTATION_H
