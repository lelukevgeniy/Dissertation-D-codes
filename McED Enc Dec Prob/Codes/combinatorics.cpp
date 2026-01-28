#include "../Codes/combinatorics.h"

int64_t comb_factorial(int64_t n)
{
    int64_t f = 1;

    for(int64_t i = 1; i < n; i++)
        f *= (i+1);

    return f;
}

int64_t comb_choose(int64_t n, int64_t m)
{
    int64_t k = n - m;
    if (m > k)
        m = k;
    if (!m)
        return 1;
    int64_t akk = k = n - m + 1;
    k++;
    for (int64_t i = 2; i <= m; i++, k++)
        akk = akk / i * k + akk % i * k / i;
    return akk;
}


uint64_t comb_rm_dimension(uint8_t r, uint8_t m)
{
    uint64_t dim = 0;
    for (uint8_t i = 0; i <= r; i++)
        dim += comb_choose(m, i);
    return dim;
}
