#include "permutation.h"
#include "../Crypto/RandomGenerators/ISAAK/srandomgenerator_isaak.h"
#include <algorithm>
#include <vector>
#include <random>
#include <algorithm>

Permutation::Permutation(uint16_t n, uint64_t seed)
    :perm(nullptr),
      inverse_perm(nullptr),
      perm_len(0)
{
    SRandomGenerator_ISAAK isaak((unsigned char*) &(seed), sizeof(uint64_t));

    perm = new uint16_t[n];
    inverse_perm = new uint16_t[n];
    perm_len = n;

    /*
    for(uint16_t i = 0; i < n; i++)
    {
        perm[i] = i;
    };

    // generate permutation
    for(uint16_t i = 0; i < n; i++){
        uint16_t randIdx = isaak.generate32() % n;
        uint16_t t = perm[i];
        perm[i] = perm[randIdx];
        perm[randIdx] = t;
    };
    */

    std::srand (clock());
    std::vector<int> myvector;
    for(uint16_t i = 0; i < n; i++)
    {
        myvector.push_back(i);
    };
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle (myvector.begin(), myvector.end(), g);

    for(uint16_t i = 0; i < n; i++)
    {
        perm[i] = myvector[i];
    };

    // generate inverse permutation
    for(uint16_t i = 0; i < n; i++){
        inverse_perm[perm[i]] = i;
    };

    /*
    for(uint16_t i = 0; i < n; i++){
        std::cout << perm[i] << " ";
    };
    std::cout << std::endl;
    for(uint16_t i = 0; i < n; i++){
        std::cout << inverse_perm[i] << " ";
    };
    std::cout << std::endl;
    */
}

Permutation::~Permutation()
{
    if(perm) delete [] perm;
    if(inverse_perm) delete [] inverse_perm;
}
