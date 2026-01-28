#ifndef SRANDOMGENERATOR_H
#define SRANDOMGENERATOR_H

#include <iostream>

class SRandomGenerator
{
public:
    SRandomGenerator(unsigned char* seed_v, int seed_length_v)
        :seed(nullptr),
          seed_length(0)
    {
        if(seed_v && seed_length_v)
        {
            seed = new unsigned char[seed_length_v];
            memcpy(seed, seed_v, seed_length_v);
            seed_length = seed_length_v;
        }else{
            seed_length = 10;
            seed = new unsigned char[seed_length];
            for(int i = 0; i < seed_length; i++)
                seed[i] = rand();
        }
    }

    ~SRandomGenerator()
    {
		//std::cout << "~SRandomGenerator()" << std::endl;
		if(seed)
        {
			delete [] seed;
            seed = nullptr;
            seed_length = 0;
		}
    }

    virtual uint32_t generate32() = 0;
    virtual uint64_t generate64() = 0;

protected:
    unsigned char *seed;
    int seed_length;
};

#endif // SRANDOMGENERATOR_H
