#include "SBitsArray.h"
#include <iostream>
#include <fstream>
#include <time.h>
//#include <QFile>
#include "../Crypto/RandomGenerators/ISAAK/srandomgenerator_isaak.h"
#include "../Crypto/BlockCiphers/Kuznyechik.hpp"
#include "../Crypto/RandomGenerators/OFBGamma/srandomgenerator_ofbgamma.h"
#include <chrono>

using namespace std;

#define SEED_LENGTH 10

SBitsArray::SBitsArray()
    :buffer_for_bits(nullptr),
      size_in_bits(0),
      size_in_quint64(0)
{
    srand(clock());
}

SBitsArray::SBitsArray(const SBitsArray& copy)
    :SBitsArray()
{
	//cout << "SBitsArray(SBitsArray& copy)" << endl;
	setFromByteArray(copy.getDataAsBytes(), copy.getSizeInBits());
}

SBitsArray::SBitsArray(qint64 size_in_bits)
    :SBitsArray()
{
//	cout << "SBitsArray(qint64 size_in_bits)" << endl;
    allocateBufferForBits(size_in_bits);
}

SBitsArray::SBitsArray(unsigned char *src_buffer, qint64 size_in_bits)
    :SBitsArray()
{
    setFromByteArray(src_buffer, size_in_bits);
}

void SBitsArray::setFromByteArray(unsigned char *src_buffer, qint64 size_in_bits_val)
{
    if( size_in_bits_val > 0 )
    {
        allocateBufferForBits(size_in_bits_val);
        memcpy(buffer_for_bits, src_buffer, size_in_bits_val / 8 + (size_in_bits_val % 8 == 0 ? 0 : 1));
    }else{
        //cout << "Wrong parameters for SBitArray constructor" << endl;
    }
}

void SBitsArray::randomBitArray(std::string seed)
{
    SRandomGenerator_ISAAK generator((unsigned char*)(seed.c_str()), seed.length());

    for(int i = 0; i < size_in_quint64; i++)
        buffer_for_bits[i] = generator.generate64();
}

void SBitsArray::randomBitArray(qint64 seed)
{
    SRandomGenerator_ISAAK generator((unsigned char*)(&seed), sizeof(qint64));

    for(int i = 0; i < size_in_quint64; i++)
        buffer_for_bits[i] = generator.generate64();
}

void SBitsArray::randomBitArray()
{
    //srand(time(NULL));
    srand(chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    char seed[SEED_LENGTH + 1] = {0};
	
	for (int i = 0; i < SEED_LENGTH; i++)
	{
		seed[i] = rand() % 256;
	}

    std::string seed_str(seed);

	randomBitArray(seed_str);
}

bool SBitsArray::randomBitArrayWixedWeight(quint64 weight)
{
    int32_t c = int32_t(weight), i = 0;
    uint64_t seed = clock();
    SRandomGenerator_ISAAK isaak((unsigned char*) &(seed), sizeof(uint64_t));

    if(weight > size_in_bits)
        return false;

    toZeroAllBits();

    if(weight == 0)
        return true;

    while(c > 0){
        i = isaak.generate32() % size_in_bits;
        if(!getBit(i))
        {
            setBit(i);
            c = c - 1;
        };
    };
    return true;
}

void SBitsArray::randomBitArrayWixedWeight_fast(quint64 weight)
{
    int32_t c = int32_t(weight), i = 0;
    uint64_t seed = clock();
    SRandomGenerator_ISAAK isaak((unsigned char*) &(seed), sizeof(uint64_t));
    toZeroAllBits();
    while(c > 0){
        i = isaak.generate32() % size_in_bits;
        if(!getBit_fast(i))
        {
            setBit_fast(i);
            c = c - 1;
        };
    };
}

bool SBitsArray::randomBitArrayWeightNoMore(quint64 weight)
{
    srand(time(NULL));
    int w = rand()%(weight + 1);

    return randomBitArrayWixedWeight(w);
}

uint32_t SBitsArray::addNoiseWithErrorProbability(float p, uint64_t seed)
{
    uint32_t inversion_counter = 0;
    srand(seed);
    for(int i = 0; i < size_in_bits; i++)
    {
        if((rand()%10000) < (p*10000))
        {
            invertBit(i);
            inversion_counter++;
        }
    }
    return inversion_counter;
}

uint32_t SBitsArray::pass_through_bsc(float p, uint64_t seed)
{
    return addNoiseWithErrorProbability(p, seed);
}

uint32_t SBitsArray::pass_through_bsc_secure_rand_isaak(float p, uint64_t seed)
{
    SRandomGenerator_ISAAK isaak((unsigned char*) &(seed), sizeof(uint64_t));

    uint32_t inversion_counter = 0;
    for(int i = 0; i < size_in_bits; i++)
    {
        if((isaak.generate32()%10000) < (p*10000))
        {
            invertBit_fast(i);
            inversion_counter++;
        }
    }
    return inversion_counter;
}

uint32_t SBitsArray::pass_through_bsc_secure_rand_OFB(float p, uint64_t seed)
{
    BYTE key_bytes[40] = {  '1','1','1','1','1','1','1','1','1','1',
                            '2','2','2','2','2','2','2','2','2','2',
                            '3','3','3','3','3','3','3','3','3','3',
                            '4','4','4','4','4','4','4','4','4','4'};
    ByteBlock key(key_bytes, 32);
    Kuznyechik cipher(key);
    SRandomGenerator_OFBGamma<Kuznyechik> ofb_gamma(cipher, (unsigned char*) &(seed), sizeof(uint64_t));

    uint32_t inversion_counter = 0;
    for(int i = 0; i < size_in_bits; i++)
    {
        if((ofb_gamma.generate32()%10000) < (p*10000))
        {
            invertBit(i);
            inversion_counter++;
        }
    }
    return inversion_counter;
}

bool SBitsArray::getBit(qint64 bit_index) const
{
    quint64 mask = ((quint64)1) << (IN_QINT64_INDEX(bit_index));

    if((buffer_for_bits == nullptr) || bit_index > (size_in_bits-1))
        return false;

    return (buffer_for_bits[QINT64_INDEX(bit_index)] & mask) == 0 ? false : true;
}

bool SBitsArray::getBit_fast(qint64 bit_index) const
{
    quint64 mask = ((quint64)1) << (IN_QINT64_INDEX(bit_index));

    return (buffer_for_bits[QINT64_INDEX(bit_index)] & mask) == 0 ? false : true;
}

bool SBitsArray::setBitValue(qint64 bit_index, bool value)
{
    quint64 mask = ((quint64)1) << (IN_QINT64_INDEX(bit_index));

    if((buffer_for_bits == nullptr) || (bit_index > (size_in_bits-1)))
        return false;

    if(value)
    {
        buffer_for_bits[QINT64_INDEX(bit_index)] = buffer_for_bits[QINT64_INDEX(bit_index)] | mask;
    }else{
        buffer_for_bits[QINT64_INDEX(bit_index)] = buffer_for_bits[QINT64_INDEX(bit_index)] & (~mask);
    }
    return true;
}

void SBitsArray::setBitValue_fast(qint64 bit_index, bool value)
{
    quint64 mask = ((quint64)1) << (IN_QINT64_INDEX(bit_index));

    if(value)
    {
        buffer_for_bits[QINT64_INDEX(bit_index)] = buffer_for_bits[QINT64_INDEX(bit_index)] | mask;
    }else{
        buffer_for_bits[QINT64_INDEX(bit_index)] = buffer_for_bits[QINT64_INDEX(bit_index)] & (~mask);
    }
}

bool SBitsArray::setBitsValueFromQint64(qint64 start_bit_index, qint16 number_of_low_bits, qint64 value)
{
    if((buffer_for_bits == nullptr) || (start_bit_index + number_of_low_bits) > size_in_bits)
        return false;

    for(int i = 0; i < number_of_low_bits; i++)
    {
        bool bool_value = (bool) ((value >> i) & 1);
        quint64 mask = ((quint64)1) << (IN_QINT64_INDEX(start_bit_index + i));
        if(bool_value)
        {
            buffer_for_bits[QINT64_INDEX(start_bit_index + i)] = buffer_for_bits[QINT64_INDEX(start_bit_index + i)] | mask;
        }else{
            buffer_for_bits[QINT64_INDEX(start_bit_index + i)] = buffer_for_bits[QINT64_INDEX(start_bit_index + i)] & (~mask);
        }
    }

    return true;
}

void SBitsArray::setBitsValueFromQint64_fast(qint64 value)
{
    buffer_for_bits[0] = value;
}

bool SBitsArray::setBit(const qint64 bit_index)
{
    return setBitValue(bit_index, true);
}

void SBitsArray::setBit_fast(const qint64 bit_index)
{
    setBitValue_fast(bit_index, true);
}

bool SBitsArray::unsetBit(const qint64 bit_index)
{
    return setBitValue(bit_index, false);
}

bool SBitsArray::xorBit(qint64 bit_index, bool value)
{
    if((buffer_for_bits == nullptr) || bit_index > (size_in_bits-1))
        return false;
    else{
        if(value)
        {
            quint64 mask = ((quint64)1) << (IN_QINT64_INDEX(bit_index));
            buffer_for_bits[QINT64_INDEX(bit_index)] = buffer_for_bits[QINT64_INDEX(bit_index)] ^ mask;
        }
        return true;
    }
}

void SBitsArray::xorBit_fast(qint64 bit_index, bool value)
{
    quint64 mask = ((quint64)1) << (IN_QINT64_INDEX(bit_index));
    buffer_for_bits[QINT64_INDEX(bit_index)] = buffer_for_bits[QINT64_INDEX(bit_index)] ^ mask;
}

bool SBitsArray::orBit(qint64 bit_index, bool value)
{
    if((buffer_for_bits == nullptr) || bit_index > (size_in_bits-1))
        return false;
    else{
        if(value)
        {
            setBit(bit_index);
        }
        return true;
    }
}

bool SBitsArray::andBit(qint64 bit_index, bool value)
{
    if((buffer_for_bits == nullptr) || bit_index > (size_in_bits-1))
        return false;
    else{
        if(!value)
        {
            unsetBit(bit_index);
        }
        return true;
    }
}


bool SBitsArray::invertBit(qint64 bit_index)
{
    return xorBit(bit_index, true);
}

void SBitsArray::invertBit_fast(qint64 bit_index)
{
    xorBit_fast(bit_index, true);
}

qint64 SBitsArray::getIndexOfFirstNonzero()
{
    qint64 index = -1;

    for(quint64 i = 0; i < size_in_bits; i++)
        if(getBit_fast(i) == true)
        {
            index = i;
            break;
        };

    return index;
}

void SBitsArray::clearBufferForBits()
{
    if(buffer_for_bits != nullptr)
    {
        delete [] buffer_for_bits;
        buffer_for_bits = nullptr;
        size_in_bits = 0;
        size_in_quint64 = 0;
    }
}

bool SBitsArray::allocateBufferForBits(qint64 number_of_bits)
{
    bool result;

    if (number_of_bits <= 0)
        return false;

    clearBufferForBits();

    size_in_quint64 = IN_QINT64_INDEX(number_of_bits) == 0 ? QINT64_INDEX(number_of_bits) : QINT64_INDEX(number_of_bits) + 1;
    //cout << size_in_quint64 << endl;
    size_in_bits = number_of_bits;

    buffer_for_bits = new quint64[static_cast<unsigned int>(size_in_quint64)];

    if(buffer_for_bits != nullptr)
    {
        result = true;
        // zeroing byffer
        memset(buffer_for_bits, 0, static_cast<unsigned int>(size_in_quint64) * sizeof(quint64));
    }else{
        result = false;
        size_in_bits = 0;
        size_in_quint64 = 0;
    }
    return result;
}

bool SBitsArray::allocateBufferForPotentialBits(qint64 number_of_bits, qint64 potentianl_number_of_bits)
{
    bool result;

    if (number_of_bits <= 0)
        return false;

    clearBufferForBits();

    size_in_quint64 = IN_QINT64_INDEX(potentianl_number_of_bits) == 0 ? QINT64_INDEX(potentianl_number_of_bits) : QINT64_INDEX(potentianl_number_of_bits) + 1;
    //cout << size_in_quint64 << endl;
    size_in_bits = number_of_bits;

    buffer_for_bits = new quint64[static_cast<unsigned int>(size_in_quint64)];

    if(buffer_for_bits != nullptr)
    {
        result = true;
        // zeroing byffer
        memset(buffer_for_bits, 0, static_cast<unsigned int>(size_in_quint64) * sizeof(quint64));
    }else{
        result = false;
        size_in_bits = 0;
        size_in_quint64 = 0;
    }
    return result;
}

ostream& operator<<(ostream& os, const SBitsArray& bit_array)
{
    for(int i = 0; i < bit_array.size_in_bits; i++)
        os << (bit_array.getBit(i) == true ? '1' : '0');
    return os;
}

void SBitsArray::concat(SBitsArray& bit_array)
{
    quint64 *new_buffer = new quint64[SIZE_IN_QINT64(getSizeInBits() + bit_array.getSizeInBits())];

    //quint64 *new_buffer = new quint64[getSizeInQuint64() + bit_array.getSizeInQuint64()];
    // copy bytes
    memcpy(new_buffer, buffer_for_bits, sizeof(quint64)*getSizeInQuint64());
    // free old data
    delete [] buffer_for_bits;
    // assign new buffer

    quint64 old_size_in_bits = getSizeInBits();

    buffer_for_bits = new_buffer;
    size_in_quint64 = SIZE_IN_QINT64(getSizeInBits() + bit_array.getSizeInBits());
    size_in_bits = getSizeInBits() + bit_array.getSizeInBits();

    for(int i = old_size_in_bits; i < old_size_in_bits + bit_array.getSizeInBits(); i++)
    {
        setBitValue(i, bit_array.getBit(i-old_size_in_bits));
    }
}

void SBitsArray::concat_with_bit(bool bit)
{
    if( ((size_in_bits % (8*sizeof(quint64))) != 0) || (size_in_bits == 0))
    {
        quint64 *new_buffer = new quint64[SIZE_IN_QINT64(getSizeInBits() + 1)];
        memcpy(new_buffer, buffer_for_bits, sizeof(quint64)*getSizeInQuint64());
        // free old data
        delete [] buffer_for_bits;
        // assign new buffer

        quint64 old_size_in_bits = getSizeInBits();

        buffer_for_bits = new_buffer;
        size_in_quint64 = SIZE_IN_QINT64(getSizeInBits() + 1);
        size_in_bits = getSizeInBits() + 1;
        setBitValue(old_size_in_bits, bit);
    }else{
        setBitValue(size_in_bits, bit);
        size_in_bits++;
    }
}

void SBitsArray::concat_with_bit_fast(bool bit)
{
    if( ((size_in_bits % (8*sizeof(quint64))) != 0) || (size_in_bits == 0))
    {
        quint64 *new_buffer = new quint64[SIZE_IN_QINT64(getSizeInBits() + 1)];
        memcpy(new_buffer, buffer_for_bits, sizeof(quint64)*getSizeInQuint64());
        // free old data
        delete [] buffer_for_bits;
        // assign new buffer

        quint64 old_size_in_bits = getSizeInBits();

        buffer_for_bits = new_buffer;
        size_in_quint64 = SIZE_IN_QINT64(getSizeInBits() + 1);
        size_in_bits = getSizeInBits() + 1;
        setBitValue(old_size_in_bits, bit);
    }else{
        setBitValue_fast(size_in_bits, bit);
        size_in_bits++;
    }
}

unsigned char *SBitsArray::getCopyAsBytes()
{
    unsigned char *copy = new unsigned char[static_cast<unsigned int>(sizeof(quint64)*size_in_quint64)];
    memcpy(copy, reinterpret_cast<unsigned char*>(buffer_for_bits), static_cast<unsigned int>(sizeof(quint64)*size_in_quint64));
    return copy;
}

unsigned char *SBitsArray::getDataAsBytes() const
{
    return reinterpret_cast<unsigned char*>(buffer_for_bits);
}

void SBitsArray::setSubvectorInPosition_fast(quint64 pos, SBitsArray& subvector)
{
    if( (static_cast<quint64>(pos + subvector.getSizeInBits())) > size_in_bits)
    {
        cout << "Place of subvector is out of range"  << endl;
        return;
    }

    for(int i = pos; i < pos + subvector.getSizeInBits(); i++)
    {
        setBitValue_fast(i, subvector.getBit_fast(i-pos));
    }
}

bool SBitsArray::setSubvectorInPosition(quint64 pos, SBitsArray& subvector)
{
    if( (static_cast<quint64>(pos + subvector.getSizeInBits())) > size_in_bits)
    {
        cout << "Place of subvector is out of range"  << endl;
        return false;
    }

    for(int i = pos; i < pos + subvector.getSizeInBits(); i++)
    {
        setBitValue(i, subvector.getBit(i-pos));
    }
    return true;
}

SBitsArray* SBitsArray::getSubvectorInPosition(quint64 pos, quint64 len)
{
    if( ((static_cast<quint64>(pos + len)) > size_in_bits) || (len == 0))
    {
        return nullptr;
    }

    SBitsArray *subvector = new SBitsArray(len);

    for(int i = 0; i < len; i++)
    {
        subvector->setBitValue(i, getBit(pos + i));
    }

    return subvector;
}

SBitsArray SBitsArray::getSubvectorInPosition_fast(quint64 pos, quint64 len) const
{
    SBitsArray subvector(len);
    for(int i = 0; i < len; i++)
    {
        subvector.setBitValue_fast(i, getBit_fast(pos + i));
    }

    return subvector;
}

bool SBitsArray::deleteSubVectorInPosition(quint64 pos, quint64 len)
{
    bool result = true;
    SBitsArray *return_value = cutSubvectorInPosition(pos, len);

    if(return_value)
    {
        delete return_value;
    }else{
        result = false;
    }

    return result;
}

SBitsArray* SBitsArray::cutSubvectorInPosition(quint64 pos, quint64 len)
{
    SBitsArray *prefix, *postfix, *return_value;

    prefix = getSubvectorInPosition(0, pos);
    return_value = getSubvectorInPosition(pos, len);
    postfix = getSubvectorInPosition(pos + len, size_in_bits - (pos + len));

    if(return_value == nullptr)
    {
        if(prefix != nullptr)
            delete prefix;

        if(postfix != nullptr)
            delete postfix;

        return nullptr;
    }

    if((prefix == nullptr) && ((postfix != nullptr)))
    {
        prefix = postfix;
    }else{
        if((prefix != nullptr) && ((postfix == nullptr)))
        {
            prefix = prefix;
        }else{
            if((prefix != nullptr) && ((postfix != nullptr)))
            {
                prefix->concat(*postfix);
                delete postfix;
            }
        }
    }

    clearBufferForBits();
    setFromByteArray(prefix->getDataAsBytes(), prefix->getSizeInBits());
    delete prefix;

    return return_value;
}

qint64 SBitsArray::getHammingDistance(SBitsArray& compared_bit_array)
{
    int distance = -1;
    if(size_in_bits == compared_bit_array.getSizeInBits())
    {
        distance++;
        for(int i = 0; i < size_in_bits; i++)
            if(getBit(i) != compared_bit_array.getBit(i))
                distance++;
    }

    return distance;
}

qint64 SBitsArray::getHammingWeight()
{
    int weight = 0;
    for(int i = 0; i < size_in_bits; i++)
        if(getBit(i)) weight++;

    return weight;
}

qint64 SBitsArray::getHammingWeight_fast()
{
    int weight = 0;

    if(size_in_bits%(8*sizeof(uint64_t)) == 0)
    {
        for(uint64_t i = 0; i < size_in_quint64; i++)
            weight += hamming_weight_of_uint64(buffer_for_bits[i]);
    }else{

        for(uint64_t i = 0; i < size_in_quint64 - 1; i++)
            weight += hamming_weight_of_uint64(buffer_for_bits[i]);

        weight += hamming_weight_of_uint64(buffer_for_bits[size_in_quint64 - 1] & ((((uint64_t)1) << (size_in_bits%(8*sizeof(uint64_t)))) - 1));
    }

    return weight;
}

quint64 SBitsArray::getLowestBitsAsInteger(int number_of_lowest_bits)
{
    quint64 value = 0;

    for(int i = 0; i < min(number_of_lowest_bits%64, (int)size_in_bits); i++)
    {
        if(getBit(i))
            value += (1<<i);
    }

    return value;
}

quint64 SBitsArray::getLowestBitsAsInteger_fast(int number_of_lowest_bits)
{
    uint64_t min_val = min(number_of_lowest_bits%64, (int)size_in_bits);
    return buffer_for_bits[0] & ((uint64_t(1) << min_val)-1);
}

quint64 SBitsArray::getNBitsFromPositionAsInteger(int position, int number_of_bits) const
{
    quint64 value = 0;

    if(position >= size_in_bits)
    {
        cout << "Wrong position" << endl;
        return value;
    }

    for(int i = position; i < min(position + number_of_bits, (int)size_in_bits); i++)
    {
        if(getBit(i))
            value += (1<<(i-position));
    }

    return value;
}

bool SBitsArray::xorWithSBitArray(const SBitsArray& xored_bit_array)
{
    bool result = false;

    if(size_in_bits == xored_bit_array.getSizeInBits())
    {
        for(int i = 0; i < size_in_bits; i++)
            xorBit(i, xored_bit_array.getBit(i));
        result = true;
    }

    return result;
}

void SBitsArray::xorWithSBitArray_fast(const SBitsArray& xored_bit_array)
{
    for(int i = 0; i < size_in_quint64; i++)
    {
        buffer_for_bits[i] ^= xored_bit_array.buffer_for_bits[i];
    }
}

void SBitsArray::orWithSBitArray_fast(const SBitsArray& xored_bit_array)
{
    for(int i = 0; i < size_in_quint64; i++)
    {
        buffer_for_bits[i] |= xored_bit_array.buffer_for_bits[i];
    }
}

void SBitsArray::andWithSBitArray_fast(const SBitsArray& xored_bit_array)
{
    for(int i = 0; i < size_in_quint64; i++)
    {
        buffer_for_bits[i] &= xored_bit_array.buffer_for_bits[i];
    }
}

SBitsArray SBitsArray::xorWithSBitArray(const SBitsArray& bit_array1, const SBitsArray& bit_array2)
{
    SBitsArray result(bit_array1);
    result.xorWithSBitArray(bit_array2);
    return result;
}

SBitsArray SBitsArray::xorWithSBitArray_fast(const SBitsArray& bit_array1, const SBitsArray& bit_array2)
{
    SBitsArray result(bit_array1);
    result.xorWithSBitArray_fast(bit_array2);
    return result;
}

SBitsArray SBitsArray::andWithSBitArray(const SBitsArray& bit_array1, const SBitsArray& bit_array2)
{
    SBitsArray result(bit_array1);
    result.andWithSBitArray(bit_array2);
    return result;
}

SBitsArray SBitsArray::andWithSBitArray_fast(const SBitsArray& bit_array1, const SBitsArray& bit_array2)
{
    SBitsArray result(bit_array1);
    result.andWithSBitArray_fast(bit_array2);
    return result;
}

SBitsArray SBitsArray::orWithSBitArray(const SBitsArray& bit_array1, const SBitsArray& bit_array2)
{
    SBitsArray result(bit_array1);
    result.orWithSBitArray(bit_array2);
    return result;
}

SBitsArray SBitsArray::orWithSBitArray_fast(const SBitsArray& bit_array1, const SBitsArray& bit_array2)
{
    SBitsArray result(bit_array1);
    result.orWithSBitArray_fast(bit_array2);
    return result;
}

bool SBitsArray::xorWithSBitArrayFromPosition(const SBitsArray& xored_bit_array, int position)
{
    bool result = false;

    if(position + xored_bit_array.getSizeInBits() <= size_in_bits)
    {
        for(int i = position; i < position + xored_bit_array.getSizeInBits(); i++)
            xorBit(i, xored_bit_array.getBit(i-position));
        result = true;
    }else{
        cout << "XOR operation is out of range" << endl;
    }

    return result;
}

bool SBitsArray::orWithSBitArray(const SBitsArray& ored_bit_array)
{
    bool result = false;

    if(size_in_bits == ored_bit_array.getSizeInBits())
    {
        for(int i = 0; i < size_in_bits; i++)
            orBit(i, ored_bit_array.getBit(i));
        result = true;
    }

    return result;
}

bool SBitsArray::andWithSBitArray(const SBitsArray& anded_bit_array)
{
    bool result = false;

    if(size_in_bits == anded_bit_array.getSizeInBits())
    {
        for(int i = 0; i < size_in_bits; i++)
            andBit(i, anded_bit_array.getBit(i));
        result = true;
    }

    return result;
}

void SBitsArray::invertAllBits()
{
    for(int i = 0; i < size_in_bits; i++)
        invertBit(i);
}

void SBitsArray::toZeroAllBits()
{
    for(int i = 0; i < size_in_quint64; i++)
        buffer_for_bits[i] = 0;
}

void SBitsArray::toOneAllBits()
{
    for(int i = 0; i < size_in_quint64; i++)
    {
        buffer_for_bits[i] = 0;
        buffer_for_bits[i] = ~(buffer_for_bits[i]);
    }
}

bool SBitsArray::getHammingWeightMod2()
{
    return getHammingWeight()%2;
}

bool SBitsArray::getHammingWeightMod2_fast()
{
    return getHammingWeight_fast()%2;
}

bool SBitsArray::insertSubvectorInPosition(quint64 pos, SBitsArray& subvector)
{
    if( (static_cast<quint64>(pos)) > size_in_bits)
    {
        cout << "Place of subvector is out of range"  << endl;
        return false;
    }

    if(pos == 0)
    {
        SBitsArray *result = subvector.getSubvectorInPosition(0, subvector.getSizeInBits());
        result->concat(*this);
        clearBufferForBits();
        setFromByteArray(result->getDataAsBytes(), result->getSizeInBits());
        delete result;
     }else{
        if(pos == getSizeInBits())
        {
            concat(subvector);
         }else{
            SBitsArray *prefix = getSubvectorInPosition(0, pos);
            SBitsArray *postfix = getSubvectorInPosition(pos, size_in_bits - pos);
            prefix->concat(subvector);
            prefix->concat(*postfix);
            clearBufferForBits();
            setFromByteArray(prefix->getDataAsBytes(), prefix->getSizeInBits());
            delete prefix;
            delete postfix;
        }
    }
    return true;
}

void SBitsArray::cyclicShiftLeft(int shift)
{
    int normalized_shift = shift % size_in_bits;
    if(normalized_shift > 0)
    {
        SBitsArray *left_tail = cutSubvectorInPosition(0, normalized_shift);
        concat(*left_tail);
        delete left_tail;
    }
}

void SBitsArray::cyclicShiftRight(int shift)
{
    int normalized_shift = shift % size_in_bits;
    if(normalized_shift > 0)
    {
        SBitsArray *right_tail = cutSubvectorInPosition(size_in_bits - normalized_shift, normalized_shift);
        right_tail->concat(*this);
        clearBufferForBits();
        setFromByteArray(right_tail->getDataAsBytes(), right_tail->getSizeInBits());
        delete right_tail;
    }
}

void SBitsArray::permute(Permutation &perm)
{
    SBitsArray clone(size_in_bits);
    clone.fast_copy_without_check(*this);
    for(uint16_t i = 0; i < size_in_bits; i++)
    {
        setBitValue_fast(i, clone.getBit_fast(perm.get(i)));
    }
}

void SBitsArray::inverse_permute(Permutation &perm)
{
    SBitsArray clone(size_in_bits);
    clone.fast_copy_without_check(*this);
    for(uint16_t i = 0; i < size_in_bits; i++)
    {
        setBitValue_fast(i, clone.getBit_fast(perm.get_inverse(i)));
    }
}

bool operator==(const SBitsArray& bit_array_left, const SBitsArray& bit_array_right)
{
    bool thesame = true;

    if(bit_array_left.size_in_bits != bit_array_right.size_in_bits)
    {
        return false;
    }

    for(int i = 0; i < bit_array_left.size_in_quint64 - 1; i++)
        if(bit_array_left.buffer_for_bits[i] ^ bit_array_right.buffer_for_bits[i])
        {
            thesame = false;
        }

    if(thesame)
    {
        if( (bit_array_left.buffer_for_bits[bit_array_left.size_in_quint64 - 1] ^ bit_array_right.buffer_for_bits[bit_array_left.size_in_quint64 - 1])
                & ((((quint64)1) << (bit_array_left.size_in_bits%64)) - 1) )
        {
            thesame = false;
        }
    }
    return thesame;
}

SBitsArray& SBitsArray::operator=(const SBitsArray& bit_array)
{
    if(this == &bit_array)
        return *this;

    if(size_in_bits != bit_array.size_in_bits)
    {
        cout << "Dimensions do not match" << endl;
        return *this;
    }

    if(size_in_bits > 0)
    {
        clearBufferForBits();
        setFromByteArray(bit_array.getDataAsBytes(), bit_array.size_in_bits);
    }
    return *this;
}

SBitsArray& SBitsArray::fast_copy(const SBitsArray& copy)
{
    if(size_in_bits != copy.size_in_bits)
    {
        cout << "Dimensions do not match" << endl;
        return *this;
    }

    memcpy((uint8_t*)buffer_for_bits, (uint8_t*)copy.buffer_for_bits, sizeof(uint64_t)*size_in_quint64);

    return *this;
}

SBitsArray& SBitsArray::fast_copy_without_check(const SBitsArray& copy)
{
    memcpy((uint8_t*)buffer_for_bits, (uint8_t*)copy.buffer_for_bits, sizeof(uint64_t)*size_in_quint64);
    return *this;
}

SBitsArray SBitsArray::getProjection(quint64* positions, quint64 number_of_positions)
{
    SBitsArray result(number_of_positions);
    uint64_t counter = 0;
    for(counter = 0; counter < number_of_positions; counter++)
    {
        result.setBitValue(counter, getBit(positions[counter]));
    }
    return result;
}

SBitsArray SBitsArray::getProjection_fast(quint64* positions, quint64 number_of_positions)
{
    SBitsArray result(number_of_positions);
    uint64_t counter = 0;
    for(counter = 0; counter < number_of_positions; counter++)
    {
        result.setBitValue_fast(counter, getBit_fast(positions[counter]));
    }
    return result;
}

qint64 SBitsArray::getNumberOfSignificantBytes() const
{
    return  (size_in_bits % 8) == 0 ? size_in_bits / 8 : size_in_bits / 8 + 1;
}

qint64 SBitsArray::getNumberOfSignificantBytes(quint64 size_in_bits_c_v)
{
    return  (size_in_bits_c_v % 8) == 0 ? size_in_bits_c_v / 8 : size_in_bits_c_v / 8 + 1;
}

bool SBitsArray::save_to_file(std::string file_name)
{
    bool result = false;

    if(!file_name.empty())
    {
        ofstream f(file_name, std::ios::out | std::ios::binary | std::ios::app);
        if (f.is_open()){
            f.write((const char*)(buffer_for_bits), (std::streamsize)getNumberOfSignificantBytes());
            f.close();
            result = true;
        }
    }

    return result;
}

uint8_t SBitsArray::hamming_weight_of_uint64(uint64_t val)
{
    uint8_t count = 0;
    while (val) {
        val &= (val - 1);
        count++;
    };
    return count;
}
