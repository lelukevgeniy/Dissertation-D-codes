#ifndef SBITSARRAY_H
#define SBITSARRAY_H

#include "d_types.h"
#include <ostream>
#include <iostream>
#include "permutation.h"
using namespace std;

#define QINT64_INDEX(a) ((a) / (sizeof(quint64)*8))
#define IN_QINT64_INDEX(a) ((a) % (sizeof(quint64)*8))
#define SIZE_IN_QINT64(a) (QINT64_INDEX(a)+ (IN_QINT64_INDEX(a) == 0 ? 0 : 1))

class SBitsArray{
public:
    SBitsArray();
    SBitsArray(const SBitsArray& copy);
    SBitsArray(qint64 size_in_bits);
    SBitsArray(unsigned char *src_buffer, qint64 size_in_bits);
	~SBitsArray() {
		clearBufferForBits();
	}

    SBitsArray& fast_copy(const SBitsArray& copy);
    SBitsArray& fast_copy_without_check(const SBitsArray& copy);

    void setFromByteArray(unsigned char *src_buffer, qint64 size_in_bits);
    void randomBitArray(std::string seed);
    void randomBitArray(qint64 seed);
	void randomBitArray();
    bool randomBitArrayWixedWeight(quint64 weight);
    void randomBitArrayWixedWeight_fast(quint64 weight);
    bool randomBitArrayWeightNoMore(quint64 weight);
    uint32_t addNoiseWithErrorProbability(float p, uint64_t seed);
    uint32_t pass_through_bsc(float p, uint64_t seed);
    uint32_t pass_through_bsc_secure_rand_isaak(float p, uint64_t seed);
    uint32_t pass_through_bsc_secure_rand_OFB(float p, uint64_t seed);
    void clearBufferForBits();
    bool allocateBufferForBits(qint64 number_of_bits);
    bool allocateBufferForPotentialBits(qint64 number_of_bits, qint64 potentianl_number_of_bits);

    // single bit operations
    bool getBit(qint64 bit_index) const;
    bool getBit_fast(qint64 bit_index) const;
    bool setBitValue(qint64 bit_index, bool value);
    void setBitValue_fast(qint64 bit_index, bool value);
    bool setBitsValueFromQint64(qint64 start_bit_index, qint16 number_of_low_bits, qint64 value);
    void setBitsValueFromQint64_fast(qint64 value);
    bool setBit(const qint64 bit_index);
    void setBit_fast(const qint64 bit_index);
    bool unsetBit(const qint64 bit_index);
    bool xorBit(qint64 bit_index, bool value);
    void xorBit_fast(qint64 bit_index, bool value);
    bool orBit(qint64 bit_index, bool value);
    bool andBit(qint64 bit_index, bool value);
    bool invertBit(qint64 bit_index);
    void invertBit_fast(qint64 bit_index);
    qint64 getIndexOfFirstNonzero();

    // operations with arrays
    void concat(SBitsArray& bit_array);
    void concat_with_bit(bool bit);
    void concat_with_bit_fast(bool bit);
    bool setSubvectorInPosition(quint64 pos, SBitsArray& subvector);
    void setSubvectorInPosition_fast(quint64 pos, SBitsArray& subvector);
    bool deleteSubVectorInPosition(quint64 pos, quint64 len);
    bool insertSubvectorInPosition(quint64 pos, SBitsArray& subvector);
    SBitsArray* getSubvectorInPosition(quint64 pos, quint64 len);
    SBitsArray getSubvectorInPosition_fast(quint64 pos, quint64 len) const;
    SBitsArray* cutSubvectorInPosition(quint64 pos, quint64 len);
    SBitsArray getProjection(quint64* positions, quint64 number_of_positions);
    SBitsArray getProjection_fast(quint64* positions, quint64 number_of_positions);
    qint64 getHammingDistance(SBitsArray& compared_bit_array);
    qint64 getHammingWeight();
    qint64 getHammingWeight_fast();
    quint64 getLowestBitsAsInteger(int number_of_lowest_bits);
    quint64 getLowestBitsAsInteger_fast(int number_of_lowest_bits);
    quint64 getNBitsFromPositionAsInteger(int position, int number_of_bits) const;
    bool xorWithSBitArray(const SBitsArray& xored_bit_array);
    bool orWithSBitArray(const SBitsArray& ored_bit_array);
    bool andWithSBitArray(const SBitsArray& anded_bit_array);
    void xorWithSBitArray_fast(const SBitsArray& xored_bit_array);
    void orWithSBitArray_fast(const SBitsArray& ored_bit_array);
    void andWithSBitArray_fast(const SBitsArray& anded_bit_array);
    static SBitsArray xorWithSBitArray(const SBitsArray& bit_array1, const SBitsArray& bit_array2);
    static SBitsArray xorWithSBitArray_fast(const SBitsArray& bit_array1, const SBitsArray& bit_array2);
    static SBitsArray andWithSBitArray(const SBitsArray& bit_array1, const SBitsArray& bit_array2);
    static SBitsArray andWithSBitArray_fast(const SBitsArray& bit_array1, const SBitsArray& bit_array2);
    static SBitsArray orWithSBitArray(const SBitsArray& bit_array1, const SBitsArray& bit_array2);
    static SBitsArray orWithSBitArray_fast(const SBitsArray& bit_array1, const SBitsArray& bit_array2);
    bool xorWithSBitArrayFromPosition(const SBitsArray& xored_bit_array, int position);
    bool getHammingWeightMod2();
    bool getHammingWeightMod2_fast();
    void invertAllBits();
    void toZeroAllBits();
    void toOneAllBits();
    void cyclicShiftLeft(int shift);
    void cyclicShiftRight(int shift);
    void permute(Permutation &perm);
    void inverse_permute(Permutation &perm);

    inline qint64 getSizeInBits() const {return size_in_bits;}
    inline qint64 getSizeInQuint64() const {return size_in_quint64;}
    qint64 getNumberOfSignificantBytes() const;
    static qint64 getNumberOfSignificantBytes(quint64 size_in_bits_c_v);
    unsigned char *getCopyAsBytes();
    unsigned char *getDataAsBytes() const;
    bool save_to_file(std::string file_name);

private:
    quint64 *buffer_for_bits;
    qint64 size_in_bits;
    qint64 size_in_quint64;
    uint8_t hamming_weight_of_uint64(uint64_t val);

public:
    // LSBs on the left, MSBs on the right

    // operators
    friend ostream& operator<<(ostream& os, const SBitsArray& bit_array);
    friend bool operator==(const SBitsArray& bit_array_left, const SBitsArray& bit_array_right);
    SBitsArray& operator=(const SBitsArray& bit_array);
};

#endif // SBITSARRAY_H
