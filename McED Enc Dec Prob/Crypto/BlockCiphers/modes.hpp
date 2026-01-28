#include <algorithm>
#include <functional>
#include <iostream>

typedef union {
    BYTE a[2];
	unsigned short b;
} GAMMA;

template <typename CipherType>
CFB_Mode<CipherType>::CFB_Mode(const CipherType & alg, const ByteBlock & init_vec) :
    algorithm(alg), iv(init_vec.deep_copy()), current_iv(init_vec.deep_copy()), current_offset(init_vec.size()/2)
{
    // nothing
}

template <typename CipherType>
void CFB_Mode<CipherType>::encrypt(const ByteBlock & src, ByteBlock & dst) const {
    auto blocks = split_blocks(src, CipherType::block_lenght);
    ByteBlock tmp;

    algorithm.encrypt(iv, tmp);
    xor_blocks(tmp, tmp, blocks[0]);
    blocks[0] = std::move(tmp);
    for(int i = 1; i < blocks.size(); i++) {
        algorithm.encrypt(blocks[i-1], tmp);
        xor_blocks(tmp, tmp, blocks[i]);
        blocks[i] = std::move(tmp);
    }
    dst = join_blocks(blocks);
}

template <typename CipherType>
unsigned short CFB_Mode<CipherType>::next_gamma_value() {
	ByteBlock tmp;
	GAMMA ret;
	
	if (2 * current_offset == current_iv.size())
	{
		algorithm.encrypt(current_iv, tmp);
		current_iv = std::move(tmp);
		current_offset = 0;
	}

	ret.a[0] = current_iv[2 * current_offset];
	ret.a[1] = current_iv[2 * current_offset + 1];
	current_offset++;

	return ret.b;
}

template <typename CipherType>
void CFB_Mode<CipherType>::decrypt_with_iv(const ByteBlock & src, ByteBlock & dst, const ByteBlock & iv_) const {
    auto blocks = split_blocks(src, CipherType::block_lenght);
	ByteBlock tmp;

	algorithm.encrypt(iv_, tmp);
	xor_blocks(tmp, blocks[0], tmp);
	swap(tmp, blocks[0]);
	for(int i = 1; i < blocks.size(); i++) {
		algorithm.encrypt(tmp, tmp);
		xor_blocks(tmp, blocks[i], tmp);
		swap(tmp, blocks[i]);
	}
	dst = join_blocks(blocks);
}

template <typename CipherType>
void CFB_Mode<CipherType>::decrypt(const ByteBlock & src, ByteBlock & dst) const {
	decrypt_with_iv(src, dst, iv);
}
