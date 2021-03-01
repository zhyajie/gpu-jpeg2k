/* 
Copyright 2009-2013 Poznan Supercomputing and Networking Center

Authors:
Milosz Ciznicki miloszc@man.poznan.pl

GPU JPEG2K is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GPU JPEG2K is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with GPU JPEG2K. If not, see <http://www.gnu.org/licenses/>.
*/
//#define DEBUG_MQ

#ifdef DEBUG_MQ
 int l = 0;
 int Cstates[1200000];
#endif

#include "cpu_coeff_coder2.h"
#include "../../threads/ThreadPool.h"
extern "C" {
	#include "../../misc/memory_management.cuh"
	#include "../../print_info/print_info.h"
}
namespace CPU_JPEG2K
{
	#include "cpu_mq-coder.h"

 void SetMaskedBits(unsigned int &reg, unsigned int mask, unsigned int bits)
{
	reg = (reg & ~mask) | (bits & mask);
}

 void SetNthBit(unsigned int &reg, unsigned int n)
{
	SetMaskedBits(reg, 1 << n, 1 << n);
}

 void ResetNthBit(unsigned int &reg, unsigned int n)
{
	SetMaskedBits(reg, 1 << n, 0);
}

typedef struct
{
	CoefficientState tl;
	CoefficientState t;
	CoefficientState tr;
	
	CoefficientState l;
	CoefficientState c;
	CoefficientState r;
	
	CoefficientState bl;
	CoefficientState b;
	CoefficientState br;

	short pos;
} CtxWindow;

 void debug_print(float *val, int tid)
{
//	if(tid == 3)
//		printf("dist:%f\n", *val);
}


 void down(CodeBlockAdditionalInfo &info, CtxWindow &window, CoefficientState *coeffs)
{
	window.tr = coeffs[window.pos + 1 - info.width];
	window.r = coeffs[window.pos + 1];
	window.br = coeffs[window.pos + 1 + info.width];
}

 void up(CtxWindow &window, CoefficientState *coeffs)
{
	coeffs[window.pos - 1] = window.l;
}

 void shift(CtxWindow &window)
{
	window.tl = window.t; window.t = window.tr; window.tr = 0; // top layer
	window.l = window.c; window.c = window.r; window.r = 0; // middle layer
	window.bl = window.b; window.b = window.br; window.br = 0; // bottom layer
	window.pos += 1;
}

typedef int CtxReg;

#define TRIMASK 0x249 //((1 << 0) | (1 << 3) | (1 << 6) | (1 << 9))

 CtxReg buildCtxReg(CtxWindow &window, unsigned char bitoffset)
{
	CtxReg reg = 0;

	reg |= ((window.tl >> (bitoffset + 9)) & 1) << 0;
	reg |= ((window.t >> (bitoffset + 9)) & 1) << 1;
	reg |= ((window.tr >> (bitoffset + 9)) & 1) << 2;
	reg |= ((window.l >> (bitoffset)) & TRIMASK) << 3;
	reg |= ((window.c >> (bitoffset)) & TRIMASK) << 4;
	reg |= ((window.r >> (bitoffset)) & TRIMASK) << 5;
	reg |= ((window.bl >> (bitoffset)) & 1) << 15;
	reg |= ((window.b >> (bitoffset)) & 1) << 16;
	reg |= ((window.br >> (bitoffset)) & 1) << 17;
	
	return reg;
}

/* void updateWindow(CtxWindow &window, CtxReg reg, unsigned char bitoffset)
{
	window.tl = (window.tl & ~(1 << (bitoffset + 9)) | ((reg >> 0) & 1) << (bitoffset + 9);
	window.t = (window.t & ~(1 << (bitoffset + 9)) | ((reg >> 1) & 1) << (bitoffset + 9);
	window.tr = (window.tr & ~(1 << (bitoffset + 9)) | ((reg >> 2) & 1) << (bitoffset + 9);
	window.l = (window.l & ~(TRIMASK << (bitoffset)) | ((reg >> 3) & TRIMASK) << (bitoffset);
	window.c = (window.c & ~(TRIMASK << (bitoffset)) | ((reg >> 4) & TRIMASK) << (bitoffset);
	window.r = (window.r & ~(TRIMASK << (bitoffset)) | ((reg >> 5) & TRIMASK) << (bitoffset);
	window.bl = (window.bl & ~(1 << (bitoffset)) | ((reg >> 0) & 15) << (bitoffset);
	window.b = (window.b & ~(1 << (bitoffset)) | ((reg >> 0) & 16) << (bitoffset);
	window.br = (window.br & ~(1 << (bitoffset)) | ((reg >> 0) & 17) << (bitoffset);
}*/

	 unsigned char SPCXLUT[3][512] = {
		{
			0, 1, 3, 3, 1, 2, 3, 3, 5, 6, 7, 7, 6, 6, 7, 7, 0, 1, 3, 3,
			1, 2, 3, 3, 5, 6, 7, 7, 6, 6, 7, 7, 5, 6, 7, 7, 6, 6, 7, 7,
			8, 8, 8, 8, 8, 8, 8, 8, 5, 6, 7, 7, 6, 6, 7, 7, 8, 8, 8, 8,
			8, 8, 8, 8, 1, 2, 3, 3, 2, 2, 3, 3, 6, 6, 7, 7, 6, 6, 7, 7,
			1, 2, 3, 3, 2, 2, 3, 3, 6, 6, 7, 7, 6, 6, 7, 7, 6, 6, 7, 7,
			6, 6, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 6, 6, 7, 7, 6, 6, 7, 7,
			8, 8, 8, 8, 8, 8, 8, 8, 3, 3, 4, 4, 3, 3, 4, 4, 7, 7, 7, 7,
			7, 7, 7, 7, 3, 3, 4, 4, 3, 3, 4, 4, 7, 7, 7, 7, 7, 7, 7, 7,
			7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7,
			7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 3, 3, 4, 4, 3, 3, 4, 4,
			7, 7, 7, 7, 7, 7, 7, 7, 3, 3, 4, 4, 3, 3, 4, 4, 7, 7, 7, 7,
			7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8,
			7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 1, 2, 3, 3,
			2, 2, 3, 3, 6, 6, 7, 7, 6, 6, 7, 7, 1, 2, 3, 3, 2, 2, 3, 3,
			6, 6, 7, 7, 6, 6, 7, 7, 6, 6, 7, 7, 6, 6, 7, 7, 8, 8, 8, 8,
			8, 8, 8, 8, 6, 6, 7, 7, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8,
			2, 2, 3, 3, 2, 2, 3, 3, 6, 6, 7, 7, 6, 6, 7, 7, 2, 2, 3, 3,
			2, 2, 3, 3, 6, 6, 7, 7, 6, 6, 7, 7, 6, 6, 7, 7, 6, 6, 7, 7,
			8, 8, 8, 8, 8, 8, 8, 8, 6, 6, 7, 7, 6, 6, 7, 7, 8, 8, 8, 8,
			8, 8, 8, 8, 3, 3, 4, 4, 3, 3, 4, 4, 7, 7, 7, 7, 7, 7, 7, 7,
			3, 3, 4, 4, 3, 3, 4, 4, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
			7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7,
			8, 8, 8, 8, 8, 8, 8, 8, 3, 3, 4, 4, 3, 3, 4, 4, 7, 7, 7, 7,
			7, 7, 7, 7, 3, 3, 4, 4, 3, 3, 4, 4, 7, 7, 7, 7, 7, 7, 7, 7,
			7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7,
			7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8
		},
		{
			0, 1, 5, 6, 1, 2, 6, 6, 3, 3, 7, 7, 3, 3, 7, 7, 0, 1, 5, 6,
			1, 2, 6, 6, 3, 3, 7, 7, 3, 3, 7, 7, 3, 3, 7, 7, 3, 3, 7, 7,
			4, 4, 7, 7, 4, 4, 7, 7, 3, 3, 7, 7, 3, 3, 7, 7, 4, 4, 7, 7,
			4, 4, 7, 7, 1, 2, 6, 6, 2, 2, 6, 6, 3, 3, 7, 7, 3, 3, 7, 7,
			1, 2, 6, 6, 2, 2, 6, 6, 3, 3, 7, 7, 3, 3, 7, 7, 3, 3, 7, 7,
			3, 3, 7, 7, 4, 4, 7, 7, 4, 4, 7, 7, 3, 3, 7, 7, 3, 3, 7, 7,
			4, 4, 7, 7, 4, 4, 7, 7, 5, 6, 8, 8, 6, 6, 8, 8, 7, 7, 8, 8,
			7, 7, 8, 8, 5, 6, 8, 8, 6, 6, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8,
			7, 7, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8,
			7, 7, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8, 6, 6, 8, 8, 6, 6, 8, 8,
			7, 7, 8, 8, 7, 7, 8, 8, 6, 6, 8, 8, 6, 6, 8, 8, 7, 7, 8, 8,
			7, 7, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8,
			7, 7, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8, 1, 2, 6, 6,
			2, 2, 6, 6, 3, 3, 7, 7, 3, 3, 7, 7, 1, 2, 6, 6, 2, 2, 6, 6,
			3, 3, 7, 7, 3, 3, 7, 7, 3, 3, 7, 7, 3, 3, 7, 7, 4, 4, 7, 7,
			4, 4, 7, 7, 3, 3, 7, 7, 3, 3, 7, 7, 4, 4, 7, 7, 4, 4, 7, 7,
			2, 2, 6, 6, 2, 2, 6, 6, 3, 3, 7, 7, 3, 3, 7, 7, 2, 2, 6, 6,
			2, 2, 6, 6, 3, 3, 7, 7, 3, 3, 7, 7, 3, 3, 7, 7, 3, 3, 7, 7,
			4, 4, 7, 7, 4, 4, 7, 7, 3, 3, 7, 7, 3, 3, 7, 7, 4, 4, 7, 7,
			4, 4, 7, 7, 6, 6, 8, 8, 6, 6, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8,
			6, 6, 8, 8, 6, 6, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8,
			7, 7, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8,
			7, 7, 8, 8, 7, 7, 8, 8, 6, 6, 8, 8, 6, 6, 8, 8, 7, 7, 8, 8,
			7, 7, 8, 8, 6, 6, 8, 8, 6, 6, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8,
			7, 7, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8,
			7, 7, 8, 8, 7, 7, 8, 8, 7, 7, 8, 8
		},
		{
			0, 3, 1, 4, 3, 6, 4, 7, 1, 4, 2, 5, 4, 7, 5, 7, 0, 3, 1, 4,
			3, 6, 4, 7, 1, 4, 2, 5, 4, 7, 5, 7, 1, 4, 2, 5, 4, 7, 5, 7,
			2, 5, 2, 5, 5, 7, 5, 7, 1, 4, 2, 5, 4, 7, 5, 7, 2, 5, 2, 5,
			5, 7, 5, 7, 3, 6, 4, 7, 6, 8, 7, 8, 4, 7, 5, 7, 7, 8, 7, 8,
			3, 6, 4, 7, 6, 8, 7, 8, 4, 7, 5, 7, 7, 8, 7, 8, 4, 7, 5, 7,
			7, 8, 7, 8, 5, 7, 5, 7, 7, 8, 7, 8, 4, 7, 5, 7, 7, 8, 7, 8,
			5, 7, 5, 7, 7, 8, 7, 8, 1, 4, 2, 5, 4, 7, 5, 7, 2, 5, 2, 5,
			5, 7, 5, 7, 1, 4, 2, 5, 4, 7, 5, 7, 2, 5, 2, 5, 5, 7, 5, 7,
			2, 5, 2, 5, 5, 7, 5, 7, 2, 5, 2, 5, 5, 7, 5, 7, 2, 5, 2, 5,
			5, 7, 5, 7, 2, 5, 2, 5, 5, 7, 5, 7, 4, 7, 5, 7, 7, 8, 7, 8,
			5, 7, 5, 7, 7, 8, 7, 8, 4, 7, 5, 7, 7, 8, 7, 8, 5, 7, 5, 7,
			7, 8, 7, 8, 5, 7, 5, 7, 7, 8, 7, 8, 5, 7, 5, 7, 7, 8, 7, 8,
			5, 7, 5, 7, 7, 8, 7, 8, 5, 7, 5, 7, 7, 8, 7, 8, 3, 6, 4, 7,
			6, 8, 7, 8, 4, 7, 5, 7, 7, 8, 7, 8, 3, 6, 4, 7, 6, 8, 7, 8,
			4, 7, 5, 7, 7, 8, 7, 8, 4, 7, 5, 7, 7, 8, 7, 8, 5, 7, 5, 7,
			7, 8, 7, 8, 4, 7, 5, 7, 7, 8, 7, 8, 5, 7, 5, 7, 7, 8, 7, 8,
			6, 8, 7, 8, 8, 8, 8, 8, 7, 8, 7, 8, 8, 8, 8, 8, 6, 8, 7, 8,
			8, 8, 8, 8, 7, 8, 7, 8, 8, 8, 8, 8, 7, 8, 7, 8, 8, 8, 8, 8,
			7, 8, 7, 8, 8, 8, 8, 8, 7, 8, 7, 8, 8, 8, 8, 8, 7, 8, 7, 8,
			8, 8, 8, 8, 4, 7, 5, 7, 7, 8, 7, 8, 5, 7, 5, 7, 7, 8, 7, 8,
			4, 7, 5, 7, 7, 8, 7, 8, 5, 7, 5, 7, 7, 8, 7, 8, 5, 7, 5, 7,
			7, 8, 7, 8, 5, 7, 5, 7, 7, 8, 7, 8, 5, 7, 5, 7, 7, 8, 7, 8,
			5, 7, 5, 7, 7, 8, 7, 8, 7, 8, 7, 8, 8, 8, 8, 8, 7, 8, 7, 8,
			8, 8, 8, 8, 7, 8, 7, 8, 8, 8, 8, 8, 7, 8, 7, 8, 8, 8, 8, 8,
			7, 8, 7, 8, 8, 8, 8, 8, 7, 8, 7, 8, 8, 8, 8, 8, 7, 8, 7, 8,
			8, 8, 8, 8, 7, 8, 7, 8, 8, 8, 8, 8
		}
	};

 unsigned char getSPCX(CtxReg c, unsigned char i, unsigned char subband)
{
	return SPCXLUT[subband][(c >> (3 * i)) & 0x1FF];
}

	/* sign context in the following format
		index:
			first (MSB) bit V0 significance (1 significant, 0 insignificant)
			second bit V0 sign (0 positive, 1 negative)

			next 2 bits same for H0
			next 2 bits same for H1
			next 2 bits same for V1
			
		value:
			the response contains two pieces of information
			1. context label on the 4 least significant bits
			2. XORbit on the 5-th bit from the end (5-th least significant bit)
	*/

	 unsigned char signcxlut[256] = {
		 9,  9, 10, 26,  9,  9, 10, 26, 12, 12, 13, 11, 28, 28, 27, 29,  9,  9, 10, 26,
		 9,  9, 10, 26, 12, 12, 13, 11, 28, 28, 27, 29, 12, 12, 13, 11, 12, 12, 13, 11,
		12, 12, 13, 11,  9,  9, 10, 26, 28, 28, 27, 29, 28, 28, 27, 29,  9,  9, 10, 26,
		28, 28, 27, 29,  9,  9, 10, 26,  9,  9, 10, 26, 12, 12, 13, 11, 28, 28, 27, 29,
		 9,  9, 10, 26,  9,  9, 10, 26, 12, 12, 13, 11, 28, 28, 27, 29, 12, 12, 13, 11,
		12, 12, 13, 11, 12, 12, 13, 11,  9,  9, 10, 26, 28, 28, 27, 29, 28, 28, 27, 29,
		 9,  9, 10, 26, 28, 28, 27, 29, 10, 10, 10,  9, 10, 10, 10,  9, 13, 13, 13, 12,
		27, 27, 27, 28, 10, 10, 10,  9, 10, 10, 10,  9, 13, 13, 13, 12, 27, 27, 27, 28,
		13, 13, 13, 12, 13, 13, 13, 12, 13, 13, 13, 12, 10, 10, 10,  9, 27, 27, 27, 28,
		27, 27, 27, 28, 10, 10, 10,  9, 27, 27, 27, 28, 26, 26,  9, 26, 26, 26,  9, 26,
		11, 11, 12, 11, 29, 29, 28, 29, 26, 26,  9, 26, 26, 26,  9, 26, 11, 11, 12, 11,
		29, 29, 28, 29, 11, 11, 12, 11, 11, 11, 12, 11, 11, 11, 12, 11, 26, 26,  9, 26,
		29, 29, 28, 29, 29, 29, 28, 29, 26, 26,  9, 26, 29, 29, 28, 29
	};

 unsigned char getSICX(CtxReg sig, CtxReg sign, unsigned char i)
{
	return signcxlut[
			((sig >> (i * 3)) & 0xAA) |
			(((sign >> (i * 3)) & 0xAA) >> 1)
		];
}

 unsigned char getMRCX(CtxReg sig, CoefficientState local, unsigned char i)
{
	if((local >> (12 + 3 * i)) & 1)
		return 16;
	else if(((sig >> (3 * i)) & 0x1EF) == 0)
		return 14;
	else
		return 15;
}

template<class T>  T min(T val1, T val2)
{
	if(val1 > val2)
	{
		return val2;
	} else
	{
		return val1;
	}
}

template<class T>  T max(T val1, T val2)
{
	if(val1 < val2)
	{
		return val2;
	} else
	{
		return val1;
	}
}

 float distWeights[2][4][4] = {
{//Lossless
//		LH,      HL,      HH,     LLend
	{0.1000f, 0.1000f, 0.0500f, 1.0000f},  //level 0 = biggest subbands (unimportant)
	{0.2000f, 0.2000f, 0.1000f, 1.0000f},  //      1
	{0.4000f, 0.4000f, 0.2000f, 1.0000f},  //      2
	{0.8000f, 0.8000f, 0.4000f, 1.0000f}   //      3 = smallest, contains LL
}, {//Lossy
/*	{ 0.0010f, 0.0010f, 0.0005f, 1.0000f},
	{ 0.1000f, 0.1000f, 0.0250f, 1.0000f},
	{ 0.3000f, 0.3000f, 0.0800f, 1.0000f},
	{ 0.8000f, 0.8000f, 0.4000f, 1.0000f}*/
	{0.0100f, 0.0100f, 0.0050f, 1.0000f},
	{0.2000f, 0.2000f, 0.1000f, 1.0000f},
	{0.4000f, 0.4000f, 0.2000f, 1.0000f},
	{0.8000f, 0.8000f, 0.4000f, 1.0000f}
} };

 float getDISW(CodeBlockAdditionalInfo *info)
{
	return distWeights[info->compType][min<byte>(info->dwtLevel, 3)][info->subband] * info->stepSize * info->stepSize / ((float)(info->width * info->height));
}

class RLEncodeFunctor {
public:
	 char operator()(CtxWindow window, MQEncoder &enc)
	{
		char rest = 0;

		if((window.c & TRIMASK) == 0)
		{
			mqEncode(enc, 0, CX_RUN);
			rest = -2;
		}
		else
		{
			while(rest < 4 && ((window.c >> (3 * rest)) & 1) == 0)
				rest++;

			mqEncode(enc, 1, CX_RUN);
			mqEncode(enc, rest >> 1, CX_UNI);
			mqEncode(enc, rest & 1, CX_UNI);
		}

		return rest;
	}
};

class RLDecodeFunctor {
public:
	 char operator()(CtxWindow &window, MQDecoder &dec)
	{
		char rest = 0;

		if(mqDecode(dec, CX_RUN) == 0)
		{
			rest = -2;
		}
		else
		{
			rest = mqDecode(dec, CX_UNI) & 1;
			rest <<= 1;
			rest |= mqDecode(dec, CX_UNI) & 1;

			window.c |= 1 << (3 * rest);
		}

		return rest;
	}
};

class SigEncodeFunctor {
public:
	 void operator()(CtxWindow &window, CtxReg &sig, MQEncoder &enc, int stripId, int subband)
	{
		mqEncode(enc, (window.c >> (3 * stripId)) & 1, getSPCX(sig, stripId, subband));
	}
};

class SigDecodeFunctor {
public:
	 void operator()(CtxWindow &window, CtxReg sig, MQDecoder &dec, int stripId, int subband)
	{
		window.c |= mqDecode(dec, getSPCX(sig, stripId, subband)) << (3 * stripId);
	}
};

class SignEncodeFunctor
{
public:
	 void operator()(CtxWindow &window, CtxReg &sig, MQEncoder &enc, int stripId)
	{
		unsigned char cx = getSICX(sig, buildCtxReg(window, 13), stripId);

		mqEncode(enc, (short) (((window.c >> (13 + 3 * stripId)) & 1) ^ ((cx >> 4) & 1)), cx & 0xF);
	}
};

class SignDecodeFunctor
{
public:
	 void operator()(CtxWindow &window, CtxReg sig, MQDecoder &dec, int stripId)
	{
		unsigned char cx = getSICX(sig, buildCtxReg(window, 13), stripId);

		window.c |= (mqDecode(dec, cx & 0xF) ^ ((cx >> 4) & 1) & 1) << (13 + 3 * stripId);
	}	
};

template <class RLCodingFunctor, class SigCodingFunctor, class SignCodingFunctor, typename MQCoderStateType>
class CleanUpPassFunctor
{
public:
	 void operator()(const CodeBlockAdditionalInfo &info, CtxWindow &window, MQCoderStateType &mq, float *sum_dist, unsigned char bitplane)
	{
		char rest;

		CtxReg sig = buildCtxReg(window, 1); // significance context

		rest = -1;
		if((window.c & (TRIMASK << 14)) == 0 && sig == 0) // all contexts in stripe are equal to zero
		{
			rest = RLCodingFunctor()(window, mq);
			if(rest == -2)
				return;
		}

		for(int k = 0; k < 4; k++)
		{
			if(/*	((window.c >> ( 1 + 3 * k)) & 1) == 0 &&   // check if coefficient is non-significant (sigma)
				((window.c >> ( 2 + 3 * k)) & 1) == 0 &&   // check if coefficient hasn't been coded already (pi)
				((window.c >> (14 + 3 * k)) & 1) == 0)    // forbidden state indicating out of bounds (late sigma)*/
				((window.c >> (3 * k)) & 0x4006) == 0)
			{
				if(rest >= 0)
					rest--;
				else
					SigCodingFunctor()(window, sig, mq, k, info.subband);
			
				if((window.c >> (3 * k)) & 1) // check if magnitude is 1
				{
					*sum_dist -= (float)((1<<bitplane)*(1<<bitplane));
					//debug_print(sum_dist, threadIdx.x);
//					if(blockIdx.x * blockDim.x + threadIdx.x == 0)
//					printf("clu:%f tid:%d\n", *sum_dist, blockIdx.x * blockDim.x + threadIdx.x);
					SetNthBit(window.c, 1 + 3 * k); // set k-th significant state
					sig = buildCtxReg(window, 1); // rebuild significance register

					SignCodingFunctor()(window, sig, mq, k);
				}
			}
		}
	}
};

template <class SigCodingFunctor, class SignCodingFunctor, typename MQCoderStateType>
class SigPropPassFunctor {
public:
 void operator()(const CodeBlockAdditionalInfo &info, CtxWindow &window, MQCoderStateType &mq, float *sum_dist, unsigned char bitplane)
{
	CtxReg sig = buildCtxReg(window, 1); // build significance context register

	for(int i = 0; i < 4; i++)
	{
		// not significant with non-zero context
		if(/*	((window.c >> (1 + 3 * i)) & 1) == 0 &&
			((window.c >> (14 + 3 * i)) & 1) == 0 && // out of bounds
			getSPCX(sig, i, subband) > 0)*/
			(((window.c >> (3 * i)) & 0x4002) == 0) &&
			((sig >> (3 * i)) & 0x1EF) != 0)
		{
			SigCodingFunctor()(window, sig, mq, i, info.subband);

			// if magnitude bit is one
			if((window.c >> (3 * i)) & 1)
			{
				*sum_dist -= (float)((1<<bitplane)*(1<<bitplane));
				//debug_print(sum_dist, threadIdx.x);
//				if(blockIdx.x * blockDim.x + threadIdx.x == 0)
//				printf("sig:%f tid:%d\n", *sum_dist, blockIdx.x * blockDim.x + threadIdx.x);
				SetNthBit(window.c, 1 + (3 * i));
				sig = buildCtxReg(window, 1); // rebuild

				SignCodingFunctor()(window, sig, mq, i);
			}

			// set pi (already coded)
			SetNthBit(window.c, 2 + (3 * i));
		}
		else
			// unset pi (already coded)
			ResetNthBit(window.c, 2 + (3 * i));
	}
}
};

class MagRefEncodeFunctor {
public:
	 void operator()(MQEncoder &enc, CtxWindow &window, int stripId)
	{
		mqEncode(enc, (window.c >> (3 * stripId)) & 1, getMRCX(buildCtxReg(window, 1), window.c, stripId));
	}
};

class MagRefDecodeFunctor {
public:
	 void operator()(MQDecoder &dec, CtxWindow &window, int stripId)
	{
		window.c |= (mqDecode(dec, getMRCX(buildCtxReg(window, 1), window.c, stripId)) << (3 * stripId));
	}
};

template <class MagRefCodingFunctor, typename MQCoderStateType>
class MagRefPassFunctor {
public:
 void operator()(const CodeBlockAdditionalInfo &info, CtxWindow &window, MQCoderStateType &mq, float *sum_dist, unsigned char bitplane)
{
	for(int i = 0; i < 4; i++)
	{
		if(//csSignificant(st) && !csAlreadyCoded(st) && not out of bounds
			((window.c >> (3 * i)) & 0x4006) == 0x2)
		{
			*sum_dist -= (float)((1<<bitplane)*(1<<bitplane));
			//debug_print(sum_dist, threadIdx.x);
//			if(blockIdx.x * blockDim.x + threadIdx.x == 0)
//			printf("mgr:%f tid:%d\n", *sum_dist, blockIdx.x * blockDim.x + threadIdx.x);
			MagRefCodingFunctor()(mq, window, i);
			SetNthBit(window.c, 3 * i + 12);
		}
	}
}
};

 void initCoeffs(const CodeBlockAdditionalInfo &info, CoefficientState *coeffs)
{
	unsigned char signOffset = sizeof(int) * 8 - 1;

	for(int i = 0; i < info.width; i++)
		for(int j = 0; j < info.stripeNo; j++)
		{
			CoefficientState st = 0;
			int c;

			for(int k = 0; k < 4; k++)
				if(4 * j + k < info.height)
				{
					c = info.coefficients[(4 * j + k) * info.nominalWidth + i];
					//Cstates[l++] = (4 * j + k) * info.nominalWidth + i;
					st |= (((c >> signOffset) & 1) << (13 + 3 * k));
				}
				else
					st |= (1 << (14 + 3 * k));

			coeffs[j * info.width + i] = st;
		}
}

 void initDecodingCoeffs(const CodeBlockAdditionalInfo &info, CoefficientState *coeffs)
{
	for(int i = 0; i < info.width; i++)
		for(int j = 0; j < info.stripeNo; j++)
		{
			CoefficientState st = 0;

			for(int k = 0; k < 4; k++)
				if(4 * j + k < info.height)
					info.coefficients[(4 * j + k) * info.nominalWidth + i] = 0;
				else
					st |= (1 << (14 + 3 * k));

			coeffs[j * info.width + i] = st;
		}
}

 void uploadSigns(const CodeBlockAdditionalInfo &info, CoefficientState *coeffs)
{
	unsigned char signOffset = sizeof(int) * 8 - 1;

	for(int i = 0; i < info.width; i++)
		for(int j = 0; j < info.stripeNo; j++)
		{
			CoefficientState st = coeffs[j * info.width + i];

			for(int k = 0; k < 4; k++)
				if(((st >> (14 + 3 * k)) & 1) == 0)
					info.coefficients[(4 * j + k) * info.nominalWidth + i] |= (((st >> (13 + 3 * k)) & 1) << signOffset);

			coeffs[j * info.width + i] = st;
		}
}

 void fillMags(const CodeBlockAdditionalInfo &info, CoefficientState *coeffs, int bitplane)
{
	for(int i = 0; i < info.width; i++)
		for(int j = 0; j < info.stripeNo; j++)
		{
			CoefficientState st = coeffs[j * info.width + i];

			// clear magnitudes and already coded flags
			st &= ~(TRIMASK | (TRIMASK << 2));
			//st |= ((st & (TRIMASK << 1)) << 11);

			for(int k = 0; k < 4; k++)
				if(((st >> (14 + 3 * k)) & 1) == 0)
					st |= ((info.coefficients[(4 * j + k) * info.nominalWidth + i] >> bitplane) & 1) << (3 * k);

			coeffs[j * info.width + i] = st;
		}
}

 void uploadMags(const CodeBlockAdditionalInfo &info, CoefficientState *coeffs, int bitplane)
{
	for(int i = 0; i < info.width; i++)
		for(int j = 0; j < info.stripeNo; j++)
		{
			CoefficientState st = coeffs[j * info.width + i];

			for(int k = 0; k < 4; k++)
				if(((st >> (14 + 3 * k)) & 1) == 0)
					info.coefficients[(4 * j + k) * info.nominalWidth + i] |= (((st >> (3 * k)) & 1) << bitplane);

			// clear magnitudes and already coded flags
			st &= ~(TRIMASK | (TRIMASK << 2));

			coeffs[j * info.width + i] = st;
		}
}

 void clearWindow(CtxWindow &w)
{
	w.bl = 0;
	w.b = 0;
	w.br = 0;

	w.l = 0;
	w.c = 0;
	w.r = 0;

	w.tl = 0;
	w.t = 0;
	w.tr = 0;
}

template <class PassFunctor, typename MQCoderStateType>
 void BITPLANE_WINDOW_SCAN(CodeBlockAdditionalInfo &info, CoefficientState *coeffs, MQCoderStateType &enc, float *sum_dist, unsigned char bitplane) {
	CtxWindow window;

	window.pos = -1;

	for(int j = 0; j < info.stripeNo; j++)
	{
		clearWindow(window);
		down(info, window, coeffs);
		shift(window);
		down(info, window, coeffs);
	
		PassFunctor()(info, window, enc, sum_dist, bitplane);

		for(int k = 0; k < info.width - 2; k++)
		{
			shift(window);
			down(info, window, coeffs);
			PassFunctor()(info, window, enc, sum_dist, bitplane);
			up(window, coeffs);
		}

		shift(window);
		PassFunctor()(info, window, enc, sum_dist, bitplane);
		up(window, coeffs);
		shift(window);
		up(window, coeffs);

		window.pos--;
	}
}

class PCRD_EmptyFunctor
{
public:
	 void operator()(MQEncoder state, MQEncoder *states, unsigned char &stateId, float sum_dist, PcrdCodeblock *pcrdCodeblock)
	{
	}
};

class PCRD_CollectMQStatesFunctor
{
public:
	 void operator()(MQEncoder state, MQEncoder *states, unsigned char &stateId, float sum_dist, PcrdCodeblock *pcrdCodeblock)
	{
		states[stateId++] = state;
		pcrdCodeblock[stateId].dist = sum_dist;
	}
};

class CollectMQStatesFunctor
{
public:
	 void operator()(MQEncoder state, MQEncoder *states, unsigned char &stateId, float sum_dist, PcrdCodeblock *pcrdCodeblock)
	{
		states[stateId++] = state;
	}
};

 void decode(CoefficientState *coeffs, CodeBlockAdditionalInfo &info, byte *in)
{
	MQDecoder mqdec;
	mqInitDec(mqdec, in, info.length);

	float sum_dist = 0.0f;

	if(info.significantBits > 0)
	{
		mqResetDec(mqdec);

		initDecodingCoeffs(info, coeffs);

		BITPLANE_WINDOW_SCAN
		<CleanUpPassFunctor<RLDecodeFunctor, SigDecodeFunctor, SignDecodeFunctor, MQDecoder>, MQDecoder>
			(info, coeffs, mqdec, &sum_dist, 0);

		uploadMags(info, coeffs, 30 - info.magbits + info.significantBits);

		for(unsigned char i = 1; i < info.significantBits; i++)
		{
			BITPLANE_WINDOW_SCAN
			<SigPropPassFunctor<SigDecodeFunctor, SignDecodeFunctor, MQDecoder>, MQDecoder>
				(info, coeffs, mqdec, &sum_dist, 0);

			BITPLANE_WINDOW_SCAN
			<MagRefPassFunctor<MagRefDecodeFunctor, MQDecoder>, MQDecoder>
				(info, coeffs, mqdec, &sum_dist, 0);

			BITPLANE_WINDOW_SCAN
			<CleanUpPassFunctor<RLDecodeFunctor, SigDecodeFunctor, SignDecodeFunctor, MQDecoder>, MQDecoder>
				(info, coeffs, mqdec, &sum_dist, 0);

			uploadMags(info, coeffs, 30 - info.magbits - i + info.significantBits);
		}

		uploadSigns(info, coeffs);
		//mqDecode(mqdec, CX_UNI);
	}
	else
	{
		for(int i = 0; i < info.height; i++)
			for(int j = 0; j < info.width; j++)
				info.coefficients[i * info.nominalWidth + j] = 0;
	}
}



 void cpu_decode(CoefficientState *coeffBuffors, byte *inbuf, int maxThreadBufforLength, CodeBlockAdditionalInfo *infos, int codeBlocks)
{
	zl::ThreadPool threadPool(12);
	
    for(int i=0; i< codeBlocks;i++)
    {
        CodeBlockAdditionalInfo info = infos[i];
		threadPool.addTask(std::bind(&decode, coeffBuffors + info.magconOffset, info, inbuf + i * maxThreadBufforLength));
  //   	long int start_decode;
  //   	static long int total_time=0;
		// start_decode = start_measure();
		// decode(coeffBuffors + info.magconOffset, info, inbuf + i * maxThreadBufforLength);
		// total_time += stop_measure(start_decode);
		// printf("decode %ld us\n", stop_measure(start_decode));
		// printf("total_time %ld ms\n", total_time/1000);

    }
	 while(1)  
    {  
        //printf("there are still %d tasks need to process\n", threadPool.size());  
        if (threadPool.size() == 0)
        {  
            threadPool.stop();
            printf("Now I will exit from main\n"); 
            break;   
        }   
    }  
	
}



void launch_decode( CoefficientState *coeffBuffors, byte *inbuf, int maxThreadBufforLength, CodeBlockAdditionalInfo *infos, int codeBlocks)
{
	cpu_decode(coeffBuffors, inbuf, maxThreadBufforLength, infos, codeBlocks);
}

}
