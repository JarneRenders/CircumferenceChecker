#ifndef BITSETCHOOSER
#define BITSETCHOOSER

#ifdef USE_64_BIT
	#include "bitset64Vertices.h"
	#define BITSETSIZE 64

#elif defined(USE_128_BIT)
	#include "bitset128Vertices.h"
	#define BITSETSIZE 128

#elif defined(USE_192_BIT)
	#include "bitset192Vertices.h"
	#define BITSETSIZE 192

#elif defined(USE_256_BIT)
	#include "bitset256Vertices.h"
	#define BITSETSIZE 256

#endif

#endif
