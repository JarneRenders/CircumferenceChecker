compiler=gcc
flags=-std=gnu11 -march=native -Wall -Wno-missing-braces -O3

# The 64-bit version of this program is faster but only supports graphs up to 64 vertices.
64bit: circumferenceChecker.c libs/readGraph6.c libs/bitset.h libs/hamiltonicityMethods.c
	$(compiler) -DUSE_64_BIT -o circumferenceChecker circumferenceChecker.c libs/readGraph6.c libs/hamiltonicityMethods.c $(flags)

# There are two different implementations of the 128-bit version. The array version generally performs faster.
128bit: circumferenceChecker.c libs/readGraph6.c libs/bitset.h libs/hamiltonicityMethods.c
	$(compiler) -DUSE_128_BIT -o circumferenceChecker-128 circumferenceChecker.c libs/readGraph6.c libs/hamiltonicityMethods.c $(flags)

192bit: circumferenceChecker.c libs/readGraph6.c libs/bitset.h libs/hamiltonicityMethods.c
	$(compiler) -DUSE_192_BIT -o circumferenceChecker-192 circumferenceChecker.c libs/readGraph6.c libs/hamiltonicityMethods.c $(flags)	

256bit: circumferenceChecker.c libs/readGraph6.c libs/bitset.h libs/hamiltonicityMethods.c
	$(compiler) -DUSE_256_BIT -o circumferenceChecker-256 circumferenceChecker.c libs/readGraph6.c libs/hamiltonicityMethods.c $(flags)	

profile: circumferenceChecker.c libs/readGraph6.c libs/bitset.h libs/hamiltonicityMethods.c
	$(compiler) -DUSE_64_BIT -o circumferenceChecker-pr circumferenceChecker.c libs/readGraph6.c libs/hamiltonicityMethods.c -std=gnu11 -march=native -Wall -Wno-missing-braces -g -pg -fsanitize=address

all: 64bit 128bit 192bit 256bit 

.PHONY: clean
clean:
	rm -f circumferenceChecker circumferenceChecker-128 circumferenceChecker-128 circumferenceChecker-192 circumferenceChecker-256

