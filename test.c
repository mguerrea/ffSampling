#include "ffsampling.h"
#include "test.h"
#include <assert.h>

char *seed;
extern struct samplerZ_KAT samplerZ[];

int main()
{
	double sigmin = 1.277833697;
	int z;

	init_RCDT();
    for (int i = 0; i < 8; i++)
    {
        seed = samplerZ[i].randombytes;
	    z = SamplerZ(samplerZ[i].mu, samplerZ[i].sigma, sigmin);
	    printf("z = %d\n", z);
        assert(samplerZ[i].z == z);
    }
	return (0);
}