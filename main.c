#include "ffsampling.h"

int main()
{
	float sigmin = 1.277833697;
	int z;
	float mu = -91.90471153063714;
	float sigma = 1.7037990414754918;
	char *test = "0fc5442ff043d66e91d1eacac64ea5450a22941edc6c";
	unsigned char stream[22];

	init_RCDT();
	for (size_t count = 0; count < 22; count++) {
        sscanf(test, "%2hhx", &stream[count]);
        test += 2;
    }

	for (int i = 0; i < 22; i++)
		printf("%02x ", stream[i]);
	printf("\n");
	z = SamplerZ(mu, sigma, sigmin, stream);
	printf("z = %d\n", z);
	return (0);
}