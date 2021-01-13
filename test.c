#include "ffsampling.h"
#include <assert.h>

char *seed;

struct KAT {
    double mu;
    double sigma;
    char *randombytes;
    int z;
};

// test vectors are taken from the falcon documentation https://falcon-sign.info/falcon.pdf

struct KAT test[] = {
    {-91.90471153063714, 1.7037990414754918, "0fc5442ff043d66e91d1eacac64ea5450a22941edc6c", -92},
    {-8.322564895434937, 1.7037990414754918, "f4da0f8d8444d1a77265c2ef6f98bbbb4bee7db8d9b3", -8},
    {-19.096516109216804, 1.7035823083824078, "db47f6d7fb9b19f25c36d6b9334d477a8bc0be68145d", -20},
    {-11.335543982423326, 1.7035823083824078, "ae41b4f5209665c74d00dcc1a8168a7bb516b3190cb42c1ded26cd52aed770eca7dd334e0547bcc3c163ce0b", -12},
    {7.9386734193997555, 1.6984647769450156, "31054166c1012780c603ae9b833cec73f2f41ca5807cc89c92158834632f9b1555", 8},
    {-28.990850086867255, 1.6984647769450156, "737e9d68a50a06dbbc6477", -30},
    {-9.071257914091655, 1.6980782114808988, "a98ddd14bf0bf22061d632", -10},
    {-43.88754568839566, 1.6980782114808988, "3cbf6818a68f7ab9991514", -41}
};

int main()
{
	double sigmin = 1.277833697;
	int z;

	init_RCDT();
    for (int i = 0; i < 8; i++)
    {
        seed = test[i].randombytes;
	    z = SamplerZ(test[i].mu, test[i].sigma, sigmin);
	    printf("z = %d\n", z);
        assert(test[i].z == z);
    }
	return (0);
}