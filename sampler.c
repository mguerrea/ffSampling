#include "ffsampling.h"

extern mpz_t RCDT[19];
extern uint64_t C[13];

int BaseSampler(unsigned char *stream)
{
	mpz_t u;
	gmp_randstate_t state;
	unsigned char buff[9];
	mpz_t seed;
	int z = 0;

	mpz_inits(u, seed, NULL);
	if (stream == NULL)
	{
	int fd = open("/dev/urandom", O_RDONLY);
	read(fd, buff, 9);
	close(fd);
	mpz_import(seed, 9, 1, 1, 0, 0, buff);
	}
	else
		mpz_import(seed, 9, 1, 1, 0, 0, stream);
	gmp_randinit_mt(state);
	gmp_randseed(state, seed);
	mpz_urandomb(u, state, 72);
	for (int i = 0; i < 18; i++)
		z += (mpz_cmp(RCDT[i], u) > 0);
	mpz_clears(u, seed, NULL);
	gmp_randclear(state);
	return (z);
}

uint64_t ApproxExp(float x, float ccs)
{
	uint64_t y = C[0];
	uint64_t z = (uint64_t)((float)(1UL << 63) * x);
	mpz_t tmp;
	uint64_t raw[2];

	mpz_init2(tmp, 126);
	for (int i = 1; i < 13; i++)
	{
		mpz_set_ui(tmp, y);
		mpz_mul_ui(tmp, tmp, z);
		mpz_export(raw, NULL, 1, 8, 0, 0, tmp);
		y = C[i] - ((raw[0] << 1) + (raw[1] >> 63));
	}
	z = (uint64_t)((float)(1UL << 63) * ccs);
	mpz_set_ui(tmp, y);
	mpz_mul_ui(tmp, tmp, z);
	mpz_export(raw, NULL, 1, 8, 0, 0, tmp);
	mpz_clear(tmp);
	return ((raw[0] << 1) + (raw[1] >> 63));
}

int BerExp(float x, float ccs, unsigned char *stream)
{
	unsigned char c;
	int w;
	int s = x / LN2;
	int r = x - s * LN2;

	s = MIN(s, 63);
	uint64_t z = (2 * ApproxExp(x, ccs) - 1) >> s;
	int i = 64;
	int fd = open("/dev/urandom", O_RDONLY);
	do {
		i = i - 8;
		if (stream)
			c = (*stream)++;
		else
			read(fd, &c, 1);
		w = c - ((z >> i) & 0xff);
	} while(w == 0 && i > 0);
	close(fd);
	return (w < 0);
}

int SamplerZ(float mu, float sigma, float sigmin, unsigned char *stream)
{
	float r = mu - (int)mu;
	float ccs = sigmin/sigma;
	int z, z0, b;
	float x;

	int fd = open("/dev/urandom", O_RDONLY);
	while (1)
	{
		printf("iter\n");
		z0 = BaseSampler(stream);
		if (stream)
		{
			b = stream[9];
			stream += 10;
		}
		else
			read(fd, &b, 1);
		b = b & 0x1;
		z = b + (2 * b - 1) * z0;
		x = (z - r)*(z - r) / (2*sigma*sigma) - z0*z0 / (2*SIGMAX*SIGMAX);
		if (BerExp(x, ccs, stream) == 1)
			return (z + (int)mu);
	}
	close(fd);
}