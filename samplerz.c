#include "ffsampling.h"

extern mpz_t RCDT[19];
extern uint64_t C[13];

int BaseSampler()
{
	mpz_t u;
	unsigned char buff[9];
	int z = 0;

	mpz_init(u);
	random_bytes(9, buff);
	mpz_import(u, 9, 1, 1, 0, 0, buff);
	for (int i = 0; i < 18; i++)
		z += (mpz_cmp(RCDT[i], u) > 0);
	mpz_clear(u);
	return (z);
}

uint64_t ApproxExp(double x, double ccs)
{
	uint64_t y = C[0];
	uint64_t z = (uint64_t)((double)(1UL << 63) * x);
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
	z = (uint64_t)((double)(1UL << 63) * ccs);
	mpz_set_ui(tmp, y);
	mpz_mul_ui(tmp, tmp, z);
	mpz_export(raw, NULL, 1, 8, 0, 0, tmp);
	mpz_clear(tmp);
	return ((raw[0] << 1) + (raw[1] >> 63));
}

int BerExp(double x, double ccs)
{
	unsigned char c;
	int w;
	int s = x / LN2;
	int r = x - s * LN2;

	s = MIN(s, 63);
	uint64_t z = (2 * ApproxExp(x, ccs) - 1) >> s;
	int i = 64;
	do {
		i = i - 8;
		random_bytes(1, &c);
		w = c - ((z >> i) & 0xff);
	} while(w == 0 && i > 0);
	return (w < 0);
}

int SamplerZ(double mu, double sigma, double sigmin)
{
	double r = mu - (int)mu + (mu < 0);
	double ccs = sigmin/sigma;
	int z, z0, b;
	double x;
	double dss = 1 / (2 * sigma * sigma);

	int fd = open("/dev/urandom", O_RDONLY);
	while (1)
	{
		z0 = BaseSampler();
		random_bytes(1, (unsigned char *)&b);
		b = b & 0x1;
		z = b + (2 * b - 1) * z0;
		x = ((z - r)*(z - r)) * dss - (z0 *z0) * INV_2SIGMA2;
		if (BerExp(x, ccs))
			return (z + (int)mu - (mu < 0));
	}
	close(fd);
}