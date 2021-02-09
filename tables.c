#include "ffsampling.h"

const uint64_t C[13] = {
    0x00000004741183A3,
    0x00000036548CFC06,
    0x0000024FDCBF140A,
    0x0000171D939DE045,
    0x0000D00CF58F6F84,
    0x000680681CF796E3,
    0x002D82D8305B0FEA,
    0x011111110E066FD0,
    0x0555555555070F00,
    0x155555555581FF00,
    0x400000000002B400,
    0x7FFFFFFFFFFF4800,
    0x8000000000000000};

mpz_t RCDT[19];

void init_RCDT()
{
    mpz_init_set_str(RCDT[0], "3024686241123004913666", 10);
    mpz_init_set_str(RCDT[1], "1564742784480091954050", 10);
    mpz_init_set_str(RCDT[2], "636254429462080897535", 10);
    mpz_init_set_str(RCDT[3], "199560484645026482916", 10);
    mpz_init_set_str(RCDT[4], "47667343854657281903", 10);
    mpz_init_set_str(RCDT[5], "8595902006365044063", 10);
    mpz_init_set_str(RCDT[6], "1163297957344668388", 10);
    mpz_init_set_str(RCDT[7], "117656387352093658", 10);
    mpz_init_set_str(RCDT[8], "8867391802663976", 10);
    mpz_init_set_str(RCDT[9], "496969357462633", 10);
    mpz_init_set_str(RCDT[10], "20680885154299", 10);
    mpz_init_set_str(RCDT[11], "638331848991", 10);
    mpz_init_set_str(RCDT[12], "14602316184", 10);
    mpz_init_set_str(RCDT[13], "247426747", 10);
    mpz_init_set_str(RCDT[14], "3104126", 10);
    mpz_init_set_str(RCDT[15], "28824", 10);
    mpz_init_set_str(RCDT[16], "198", 10);
    mpz_init_set_str(RCDT[17], "1", 10);
	mpz_init_set_str(RCDT[18], "0", 10);
}

void free_RCDT()
{
    for (int i = 0; i < 19; i++)
        mpz_clear(RCDT[i]);
}