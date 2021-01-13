#include "ffsampling.h"

#ifndef TEST

void random_bytes(int nb, unsigned char *buff)
{
    int fd = open("/dev/urandom", O_RDONLY);
    if (fd == -1)
        exit(1);
    read(fd, buff, nb);
    close(fd);
}

#else

extern char *seed;

void random_bytes(int nb, unsigned char *buff)
{
    for (size_t count = 0; count < nb; count++) {
        sscanf(seed, "%2hhx", &buff[count]);
        seed += 2;
    }
}

#endif