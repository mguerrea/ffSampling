#include "ffsampling.h"
#include <math.h>

/*
 * Map a string of arbitrary length to a point in the plane
 * @param message a string of arbitrary length
 * @param n the dimension of the ring
 * @return a point represented by a polynomial (coefficients)
*/
t_pol HashToPoint(char *message, int n)
{
	unsigned long hash = 5381;
	char *str = message;
    int c;
	t_pol point = {.len = n};
	point.coeffs = malloc(sizeof(double) * point.len);

	for (int i = 0; i < n; i++)
	{
		str = message;
    	while ((c = *str++))
        	hash = (((hash << 5) + hash) + c) % Q;
		point.coeffs[i] = hash;
	}
	return (point);
}

/*
 * Computes the gaussian sampling of t over a lattice
 * @param t a vector of two polynomials (FFT)
 * @param T a Falcon tree computed during key generation
 * @param params Falcon parameters, depens of the dimension
 * @return a vector of two polynomials (FFF)
*/
t_pol_fft *ffSampling(t_pol_fft t[2], t_tree *T, t_params params)
{
    t_pol_fft *z1, *z0, *z = malloc(sizeof(t_pol_fft) * 2);

    if (t[0].len == 1)
    {
        double sigma = creal(T->value.coeffs[0]);
        z[0] = new_pol(1);
        z[0].coeffs[0] = SamplerZ(creal(t[0].coeffs[0]), sigma, params.sigmin);
        z[1] = new_pol(1);
        z[1].coeffs[0] = SamplerZ(creal(t[1].coeffs[0]), sigma, params.sigmin);
        return(z);
    }
    t_pol_fft l = T->value;
    t_tree *T0 = T->leftchild, *T1 = T->rightchild;
    t_pol_fft t1[2], t0[2];
    split_fft(t[1], &(t1[0]), &(t1[1]));
    
    z1 = ffSampling(t1, T1, params);
    free(t1[0].coeffs);
    free(t1[1].coeffs);
    z[1] = merge_fft(z1[0], z1[1]);

    t_pol_fft t_prime = new_pol(t[0].len);
    sub_fft(&t_prime, t[1], z[1]);
    mul_fft(&t_prime, t_prime, l);
    add_fft(&t_prime, t_prime, t[0]);
    split_fft(t_prime, &(t0[0]), &(t0[1]));
    free(t_prime.coeffs);
    
    z0 = ffSampling(t0, T0, params);
    free(t0[0].coeffs);
    free(t0[1].coeffs);
    z[0] = merge_fft(z0[0], z0[1]);

    return (z);
}

/*
 * Computes a pseudo signature to illustrate the use of ffSampling
 * @param message a string representing the message to be signed
 * @param key a secret key generated with gen_sk
 * @param params Falcon parameters, depens of the dimension
 * @return a polynomial (coefficients) representing a pseudo signature
*/
t_pol pseudo_sign(char *message, t_sk key, t_params params)
{
    init_RCDT();
    t_pol point = HashToPoint(message, params.n);
    t_pol_fft point_fft = fft(point);
    t_pol_fft a = key.basis[0][0];
    t_pol_fft b = key.basis[0][1];
    t_pol_fft c = key.basis[1][0];
    t_pol_fft d = key.basis[1][1];
    t_pol_fft t[2], v[2];
    t_pol s[2];

    t[0] = new_pol(point_fft.len);
    t[1] = new_pol(point_fft.len);
    v[0] = new_pol(point_fft.len);
    v[1] = new_pol(point_fft.len);

    // we compute t which is a preimage of message (but not a short one)
    for (int i = 0; i < t[0].len; i++)
    {
        t[0].coeffs[i] = (point_fft.coeffs[i] * d.coeffs[i]) / (Q);
        t[1].coeffs[i] = (-point_fft.coeffs[i] * b.coeffs[i]) / (Q);
    }

    while (1)
    {

    // ffSampling gives us a random point z on our lattice close to t
    t_pol_fft *z = ffSampling(t, key.T, params);

    vect_mat_mul(v, z, key.basis);

    s[0] = ifft(v[0]);
    s[1] = ifft(v[1]);

    for (int i = 0; i < v[1].len; i++)
    {
        s[0].coeffs[i] = point.coeffs[i] - round(s[0].coeffs[i]);
        s[1].coeffs[i] = -round(s[1].coeffs[i]);
    }

    double norm = 0;
    for (int i = 0; i < s[0].len; i++)
        norm = norm + s[0].coeffs[i]*s[0].coeffs[i] + s[1].coeffs[i]*s[1].coeffs[i];
    if (norm <= params.bound)
        break;
    }
    
    // s is the distance between the message and the point z
    // s is short and we get s[0] + s[1] * h = message % q
    // we return s[1] which is the signature
    return (s[1]);
}

/*
 * Verify the pseudo signature created by pseudo sign.
 * If the norm of the signature is too large then the signature is rejected.
 * @param message a string representing the message which has been signed
 * @param sig a polynomial (coefficients) representing the signature of the message
 * @param h the public key used to verify the signature (generated in gen_sk)
 * @param params Falcon parameters, depens of the dimension
 * @return 1 if the signature is valid, 0 otherwise
*/
int pseudo_verify(char *message, t_pol sig, t_pol h, t_params params)
{
	double norm = 0;
    t_pol point = HashToPoint(message, params.n);
	t_pol tmp = mul_zq(sig, h);
	t_pol s1 = sub_zq(point, tmp);
    // we normalize the s1 coefficients around q/2
	for (int i = 0; i < s1.len; i++)
		s1.coeffs[i] = ((int)s1.coeffs[i] + (Q >> 1)) % Q - (Q >> 1);
	// we compute the norm of the signature
    for (int i = 0; i < s1.len; i++)
	{
		norm += sig.coeffs[i]*sig.coeffs[i];
		norm += s1.coeffs[i]*s1.coeffs[i];
	}
	if (norm <= params.bound)
		return (1);
	return (0);
}