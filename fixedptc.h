#ifndef _FIXEDPTC_H_
#define _FIXEDPTC_H_

/*
 * fixedptc.h is a 32-bit or 64-bit fixed point numeric library.
 *
 * The symbol FIXEDPT_BITS, if defined before this library header file
 * is included, determines the number of bits in the data type (its "width").
 * The default width is 32-bit (FIXEDPT_BITS=32) and it can be used
 * on any recent C99 compiler.
 *
 * Changes were made in order to support 64-bit precision (FIXEDPT_BITS=64)
 * without the support for __int128_t and __uint128_t - 128-bit "long long"
 * types.
 * fixedpt_xmul & fixedpt_xdiv are no longer available when FIXEDPT_BITS=64.
 * fixedpt_str & fixedpt_cstr are no longer available when FIXEDPT_BITS=64.
 *
 * Changes were made in order to support the omission of stdint.h include 
 *
 * The FIXEDPT_WBITS symbols governs how many bits are dedicated to the
 * "whole" part of the number (to the left of the decimal point). The larger
 * this width is, the larger the numbers which can be stored in the fixedpt
 * number. The rest of the bits (available in the FIXEDPT_FBITS symbol) are
 * dedicated to the fraction part of the number (to the right of the decimal
 * point).
 *
 * Since the number of bits in both cases is relatively low, many complex
 * functions (more complex than div & mul) take a large hit on the precision
 * of the end result because errors in precision accumulate.
 * This loss of precision can be lessened by increasing the number of
 * bits dedicated to the fraction part, but at the loss of range.
 *
 * Adventurous users might utilize this library to build two data types:
 * one which has the range, and one which has the precision, and carefully
 * convert between them (including adding two number of each type to produce
 * a simulated type with a larger range and precision).
 *
 * The ideas and algorithms have been cherry-picked from a large number
 * of previous implementations available on the Internet.
 * Tim Hartrick has contributed cleanup and 64-bit support patches.
 *
 * == Special notes for the 32-bit precision ==
 * Signed 32-bit fixed point numeric library for the 24.8 format.
 * The specific limits are -8388608.999... to 8388607.999... and the
 * most precise number is 0.00390625. In practice, you should not count
 * on working with numbers larger than a million or to the precision
 * of more than 2 decimal places. Make peace with the fact that PI
 * is 3.14 here. :)
 */

/*-
 * Copyright (c) 2010-2012 Ivan Voras <ivoras@freebsd.org>
 * Copyright (c) 2012 Tim Hartrick <tim@edgecast.com>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#ifndef FIXEDPT_BITS
#define FIXEDPT_BITS	32
#endif

#ifndef OMIT_STDINT
#include <stdint.h>
#endif

#if FIXEDPT_BITS == 32
typedef int32_t fixedpt;
typedef	int64_t	fixedptd;
typedef	uint32_t fixedptu;
typedef	uint64_t fixedptud;
#elif FIXEDPT_BITS == 64
typedef int64_t fixedpt;
typedef	int64_t fixedptd;
typedef	uint64_t fixedptu;
typedef	int64_t fixedptud;
#else
#error "FIXEDPT_BITS must be equal to 32 or 64"
#endif

#ifndef FIXEDPT_WBITS
#define FIXEDPT_WBITS	24
#endif

#if FIXEDPT_WBITS >= FIXEDPT_BITS
#error "FIXEDPT_WBITS must be less than or equal to FIXEDPT_BITS"
#endif

#define FIXEDPT_VCSID "$Id$"

#define FIXEDPT_FBITS	(FIXEDPT_BITS - FIXEDPT_WBITS)
#define FIXEDPT_FMASK	(((fixedpt)1 << FIXEDPT_FBITS) - 1)

#define fixedpt_rconst(R) ((fixedpt)((R) * FIXEDPT_ONE + ((R) >= 0 ? 0.5 : -0.5)))
#define fixedpt_fromint(I) ((fixedptd)(I) << FIXEDPT_FBITS)
#define fixedpt_toint(F) ((F) >> FIXEDPT_FBITS)
#define fixedpt_add(A,B) ((A) + (B))
#define fixedpt_sub(A,B) ((A) - (B))
#if FIXEDPT_BITS != 64
#define fixedpt_xmul(A,B)						\
	((fixedpt)(((fixedptd)(A) * (fixedptd)(B)) >> FIXEDPT_FBITS))
#define fixedpt_xdiv(A,B)						\
	((fixedpt)(((fixedptd)(A) << FIXEDPT_FBITS) / (fixedptd)(B)))
#endif
#define fixedpt_fracpart(A) ((fixedpt)(A) & FIXEDPT_FMASK)

#define FIXEDPT_ONE	((fixedpt)((fixedpt)1 << FIXEDPT_FBITS))
#define FIXEDPT_ONE_HALF (FIXEDPT_ONE >> 1)
#define FIXEDPT_TWO	(FIXEDPT_ONE + FIXEDPT_ONE)
#define FIXEDPT_PI	fixedpt_rconst(3.14159265358979323846)
#define FIXEDPT_TWO_PI	fixedpt_rconst(2 * 3.14159265358979323846)
#define FIXEDPT_HALF_PI	fixedpt_rconst(3.14159265358979323846 / 2)
#define FIXEDPT_E	fixedpt_rconst(2.7182818284590452354)

#define fixedpt_abs(A) ((A) < 0 ? -(A) : (A))

/* fixedpt is meant to be usable in environments without floating point support
 * (e.g. microcontrollers, kernels), so we can't use floating point types directly.
 * Putting them only in macros will effectively make them optional. */
#define fixedpt_tofloat(T) ((float) ((T)*((float)(1)/(float)(1 << FIXEDPT_FBITS))))

#if FIXEDPT_BITS == 64
/*
 * The following code is taken from the following post to allow multiplication
 * and division of 128-bit (64 hi, 64 lo) numbers by 64 bit quotient without
 * the use of 128-bit "long long" types.
 *
 * https://codereview.stackexchange.com/questions/67962/mostly-portable-128-by-64-bit-division
 *
 * Changes were made to port c++ to c.
 */

int bsr(uint64_t x)
{
	uint64_t y;
	uint64_t r;

	r = (x > 0xFFFFFFFF) << 5; x >>= r;
	y = (x > 0xFFFF) << 4; x >>= y; r |= y;
	y = (x > 0xFF) << 3; x >>= y; r |= y;
	y = (x > 0xF) << 2; x >>= y; r |= y;
	y = (x > 0x3) << 1; x >>= y; r |= y;

	return (int)(r | (x >> 1));
}

uint64_t mul(uint64_t a, uint64_t b, uint64_t *y)
{
	uint64_t a_lo = a & 0x00000000FFFFFFFF;
	uint64_t a_hi = a >> 32;

	uint64_t b_lo = b & 0x00000000FFFFFFFF;
	uint64_t b_hi = b >> 32;

	uint64_t c0 = a_lo * b_lo;
	uint64_t c1 = a_hi * b_lo;
	uint64_t c2 = a_hi * b_hi;

	uint64_t u1 = c1 + (a_lo * b_hi);
        uint64_t u0;

	if (u1 < c1) {
		c2 += 1LL << 32;
	}

	u0 = c0 + (u1 << 32);
	if (u0 < c0) {
		++c2;
	}

	*y = c2 + (u1 >> 32);

	return u0;
}

uint64_t div(uint64_t a_lo, uint64_t a_hi, uint64_t b, uint64_t *r)
{
	uint64_t p_lo;
	uint64_t p_hi;
	uint64_t q = 0;

	uint64_t r_hi = a_hi;
	uint64_t r_lo = a_lo;

        int s = 0;
        int t;

        uint64_t b_hi;
        uint64_t q_hat;
        uint64_t u_hi;
        uint64_t u_lo;
        uint64_t w_lo;
        uint64_t w_hi;

	if (0 == (b >> 63)) {

		// Normalize so quotient estimates are
		// no more than 2 in error.

		// Note: If any bits get shifted out of
		// r_hi at this point, the result would
		// overflow.

		s = 63 - bsr(b);
		t = 64 - s;

		b <<= s;
		r_hi = (r_hi << s) | (r_lo >> t);
		r_lo <<= s;
	}

	b_hi = b >> 32;

	/*
	The first full-by-half division places b
	across r_hi and r_lo, making the reduction
	step a little complicated.

	To make this easier, u_hi and u_lo will hold
	a shifted image of the remainder.

	[u_hi||    ][u_lo||    ]
	[r_hi||    ][r_lo||    ]
	[ b  ||    ]
	[p_hi||    ][p_lo||    ]
	|
	V
	[q_hi||    ]
	*/

	q_hat = r_hi / b_hi;

	p_lo = mul(b, q_hat, &p_hi);

	u_hi = r_hi >> 32;
	u_lo = (r_hi << 32) | (r_lo >> 32);

	// r -= b*q_hat
	//
	// At most 2 iterations of this...
	while (
		(p_hi > u_hi) ||
		((p_hi == u_hi) && (p_lo > u_lo))
		)
	{
		if (p_lo < b) {
			--p_hi;
		}
		p_lo -= b;
		--q_hat;
	}

	w_lo = (p_lo << 32);
	w_hi = (p_hi << 32) | (p_lo >> 32);

	if (w_lo > r_lo) {
		++w_hi;
	}

	r_lo -= w_lo;
	r_hi -= w_hi;

	q = q_hat << 32;

	/*
	The lower half of the quotient is easier,
	as b is now aligned with r_lo.

	|r_hi][r_lo||    ]
	[ b  ||    ]
	[p_hi||    ][p_lo||    ]
	|
	V
	[q_hi||q_lo]
	*/

	q_hat = ((r_hi << 32) | (r_lo >> 32)) / b_hi;

	p_lo = mul(b, q_hat, &p_hi);

	// r -= b*q_hat
	//
	// ...and at most 2 iterations of this.
	while (
		(p_hi > r_hi) ||
		((p_hi == r_hi) && (p_lo > r_lo))
		)
	{
		if (p_lo < b) {
			--p_hi;
		}
		p_lo -= b;
		--q_hat;
	}

	r_lo -= p_lo;

	q |= q_hat;

	*r = r_lo >> s;

	return q;
}
#endif

#if FIXEDPT_BITS == 64
#define fpabs(x) ((uint64_t)((x >= 0) ? x : -x))
#define fpsign(x) ((x >= 0) ? 1 : -1)
#endif

/* Multiplies two fixedpt numbers, returns the result. */
static inline fixedpt
fixedpt_mul(fixedpt A, fixedpt B)
{
#if FIXEDPT_BITS == 64
	uint64_t y;
	uint64_t res = mul(fpabs(A), fpabs(B), &y);
	return fpsign(A) * fpsign(B) *
		((y << FIXEDPT_WBITS) + (res >> FIXEDPT_FBITS));
#else
	return (((fixedptd)A * (fixedptd)B) >> FIXEDPT_FBITS);
#endif
}

/* Divides two fixedpt numbers, returns the result. */
static inline fixedpt
fixedpt_div(fixedpt A, fixedpt B)
{
#if FIXEDPT_BITS == 64
	uint64_t dummy;
	uint64_t res = div(fpabs(A) << FIXEDPT_FBITS,
		fpabs(A) >> FIXEDPT_WBITS, fpabs(B), &dummy);
	return fpsign(A) * fpsign(B) * res;
#else
	return (((fixedptd)A << FIXEDPT_FBITS) / (fixedptd)B);
#endif
}

/*
 * Note: adding and substracting fixedpt numbers can be done by using
 * the regular integer operators + and -.
 */

#if FIXEDPT_BITS != 64

/**
 * Convert the given fixedpt number to a decimal string.
 * The max_dec argument specifies how many decimal digits to the right
 * of the decimal point to generate. If set to -1, the "default" number
 * of decimal digits will be used (2 for 32-bit fixedpt width, 10 for
 * 64-bit fixedpt width); If set to -2, "all" of the digits will
 * be returned, meaning there will be invalid, bogus digits outside the
 * specified precisions.
 */
static inline void
fixedpt_str(fixedpt A, char *str, int max_dec)
{
	int ndec = 0, slen = 0;
	char tmp[12] = {0};
	fixedptud fr, ip;
	const fixedptud one = (fixedptud)1 << FIXEDPT_BITS;
	const fixedptud mask = one - 1;

	if (max_dec == -1)
#if FIXEDPT_BITS == 32
#if FIXEDPT_WBITS > 16
		max_dec = 2;
#else
		max_dec = 4;
#endif
#elif FIXEDPT_BITS == 64
		max_dec = 10;
#else
#error Invalid width
#endif
	else if (max_dec == -2)
		max_dec = 15;

	if (A < 0) {
		str[slen++] = '-';
		A *= -1;
	}

	ip = fixedpt_toint(A);
	do {
		tmp[ndec++] = '0' + ip % 10;
		ip /= 10;
	} while (ip != 0);

	while (ndec > 0)
		str[slen++] = tmp[--ndec];
	str[slen++] = '.';

	fr = (fixedpt_fracpart(A) << FIXEDPT_WBITS) & mask;
	do {
		fr = (fr & mask) * 10;

		str[slen++] = '0' + (fr >> FIXEDPT_BITS) % 10;
		ndec++;
	} while (fr != 0 && ndec < max_dec);

	if (ndec > 1 && str[slen-1] == '0')
		str[slen-1] = '\0'; /* cut off trailing 0 */
	else
		str[slen] = '\0';
}


/* Converts the given fixedpt number into a string, using a static
 * (non-threadsafe) string buffer */
static inline char*
fixedpt_cstr(const fixedpt A, const int max_dec)
{
	static char str[25];

	fixedpt_str(A, str, max_dec);
	return (str);
}

#endif


/* Returns the square root of the given number, or -1 in case of error */
static inline fixedpt
fixedpt_sqrt(fixedpt A)
{
	int invert = 0;
	int iter = FIXEDPT_FBITS;
	int l, i;

	if (A < 0)
		return (-1);
	if (A == 0 || A == FIXEDPT_ONE)
		return (A);
	if (A < FIXEDPT_ONE && A > 6) {
		invert = 1;
		A = fixedpt_div(FIXEDPT_ONE, A);
	}
	if (A > FIXEDPT_ONE) {
		int s = A;

		iter = 0;
		while (s > 0) {
			s >>= 2;
			iter++;
		}
	}

	/* Newton's iterations */
	l = (A >> 1) + 1;
	for (i = 0; i < iter; i++)
		l = (l + fixedpt_div(A, l)) >> 1;
	if (invert)
		return (fixedpt_div(FIXEDPT_ONE, l));
	return (l);
}


/* Returns the sine of the given fixedpt number. 
 * Note: the loss of precision is extraordinary! */
static inline fixedpt
fixedpt_sin(fixedpt fp)
{
	int sign = 1;
	fixedpt sqr, result;
	const fixedpt SK[2] = {
		fixedpt_rconst(7.61e-03),
		fixedpt_rconst(1.6605e-01)
	};

	fp %= 2 * FIXEDPT_PI;
	if (fp < 0)
		fp = FIXEDPT_PI * 2 + fp;
	if ((fp > FIXEDPT_HALF_PI) && (fp <= FIXEDPT_PI)) 
		fp = FIXEDPT_PI - fp;
	else if ((fp > FIXEDPT_PI) && (fp <= (FIXEDPT_PI + FIXEDPT_HALF_PI))) {
		fp = fp - FIXEDPT_PI;
		sign = -1;
	} else if (fp > (FIXEDPT_PI + FIXEDPT_HALF_PI)) {
		fp = (FIXEDPT_PI << 1) - fp;
		sign = -1;
	}
	sqr = fixedpt_mul(fp, fp);
	result = SK[0];
	result = fixedpt_mul(result, sqr);
	result -= SK[1];
	result = fixedpt_mul(result, sqr);
	result += FIXEDPT_ONE;
	result = fixedpt_mul(result, fp);
	return sign * result;
}


/* Returns the cosine of the given fixedpt number */
static inline fixedpt
fixedpt_cos(fixedpt A)
{
	return (fixedpt_sin(FIXEDPT_HALF_PI - A));
}


/* Returns the tangens of the given fixedpt number */
static inline fixedpt
fixedpt_tan(fixedpt A)
{
	return fixedpt_div(fixedpt_sin(A), fixedpt_cos(A));
}


/* Returns the value exp(x), i.e. e^x of the given fixedpt number. */
static inline fixedpt
fixedpt_exp(fixedpt fp)
{
	fixedpt xabs, k, z, R, xp;
	const fixedpt LN2 = fixedpt_rconst(0.69314718055994530942);
	const fixedpt LN2_INV = fixedpt_rconst(1.4426950408889634074);
	const fixedpt EXP_P[5] = {
		fixedpt_rconst(1.66666666666666019037e-01),
		fixedpt_rconst(-2.77777777770155933842e-03),
		fixedpt_rconst(6.61375632143793436117e-05),
		fixedpt_rconst(-1.65339022054652515390e-06),
		fixedpt_rconst(4.13813679705723846039e-08),
	};

	if (fp == 0)
		return (FIXEDPT_ONE);
	xabs = fixedpt_abs(fp);
	k = fixedpt_mul(xabs, LN2_INV);
	k += FIXEDPT_ONE_HALF;
	k &= ~FIXEDPT_FMASK;
	if (fp < 0)
		k = -k;
	fp -= fixedpt_mul(k, LN2);
	z = fixedpt_mul(fp, fp);
	/* Taylor */
	R = FIXEDPT_TWO +
	    fixedpt_mul(z, EXP_P[0] + fixedpt_mul(z, EXP_P[1] +
	    fixedpt_mul(z, EXP_P[2] + fixedpt_mul(z, EXP_P[3] +
	    fixedpt_mul(z, EXP_P[4])))));
	xp = FIXEDPT_ONE + fixedpt_div(fixedpt_mul(fp, FIXEDPT_TWO), R - fp);
	if (k < 0)
		k = FIXEDPT_ONE >> (-k >> FIXEDPT_FBITS);
	else
		k = FIXEDPT_ONE << (k >> FIXEDPT_FBITS);
	return (fixedpt_mul(k, xp));
}


/* Returns the natural logarithm of the given fixedpt number. */
static inline fixedpt
fixedpt_ln(fixedpt x)
{
	fixedpt log2, xi;
	fixedpt f, s, z, w, R;
	const fixedpt LN2 = fixedpt_rconst(0.69314718055994530942);
	const fixedpt LG[7] = {
		fixedpt_rconst(6.666666666666735130e-01),
		fixedpt_rconst(3.999999999940941908e-01),
		fixedpt_rconst(2.857142874366239149e-01),
		fixedpt_rconst(2.222219843214978396e-01),
		fixedpt_rconst(1.818357216161805012e-01),
		fixedpt_rconst(1.531383769920937332e-01),
		fixedpt_rconst(1.479819860511658591e-01)
	};

	if (x < 0)
		return (0);
	if (x == 0)
		return 0xffffffff;

	log2 = 0;
	xi = x;
	while (xi > FIXEDPT_TWO) {
		xi >>= 1;
		log2++;
	}
	f = xi - FIXEDPT_ONE;
	s = fixedpt_div(f, FIXEDPT_TWO + f);
	z = fixedpt_mul(s, s);
	w = fixedpt_mul(z, z);
	R = fixedpt_mul(w, LG[1] + fixedpt_mul(w, LG[3]
	    + fixedpt_mul(w, LG[5]))) + fixedpt_mul(z, LG[0]
	    + fixedpt_mul(w, LG[2] + fixedpt_mul(w, LG[4]
	    + fixedpt_mul(w, LG[6]))));
	return (fixedpt_mul(LN2, (log2 << FIXEDPT_FBITS)) + f
	    - fixedpt_mul(s, f - R));
}
	

/* Returns the logarithm of the given base of the given fixedpt number */
static inline fixedpt
fixedpt_log(fixedpt x, fixedpt base)
{
	return (fixedpt_div(fixedpt_ln(x), fixedpt_ln(base)));
}


/* Return the power value (n^exp) of the given fixedpt numbers */
static inline fixedpt
fixedpt_pow(fixedpt n, fixedpt exp)
{
	if (exp == 0)
		return (FIXEDPT_ONE);
	if (n < 0)
		return 0;
	return (fixedpt_exp(fixedpt_mul(fixedpt_ln(n), exp)));
}

#endif
