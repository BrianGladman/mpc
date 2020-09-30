/* shared routines for mpcheck-xxx

Copyright (C) 2020 INRIA

This file is part of GNU MPC.

GNU MPC is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

GNU MPC is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
more details.

You should have received a copy of the GNU Lesser General Public License
along with this program. If not, see http://www.gnu.org/licenses/ .
*/

#include <limits.h>
#include <sys/types.h>
#include <sys/resource.h>

unsigned long Max_err_re = 0, Max_err_im = 0;
unsigned long Max_err = 0; /* for real-valued functions like abs */
int recheck = 0;

/* y is the reference value computed by MPC */
static unsigned long
ulp_error (mpfr_t x, mpfr_t y)
{
  mpfr_t z;
  mpfr_prec_t p = mpfr_get_prec (y);
  unsigned long n, special = 0;
  mpfr_exp_t e;

  if (mpfr_cmp (x, y) == 0)
    return 0;

  assert (!mpfr_nan_p (x));
  assert (!mpfr_nan_p (y));
  
  if (mpfr_inf_p (x) || mpfr_inf_p (y))
    {
      if (mpfr_number_p (x) && mpfr_inf_p (y))
        {
          /* only MPFR returns Inf: we replace y by the maximal
             floating-point number */
          if (mpfr_sgn (y) > 0) /* +Inf */
            mpfr_nextbelow (y);
          else
            mpfr_nextabove (y);
          special = 1;
        }
      else 
        {
          /* only glibc returns Inf */
          if (mpfr_sgn (x) > 0) /* +Inf */
            mpfr_nextbelow (x);
          else
            mpfr_nextabove (x);
          special = 2;
        }
    }

  mpfr_init2 (z, p);
  mpfr_sub (z, x, y, MPFR_RNDN);
  mpfr_abs (z, z, MPFR_RNDN);
  /* divide by ulp(y) = 2^(EXP(y) - p) */
  if (!mpfr_zero_p (y))
    {
      /* avoid underflow and take into account subnormals:
         smallest number is 2^(emin-1) */
      if (mpfr_get_exp (y) < emin - 1 + p)
        e = emin - 1;
      else
        e = mpfr_get_exp (y) - p;
    }
  else
    e = emin - 1; /* we consider ulp(0) = 2^(emin-1) */
  mpfr_div_2si (z, z, e, MPFR_RNDN);
  n = mpfr_get_ui (z, MPFR_RNDZ);
  mpfr_clear (z);
  if (special == 1) /* restore original y */
    {
      if (mpfr_sgn (y) > 0)
        mpfr_nextabove (y);
      else
        mpfr_nextbelow (y);
    }
  if (special == 2) /* restore original x */
    {
      if (mpfr_sgn (x) > 0)
        mpfr_nextabove (x);
      else
        mpfr_nextbelow (x);
    }
  if (special == 0)
    return n;
  else
    return (n < ULONG_MAX) ? n + 1 : n;
}

static void
mpcheck_random (mpc_t z)
{
  unsigned long erange;
  
  mpc_urandom (z, state);

  /* set sign and exponent of the real part */
  mpz_set_ui (expz, 2 * (emax - emin + 1));
  mpz_urandomm (expz, state, expz);
  erange = mpz_get_ui (expz);
  /* bit 0 of erange gives the sign */
  if (erange & 1)
    mpfr_neg (mpc_realref (z), mpc_realref (z), MPFR_RNDN);
  mpfr_set_exp (mpc_realref (z), emin + (erange >> 1));
  /* make sure the numbers are representable in IEEE format */
  mpfr_subnormalize (mpc_realref (z), 0, MPFR_RNDN);

  /* set sign and exponent of the imaginary part */
  mpz_set_ui (expz, 2 * (emax - emin + 1));
  mpz_urandomm (expz, state, expz);
  erange = mpz_get_ui (expz);
  /* bit 0 of erange gives the sign */
  if (erange & 1)
    mpfr_neg (mpc_imagref (z), mpc_imagref (z), MPFR_RNDN);
  mpfr_set_exp (mpc_imagref (z), emin + (erange >> 1));
  /* make sure the numbers are representable in IEEE format */
  mpfr_subnormalize (mpc_imagref (z), 0, MPFR_RNDN);
}

/* y is the reference value computed with MPC:
 * return 0 if x and y are equal (including both NaN or Inf)
 * return 1 if x and y disagree (neither NaN nor Inf)
 * return 2 if x is NaN but not y
 * return 4 if x and y are infinities of different signs
 */
static int
mpcheck_mpfr_cmp (mpfr_t x, mpfr_t y)
{
  if (mpfr_number_p (x)  && mpfr_number_p (y))
    return mpfr_cmp (x, y);
  /* now either x or y is NaN or Inf */
  if (mpfr_nan_p (x) || mpfr_nan_p (y))
    {
      if (mpfr_nan_p (x) && mpfr_nan_p (y))
        return 0;
      if (mpfr_nan_p (x))
        return 2;
      fprintf (stderr, "Error, NaN and not NaN\n");
      exit (1);
    }
  /* now either x or y is Inf */
  if (mpfr_inf_p (x) && mpfr_inf_p (y))
    {
      if (mpfr_sgn (x) * mpfr_sgn (y) > 0)
        return 0;
      return 4;
    }
  return 1; /* Inf and non-Inf: differ, will be treated by ulp_error */
}

/* y is the reference value computed with MPC */
static int
mpcheck_mpc_cmp (mpc_t x, mpc_t y)
{
  return mpcheck_mpfr_cmp (mpc_realref (x), mpc_realref (y)) |
         mpcheck_mpfr_cmp (mpc_imagref (x), mpc_imagref (y));
}

static void
report_maximal_errors (void)
{
  printf ("Maximal errors: %lu (real), %lu (imag), %lu (abs/arg)\n",
          Max_err_re, Max_err_im, Max_err);
}

/* return non-zero if both a and b are NaN, or if both are not NaN
   and are equal */
static int
mpfr_agree (mpfr_t a, mpfr_t b, int inex)
{
  if (mpfr_nan_p (a) || mpfr_nan_p (b))
    return mpfr_nan_p (a) && mpfr_nan_p (b);
  else
    /* in the case of even rounding, we can't tell for sure */
    return mpfr_equal_p (a, b) || inex == 2 || inex == -2;
}

#if 0
static int
cputime (void)
{
  struct rusage rus;

  getrusage (0, &rus);
  return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
}
#endif
