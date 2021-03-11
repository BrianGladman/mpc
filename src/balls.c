/* balls -- Functions for complex ball arithmetic.

Copyright (C) 2018, 2020, 2021 INRIA

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

#include <stdio.h>
#include <math.h>
#include <fenv.h>
#include <assert.h>
#include "mpc-impl.h"

#define FE_ERROR (FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW)
#define FE_CLEARERROR feclearexcept (FE_ERROR);
#define FE_TESTERROR assert (!fetestexcept(FE_ERROR));

void mpcb_print (mpcb_srcptr op)
{
   mpc_out_str (stdout, 10, 0, op->c, MPC_RNDNN);
   printf (" %.20g\n", op->r);
}


void
mpcb_init (mpcb_ptr rop)
{
   mpc_init2 (rop->c, 2);
   rop->r = INFINITY;
}


void
mpcb_clear (mpcb_ptr rop)
{
   mpc_clear (rop->c);
}


mpfr_prec_t
mpcb_get_prec (mpcb_srcptr op)
{
   return mpc_get_prec (op->c);
}


void
mpcb_set_prec (mpcb_ptr rop, mpfr_prec_t prec)
{
   mpc_set_prec (rop->c, prec);
   rop->r = INFINITY;
}


void
mpcb_set_c (mpcb_ptr rop, mpc_srcptr op)
{
   mpc_set_prec (rop->c, MPC_MAX_PREC (op));
   mpc_set (rop->c, op, MPC_RNDNN);
   rop->r = 0.0;
}


void
mpcb_set (mpcb_ptr rop, mpcb_srcptr op)
{
   mpc_set_prec (rop->c, mpc_get_prec (op->c));
   mpc_set (rop->c, op->c, MPC_RNDNN);
   rop->r = op->r;
}


void
mpcb_init_set_c (mpcb_ptr rop, mpc_srcptr op)
{
   mpc_init2 (rop->c, MPC_MAX_PREC (op));
   mpc_set (rop->c, op, MPC_RNDNN);
   rop->r = 0.0;
}


void
mpcb_mul (mpcb_ptr z, mpcb_srcptr z1, mpcb_srcptr z2)
{
   double r;
   mpfr_prec_t p = MPC_MIN (mpcb_get_prec (z1), mpcb_get_prec (z2));
   int overlap = (z == z1 || z == z2);
   mpc_t zc;

   if (overlap)
      mpc_init2 (zc, p);
   else {
      zc [0] = z->c [0];
      mpc_set_prec (zc, p);
   }
   mpc_mul (zc, z1->c, z2->c, MPC_RNDNN);
   if (overlap)
      mpc_clear (z->c);
   z->c [0] = zc [0];

   FE_CLEARERROR
   fesetround (FE_UPWARD);
   /* generic error of multiplication */
   r = z1->r + z2->r + z1->r * z2->r;
   /* error of rounding to nearest */
   r += ldexp (1 + r, -p);
   z->r = r;
   FE_TESTERROR
}


void
mpcb_add (mpcb_ptr z, mpcb_srcptr z1, mpcb_srcptr z2)
{
   double r, denom, x, y;
   mpfr_prec_t p = MPC_MIN (mpcb_get_prec (z1), mpcb_get_prec (z2));
   int overlap = (z == z1 || z == z2);
   mpc_t zc;

   if (overlap)
      mpc_init2 (zc, p);
   else {
      zc [0] = z->c [0];
      mpc_set_prec (zc, p);
   }
   mpc_add (zc, z1->c, z2->c, MPC_RNDZZ);
      /* rounding towards 0 makes the generic error easier to compute,
         but incurs a tiny penalty for the rounding error */

   /* generic error of addition:
      r <= (|z1|*r1 + |z2|*r2) / |z1+z2|
        <= (|z1|*r1 + |z2|*r2) / |z| since we rounded towards 0 */
   FE_CLEARERROR
   fesetround (FE_TOWARDZERO);
   x = mpfr_get_d (mpc_realref (zc), MPFR_RNDZ);
   y = mpfr_get_d (mpc_imagref (zc), MPFR_RNDZ);
   denom = sqrt (x*x + y*y);
   fesetround (FE_UPWARD);
   x = mpfr_get_d (mpc_realref (z1->c), MPFR_RNDA);
   y = mpfr_get_d (mpc_imagref (z1->c), MPFR_RNDA);
   r = sqrt (x*x + y*y) * z1->r;
   x = mpfr_get_d (mpc_realref (z2->c), MPFR_RNDA);
   y = mpfr_get_d (mpc_imagref (z2->c), MPFR_RNDA);
   r += sqrt (x*x + y*y) * z2->r;
   r /= denom;
   /* error of directed rounding */
   r += ldexp (1 + r, 1-p);
   FE_TESTERROR

   if (overlap)
      mpc_clear (z->c);
   z->c [0] = zc [0];
   z->r = r;
}


void
mpcb_sqrt (mpcb_ptr z, mpcb_srcptr z1)
{
   double r;
   mpfr_prec_t p = mpcb_get_prec (z1);
   int overlap = (z == z1);

   /* Compute the error first in case there is overlap. */
   FE_CLEARERROR
   fesetround (FE_UPWARD);
   /* generic error of square root for z1->r <= 0.5:
      0.5*epsilon1 + (sqrt(2)-1) * epsilon1^2
      see eq:propsqrt in algorithms.tex, together with a Taylor
      expansion of 1/sqrt(1-epsilon1) */
   assert (z1->r <= 0.5);
   r = ldexp (z1->r, -1) + 0.415 * z1->r * z1->r;
   /* error of rounding to nearest */
   r += ldexp (1 + r, -p);
   FE_TESTERROR

   if (!overlap)
      mpcb_set_prec (z, p);
   mpc_sqrt (z->c, z1->c, MPC_RNDNN);
   z->r = r;
}


void
mpcb_div_2ui (mpcb_ptr z, mpcb_srcptr z1, unsigned long int e)
{
   mpc_div_2ui (z->c, z1->c, e, MPC_RNDNN);
   z->r = z1->r;
}


int
mpcb_can_round (mpcb_srcptr op, mpfr_prec_t prec_re, mpfr_prec_t prec_im)
{
   mpfr_srcptr re, im;
   mpfr_exp_t exp_re, exp_im, exp_err;
   int exp_int;
   double err;

   re = mpc_realref (op->c);
   im = mpc_imagref (op->c);
   /* The question makes sense only if neither the real nor the imaginary
      part of the centre are 0: Otherwise, we can either round this part
      to 0 or we cannot round. But to round to 0, we need to have an
      absolute error that is less than the smallest representable number;
      otherwise said, the precision needs to be about as big as the negative
      of the minimal exponent, which is astronomically large. */
   if (mpfr_zero_p (re) || mpfr_zero_p (im))
      return 0;

   exp_re = mpfr_get_exp (re);
   exp_im = mpfr_get_exp (im);

   fesetround (FE_UPWARD);
   /* Bound on the relative error of the real part. */
   err = (1 + ldexp (1.0, exp_im - exp_re + 1)) * op->r;
   /* Absolute error. */
   err *= fabs (mpfr_get_d (re, MPFR_RNDA));
   /* Exponent of the error as a power of 2, rounded up. */
   frexp (err, &exp_int);
   exp_err = exp_int;
   if (!mpfr_can_round (re, exp_re - exp_err, MPFR_RNDN, MPFR_RNDN, prec_re))
      return 0;

   err = (1 + ldexp (1.0, exp_re - exp_im + 1)) * op->r;
   err *= fabs (mpfr_get_d (im, MPFR_RNDA));
   frexp (err, &exp_int);
   exp_err = exp_int;
   return mpfr_can_round (im, exp_im - exp_err, MPFR_RNDN, MPFR_RNDN, prec_im);
}


void
mpcb_round (mpc_ptr rop, mpcb_srcptr op)
{
   mpfr_set (mpc_realref (rop), mpc_realref (op->c), MPFR_RNDN);
   mpfr_set (mpc_imagref (rop), mpc_imagref (op->c), MPFR_RNDN);
}

