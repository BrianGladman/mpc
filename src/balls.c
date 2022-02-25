/* balls -- Functions for complex ball arithmetic.

Copyright (C) 2018, 2020, 2021, 2022 INRIA

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

#include "mpc-impl.h"


void mpcb_out_str (FILE *f, mpcb_srcptr op)
{
   mpc_out_str (f, 10, 0, op->c, MPC_RNDNN);
   fprintf (f, " ");
   mpcr_out_str (f, op->r);
   fprintf (f, "\n");
}


void
mpcb_init (mpcb_ptr rop)
{
   mpc_init2 (rop->c, 2);
   mpcr_set_inf (rop->r);
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


static void
mpcb_set_prec (mpcb_ptr rop, mpfr_prec_t prec)
{
   mpc_set_prec (rop->c, prec);
   mpcr_set_inf (rop->r);
}


void
mpcb_set (mpcb_ptr rop, mpcb_srcptr op)
{
   mpc_set_prec (rop->c, mpc_get_prec (op->c));
   mpc_set (rop->c, op->c, MPC_RNDNN);
   mpcr_set (rop->r, op->r);
}


void
mpcb_set_c (mpcb_ptr rop, mpc_srcptr op, mpfr_prec_t prec,
   unsigned long int err_re, unsigned long int err_im)
   /* Set the precision of rop to prec and assign a ball with centre op
      to it. err_re and err_im contain potential errors in the real and
      imaginary parts of op as multiples of a half ulp. For instance,
      if the real part of op is exact, err_re should be set to 0;
      if it is the result of rounding to nearest, it should be set to 1;
      if it is the result of directed rounding, it should be set to 2.
      The radius of the ball reflects err_re and err_im and the potential
      additional rounding error that can occur when the precision of op
      is higher than prec. If the real part of op is 0, then err_re
      should be 0, since then ulp notation makes no sense, and similarly
      for the imaginary part; otherwise the radius is set to infinity.
      The implementation takes potential different precisions in the real
      and imaginary parts of op into account. */
{
   int inex;
   mpcr_t relerr_re, relerr_im;

   mpc_set_prec (rop->c, prec);
   inex = mpc_set (rop->c, op, MPC_RNDNN);

   if (   (mpfr_zero_p (mpc_realref (op)) && err_re > 0)
       || (mpfr_zero_p (mpc_imagref (op)) && err_im > 0)
       || !mpc_fin_p (op))
       mpcr_set_inf (rop->r);
   else {
      mpcr_set_ui_2si (relerr_re, err_re,
         -mpfr_get_prec (mpc_realref (op)));
         /* prop:relerror of algorithms.tex */
      if (MPC_INEX_RE (inex))
         mpcr_add_rounding_error (relerr_re, prec, MPFR_RNDN);
      mpcr_set_ui_2si (relerr_im, err_im,
         -mpfr_get_prec (mpc_imagref (op)));
      if (MPC_INEX_IM (inex))
         mpcr_add_rounding_error (relerr_im, prec, MPFR_RNDN);
      mpcr_max (rop->r, relerr_re, relerr_im);
         /* prop:comrelerror in algorithms.tex */
   }
}


void
mpcb_neg (mpcb_ptr z, mpcb_srcptr z1)
{
   mpfr_prec_t p;
   int overlap = (z == z1);

   if (!overlap) {
      p = mpcb_get_prec (z1);
      if (mpcb_get_prec (z) != p)
         mpcb_set_prec (z, p);
   }

   mpc_neg (z->c, z1->c, MPC_RNDNN); /* exact */
   mpcr_set (z->r, z1->r);
}


void
mpcb_mul (mpcb_ptr z, mpcb_srcptr z1, mpcb_srcptr z2)
{
   mpcr_t r;
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

   /* generic error of multiplication */
   mpcr_mul (r, z1->r, z2->r);
   mpcr_add (r, r, z1->r);
   mpcr_add (r, r, z2->r);
   /* error of rounding to nearest */
   mpcr_add_rounding_error (r, p, MPFR_RNDN);
   mpcr_set (z->r, r);
}


void
mpcb_sqr (mpcb_ptr z, mpcb_srcptr z1)
{
   mpcr_t r, r2;
   mpfr_prec_t p = mpcb_get_prec (z1);
   int overlap = (z == z1);

   /* Compute the error first in case there is overlap. */
   mpcr_mul_2ui (r2, z1->r, 1);
   mpcr_sqr (r, z1->r);
   mpcr_add (r, r, r2);
   mpcr_add_rounding_error (r, p, MPFR_RNDN);

   if (!overlap)
      mpcb_set_prec (z, p);
   mpc_sqr (z->c, z1->c, MPC_RNDNN);
   mpcr_set (z->r, r);
}


void
mpcb_add (mpcb_ptr z, mpcb_srcptr z1, mpcb_srcptr z2)
{
   mpcr_t r, s, denom;
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
   mpcr_mpc_abs (denom, zc, MPFR_RNDD);
   mpcr_mpc_abs (r, z1->c, MPFR_RNDU);
   mpcr_mul (r, r, z1->r);
   mpcr_mpc_abs (s, z2->c, MPFR_RNDU);
   mpcr_mul (s, s, z2->r);
   mpcr_add (r, r, s);
   mpcr_div (r, r, denom);
   /* error of directed rounding */
   mpcr_add_rounding_error (r, p, MPFR_RNDZ);

   if (overlap)
      mpc_clear (z->c);
   z->c [0] = zc [0];
   mpcr_set (z->r, r);
}


void
mpcb_sqrt (mpcb_ptr z, mpcb_srcptr z1)
{
   mpcr_t r;
   mpfr_prec_t p = mpcb_get_prec (z1);
   int overlap = (z == z1);

   /* Compute the error first in case there is overlap. */
   /* generic error of square root for z1->r <= 0.5:
      0.5*epsilon1 + (sqrt(2)-1) * epsilon1^2
      <= 0.5 * epsilon1 * (1 + epsilon1),
      see eq:propsqrt in algorithms.tex, together with a Taylor
      expansion of 1/sqrt(1-epsilon1) */
   if (!mpcr_lt_half_p (z1->r))
      mpcr_set_inf (r);
   else {
      mpcr_set_one (r);
      mpcr_add (r, r, z1->r);
      mpcr_mul (r, r, z1->r);
      mpcr_div_2ui (r, r, 1);
      /* error of rounding to nearest */
      mpcr_add_rounding_error (r, p, MPFR_RNDN);
   }

   if (!overlap)
      mpcb_set_prec (z, p);
   mpc_sqrt (z->c, z1->c, MPC_RNDNN);
   mpcr_set (z->r, r);
}


void
mpcb_div (mpcb_ptr z, mpcb_srcptr z1, mpcb_srcptr z2)
{
   mpcr_t r, s;
   mpfr_prec_t p = MPC_MIN (mpcb_get_prec (z1), mpcb_get_prec (z2));
   int overlap = (z == z1 || z == z2);
   mpc_t zc;

   if (overlap)
      mpc_init2 (zc, p);
   else {
      zc [0] = z->c [0];
      mpc_set_prec (zc, p);
   }
   mpc_div (zc, z1->c, z2->c, MPC_RNDNN);
   if (overlap)
      mpc_clear (z->c);
   z->c [0] = zc [0];

   /* generic error of division */
   mpcr_add (r, z1->r, z2->r);
   mpcr_set_one (s);
   mpcr_sub_rnd (s, s, z2->r, MPFR_RNDD);
   mpcr_div (r, r, s);
   /* error of rounding to nearest */
   mpcr_add_rounding_error (r, p, MPFR_RNDN);
   mpcr_set (z->r, r);
}


void
mpcb_div_2ui (mpcb_ptr z, mpcb_srcptr z1, unsigned long int e)
{
   mpc_div_2ui (z->c, z1->c, e, MPC_RNDNN);
   mpcr_set (z->r, z1->r);
}


int
mpcb_can_round (mpcb_srcptr op, mpfr_prec_t prec_re, mpfr_prec_t prec_im,
   mpc_rnd_t rnd)
   /* Return a boolean value indicating whether rounding the center of op
      to an mpc_t variable of precision prec_re for the real and prec_im
      for the imaginary part returns a correctly rounded result in
      direction rnd.
      If yes, then using mpcb_round with the same rounding mode sets
      a correct result with the usual MPC semantic. */
{
   mpfr_srcptr re, im;
   mpfr_exp_t exp_re, exp_im, exp_err;

   if (mpcr_inf_p (op->r))
      return 0;
   else if (mpcr_zero_p (op->r))
      return 1;

   re = mpc_realref (op->c);
   im = mpc_imagref (op->c);
   /* The question makes sense only if neither the real nor the imaginary
      part of the centre are 0: Otherwise we need to have an absolute error
      that is less than the smallest representable number; since rounding
      only once at precision p introduces an error of about 2^-p, this
      means that the precision needs to be about as big as the negative
      of the minimal exponent, which is astronomically large. */
   if (mpfr_zero_p (re) || mpfr_zero_p (im))
      return 0;

   exp_re = mpfr_get_exp (re);
   exp_im = mpfr_get_exp (im);

   /* Absolute error of the real part, as given in the proof of
      prop:comrelerror of algorithms.tex:
      |x-x~|  = |x~*theta_R - y~*theta_I|
             <= |x~ - y~| * epsilon, where epsilon is the complex
                                     relative error
             <= (|x~|+|y~|) * epsilon
             <= 2 * max (|x~|, |y~|) * epsilon
      To call mpfr_can_round, we only need the exponent in base 2,
      which is then bounded above by
                1 + max (exp_re, exp_im) + exponent (epsilon) */
   exp_err = 1 + MPC_MAX (exp_re, exp_im) + mpcr_get_exp (op->r);

   return (   mpfr_can_round (re, exp_re - exp_err, MPFR_RNDN,
                 MPC_RND_RE (rnd), prec_re)
           && mpfr_can_round (im, exp_im - exp_err, MPFR_RNDN,
                 MPC_RND_IM (rnd), prec_im));
}


void
mpcb_round (mpc_ptr rop, mpcb_srcptr op, mpc_rnd_t rnd)
   /* Set rop to the centre of op. To make sure that this corresponds
      to the MPC semantics of returning a correctly rounded result, one
      needs to call mpcb_can_round first. */
{
   mpc_set (rop, op->c, rnd);
}

