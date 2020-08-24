/* balls -- Functions for complex ball arithmetic.

Copyright (C) 2018, 2020 INRIA

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
#include "mpc-impl.h"

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
/* FIXME: For the time being, we assume that z is different from z1 and from z2 */
{
   double r;
   mpfr_prec_t p = mpc_get_prec (z1->c);

   mpcb_set_prec (z, p);
   mpc_mul (z->c, z1->c, z2->c, MPC_RNDNN);

   fesetround (FE_UPWARD);
   /* generic error of multiplication */
   r = z1->r + z2->r + z1->r * z2->r;
   /* error of rounding to nearest */
   r += ldexp (1 + r, -p);
   z->r = r;
}


void
mpcb_add (mpcb_ptr z, mpcb_srcptr z1, mpcb_srcptr z2)
/* FIXME: For the time being, we assume that z is different from z1 and from z2 */
{
   double r, denom, x, y;
   mpfr_prec_t p = mpc_get_prec (z1->c);

   mpcb_set_prec (z, p);
   mpc_add (z->c, z1->c, z2->c, MPC_RNDZZ);
      /* rounding towards 0 makes the generic error easier to compute,
         but incurs a tiny penalty for the rounding error */

   /* generic error of addition:
      r <= (|z1|*r1 + |z2|*r2) / |z1+z2|
        <= (|z1|*r1 + |z2|*r2) / |z| since we rounded towards 0 */
   fesetround (FE_TOWARDZERO);
   x = mpfr_get_d (mpc_realref (z->c), MPFR_RNDZ);
   y = mpfr_get_d (mpc_imagref (z->c), MPFR_RNDZ);
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
   z->r = r;
}

void
mpcb_sqrt (mpcb_ptr z, mpcb_srcptr z1)
/* FIXME: For the time being, we assume that z is different from z1 */
{
   double r;
   mpfr_prec_t p = mpc_get_prec (z1->c);

   mpcb_set_prec (z, p);
   mpc_sqrt (z->c, z1->c, MPC_RNDNN);

   fesetround (FE_UPWARD);
   /* generic error of square root for z->r <= 0.5:
      0.5*epsilon1 + (sqrt(2)-1) * epsilon1^2
      see eq:propsqrt in algorithms.tex, together with a Taylor
      expansion of 1/sqrt(1-epsilon1) */
   r = ldexp (z1->r, -1) + 0.415 * z1->r * z1->r;
   /* error of rounding to nearest */
   r += ldexp (1 + r, -p);
   z->r = r;
}

void
mpcb_div_2ui (mpcb_ptr z, mpcb_srcptr z1, unsigned long int e)
{
   mpc_div_2ui (z->c, z1->c, e, MPC_RNDNN);
   z->r = z1->r;
}

