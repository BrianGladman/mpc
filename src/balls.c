/* balls -- Functions for complex ball arithmetic.

Copyright (C) 2018 INRIA

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
   printf (" %g\n", op->r);
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
   double r, r1, r2;
   mpfr_prec_t p = mpc_get_prec (z1->c);

   mpcb_set_prec (z, p);
   mpc_add (z->c, z1->c, z2->c, MPC_RNDNN);

   fesetround (FE_UPWARD);
   /* absolute error of z1 and z2, divided by sqrt(2) */
   r1 = ldexp (z1->r, MPC_MAX (mpfr_get_exp (mpc_realref (z1->c)),
                               mpfr_get_exp (mpc_imagref (z1->c))));
   r2 = ldexp (z2->r, MPC_MAX (mpfr_get_exp (mpc_realref (z2->c)),
                               mpfr_get_exp (mpc_imagref (z2->c))));
   /* absolute error of z, divided by sqrt(2) */
   r = r1 + r2;
   /* relative error of z; the sqrt(2), which come from the approximation
      of complex norms by powers of 2, cancel out */
   r = ldexp (r, - (MPC_MIN (mpfr_get_exp (mpc_realref (z->c)),
                             mpfr_get_exp (mpc_imagref (z->c)))
                    - 1));
      /* The -1 comes from taking the lower bound of 1/2 for the mantissa
         of c, whereas above we used the upper bound of 1 for r1 and r2.
         So we may lose a factor of 2 here, for instance when z1==z2.
         To prevent this, we would need to work with an approximation
         of z, z1 and z2 in the type of r. */
   /* error of rounding to nearest */
   r += ldexp (1 + r, -p);
   z->r = r;
}

void
mpcb_sqrt (mpcb_ptr z, mpcb_srcptr z1)
/* FIXME: Actually compute the error instead of assuming it is the same... */
/* FIXME: For the time being, we assume that z is different from z1 */
{
   mpc_sqrt (z->c, z1->c, MPC_RNDNN);
   z->r = z1->r;
}

void
mpcb_div_2ui (mpcb_ptr z, mpcb_srcptr z1, unsigned long int e)
{
   mpc_div_2ui (z->c, z1->c, e, MPC_RNDNN);
   z->r = z1->r;
}

