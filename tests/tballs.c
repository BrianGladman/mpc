/* tballs -- test file for complex ball arithmetic.

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

#include "mpc-tests.h"

int
test_exp (void)
{
   mpfr_prec_t p;
   mpc_t c;
   mpfr_t one;
   mpcb_t z, minusone, factor, tmp;
   int i, n;

   p = 20;
   n = 4;
   mpfr_init2 (one, p);
   mpc_init2 (c, p);
   mpc_set_si_si (c, -1, -1, MPC_RNDNN);
   mpcb_init_set_c (minusone, c);
   mpfr_set_ui (one, 1, MPFR_RNDN);
   mpfr_exp (mpc_realref (c), one, MPFR_RNDU);
   mpfr_exp (mpc_imagref (c), one, MPFR_RNDD);
   mpcb_init_set_c (z, c);
   mpcb_init (tmp);
   mpcb_init (factor);

   mpcb_print (z);
   for (i = 1; i <= n; i++) {
      mpcb_add (tmp, z, minusone);
      printf (" +");
      mpcb_print (tmp);
      mpc_set_ui_ui (c, i, 0, MPC_RNDNN);
      mpcb_set_c (factor, c);
      mpcb_mul (z, tmp, factor);
      printf ("%3i ", i);
      mpcb_print (z);
   }

   mpfr_clear (one);
   mpc_clear (c);
   mpcb_clear (z);
   mpcb_clear (minusone);
   mpcb_clear (factor);
   mpcb_clear (tmp);

   return -1;
}

int
test_agm (void)
{
   mpfr_prec_t p;
   mpc_t c;
   mpcb_t a, b, a1, b1;
   int i, n;

   p = 20;
   n = 10;
   mpc_init2 (c, p);
   mpc_set_si_si (c, 1, 0, MPC_RNDNN);
   mpcb_init_set_c (a, c);
   mpc_set_si_si (c, 0, 1, MPC_RNDNN);
   mpcb_init_set_c (b, c);
   mpcb_init (a1);
   mpcb_init (b1);

   printf ("0 ");
   mpcb_print (a);
   printf ("  ");
   mpcb_print (b);
   for (i = 1; i <= n; i++) {
      mpcb_add (a1, a, b);
      mpcb_mul (b1, a, b);
      mpcb_div_2ui (a, a1, 1);
      mpcb_sqrt (b, b1);
      printf ("%i ", i);
      mpcb_print (a);
      printf ("  ");
      mpcb_print (b);
   }

   mpc_clear (c);
   mpcb_clear (a);
   mpcb_clear (b);
   mpcb_clear (a1);
   mpcb_clear (b1);

   return -1;
}

int
main (void)
{
   return test_agm ();
}

