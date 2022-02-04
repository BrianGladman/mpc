/* tballs -- test file for complex ball arithmetic.

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

#include "mpc-tests.h"

#if 0
static int
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

   mpcb_out_str (stdout, z);
   for (i = 1; i <= n; i++) {
      mpcb_add (tmp, z, minusone);
      printf (" +");
      mpcb_print (tmp);
      mpc_set_ui_ui (c, i, 0, MPC_RNDNN);
      mpcb_set_c (factor, c);
      mpcb_mul (z, tmp, factor);
      printf ("%3i ", i);
      mpcb_out_str (stdout, z);
   }

   mpfr_clear (one);
   mpc_clear (c);
   mpcb_clear (z);
   mpcb_clear (minusone);
   mpcb_clear (factor);
   mpcb_clear (tmp);

   return -1;
}
#endif


static int
test_agm (void)
{
   mpfr_prec_t p, target;
   mpc_t c;
   mpcb_t a, b, a1;
   mpc_t agma, agmb;
   mpc_rnd_t rnd = MPC_RNDDU;
   int i, n;

   p = 1000;
   target = p - 9;
   n = 20;

   mpc_init2 (c, p);
   mpc_set_si_si (c, 100, 0, MPC_RNDNN);
   mpcb_init_set_c (a, c);
   mpc_set_si_si (c, 0, 100, MPC_RNDNN);
   mpcb_init_set_c (b, c);
   mpcb_init (a1);
   mpc_init2 (agma, target);
   mpc_init2 (agmb, target);

   printf (" 0 ");
   mpcb_out_str (stdout, a);
   printf ("  ");
   mpcb_out_str (stdout, b);
   for (i = 1; i <= n; i++) {
      mpcb_add (a1, a, b);
      mpcb_mul (b, a, b);
      mpcb_div_2ui (a, a1, 1);
      mpcb_sqrt (b, b);
      printf ("%2i ", i);
      mpcb_out_str (stdout, a);
      printf ("   ");
      mpcb_out_str (stdout, b);
      if (   mpcb_can_round (a, target, target, rnd)
          && mpcb_can_round (b, target, target, rnd)) {
         mpcb_round (agma, a, rnd);
         mpcb_round (agmb, b, rnd);
         if (mpc_cmp (agma, agmb) == 0)
            return 0;
      }
   }

   mpc_clear (c);
   mpcb_clear (a);
   mpcb_clear (b);
   mpcb_clear (a1);
   mpc_clear (agma);
   mpc_clear (agmb);

   return 1;
}


int
main (void)
{
   return test_agm ();
}

