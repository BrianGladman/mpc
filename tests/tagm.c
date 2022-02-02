/* targ -- test file for mpc_agm.

Copyright (C) 2022 INRIA

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

static int
test_agm (void)
{
   mpc_t a, b, z;
   mpfr_prec_t prec = 1000;

   mpc_init2 (a, prec);
   mpc_init2 (b, prec);
   mpc_init2 (z, prec);

   /* Example where a is real and b is close to -a and which causes
      cancellation in the first step. */
   mpc_set_ui_ui (a, 1, 0, MPC_RNDNN);
   mpc_mul_2ui (a, a, 100, MPC_RNDNN);
   mpc_neg (b, a, MPC_RNDNN);
   mpfr_add_ui (mpc_realref (b), mpc_realref (b), 1, MPFR_RNDZ);
   mpfr_set_ui (mpc_imagref (b), 1, MPFR_RNDN);
   mpc_agm (z, a, b, MPC_RNDNN);

   /* Example where a is real and b is close to a and which causes a loss
      in the error analysis when switching from a relative error to an
      absolute error (in ulp) for the imaginary part. */
   mpc_set (b, a, MPC_RNDNN);
   mpfr_set_ui (mpc_imagref (b), 1, MPFR_RNDN);
   mpc_agm (z, a, b, MPC_RNDNN);

   /* Example where a is real and b is close to -a and which causes a loss
      in the error analysis when switching from a relative error to an
      absolute error (in ulp) for the real part. */
   mpc_neg (b, a, MPC_RNDNN);
   mpfr_set_ui (mpc_imagref (b), 1, MPFR_RNDN);
   mpc_agm (z, a, b, MPC_RNDNN);

   /* Non-trivial example with an angle of pi. */
   mpc_set_ui_ui (a, 1, 2, MPC_RNDNN);
   mpc_set_si_si (b, -10, -20, MPC_RNDNN);
   mpc_agm (z, a, b, MPC_RNDNN);

   /* Non-trivial example with an angle 0. */
   mpc_set_ui_ui (a, 1, 2, MPC_RNDNN);
   mpc_set_ui_ui (b, 10, 20, MPC_RNDNN);
   mpc_agm (z, a, b, MPC_RNDNN);

   return 0;
}


int
main (void)
{
  return test_agm ();
}

