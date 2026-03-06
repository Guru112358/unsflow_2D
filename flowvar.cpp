#include "flowvar.h"

void compute_primitives(primitive &var,conserved &cons,material &mat);




void compute_primitives(primitive &var,conserved &cons,material &mat)
{

    var.rho = cons.U[0];
    var.u   = cons.U[1] / var.rho;
    var.v   = cons.U[2] / var.rho;
    var.p   = (mat.gamma - 1.0) * (cons.U[3] - 0.5 * var.rho * (var.u*var.u + var.v*var.v));

}
