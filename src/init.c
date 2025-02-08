#include <R.h>
#include <Rinternals.h>
#include <Rversion.h>

extern void calc_boot_row_vals(char **filename, int *n_trans, int *n_boot, double *result);  // Declare your C function

static const R_CallMethodDef CallEntries[] = {
    {"calc_boot_row_vals", (DL_FUNC) &calc_boot_row_vals, 4},
    {NULL, NULL, 0}  // End of the list
};

void R_init_mypackage(DllInfo *dll) {
    /* Register routines */
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);

    /* Disable symbol search */
    R_useDynamicSymbols(dll, FALSE);
}
