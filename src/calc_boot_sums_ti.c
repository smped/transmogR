#include <R.h>
#include <Rinternals.h>
#include <zlib.h>
#include <float.h>
#include <errno.h>
#include <stdlib.h>
#include <stdint.h>

SEXP calc_boot_row_vals(SEXP filename_sexp, SEXP n_trans_sexp, SEXP n_boot_sexp) {
    // Protect R objects from garbage collection
    PROTECT_INDEX px;
    PROTECT_WITH_INDEX(filename_sexp, &px);

    // Convert inputs from SEXP to C types
    const char* filename = CHAR(STRING_ELT(filename_sexp, 0));
    int n_trans = INTEGER(n_trans_sexp)[0];
    int n_boot = INTEGER(n_boot_sexp)[0];

    // Input validation
    if (n_trans <= 1) {
        UNPROTECT(1);
        Rf_error("Invalid number of transcripts: n_trans=%d", n_trans);
    }

    if (n_boot <= 1) {
        UNPROTECT(1);
        Rf_error("Invalid number of bootstraps:  n_boot=%d", n_boot);
    }

    // Calculate total values needed and check for overflow
    double total_values_double = (double)n_trans * (double)n_boot;
    if (total_values_double > (double)SIZE_MAX / sizeof(double)) {
        UNPROTECT(1);
        Rf_error("Memory size overflow");
    }
    size_t total_values = (size_t)total_values_double;
    size_t memory_needed = total_values * sizeof(double);

    // Create R vector for results
    SEXP result_sexp = PROTECT(allocVector(REALSXP, n_trans));
    double *result = REAL(result_sexp);
    memset(result, 0, n_trans * sizeof(double));

    // Open the gzipped file
    errno = 0;
    gzFile file = gzopen(filename, "rb");
    if (file == NULL) {
        UNPROTECT(2);
        Rf_error("Failed to open file '%s': %s", filename, strerror(errno));
    }

    // Check file size if possible
    z_off_t initial_pos = gztell(file);
    if (gzseek(file, 0L, SEEK_END) != -1) {
        z_off_t file_size = gztell(file);
        gzseek(file, initial_pos, SEEK_SET);

        if (file_size > memory_needed) {
            warning("File '%s' contains more data than expected (%zu bytes vs %zu bytes needed). Extra data will be ignored.",
                    filename, (size_t)file_size, memory_needed);
        }
    } else {
        gzrewind(file);
    }

    // Allocate memory for matrix using R's allocator
    SEXP matrix_sexp = PROTECT(allocVector(REALSXP, total_values));
    double *matrix = REAL(matrix_sexp);

    // Read data in chunks
    const size_t chunk_size = 1024 * 1024;  // 1MB chunks
    size_t remaining = memory_needed;
    size_t offset = 0;
    char *buffer = (char *)matrix;

    while (remaining > 0) {
        size_t to_read = (remaining < chunk_size) ? remaining : chunk_size;
        int bytes_read = gzread(file, buffer + offset, to_read);

        if (bytes_read < 0) {
            int errnum;
            const char *errmsg = gzerror(file, &errnum);
            gzclose(file);
            UNPROTECT(3);
            Rf_error("Error reading from file: %s", errmsg);
        }

        if (bytes_read == 0) break;  // EOF

        remaining -= bytes_read;
        offset += bytes_read;
    }

    if (remaining > 0) {
        gzclose(file);
        UNPROTECT(3);
        Rf_error("Incomplete read: expected %zu bytes, got %zu", memory_needed, offset);
    }

    // Check for extra data
    char extra_byte;
    if (gzread(file, &extra_byte, 1) > 0) {
        warning("Additional data exists in file after expected matrix data");
    }

    gzclose(file);

    for (int i = 0; i < n_trans; i++) {
        double row_sum = 0.0;
        double row_mean;

        for (int j = 0; j < n_boot; j++) {
            row_sum += matrix[i + j * n_trans];
        }
        row_mean = row_sum / n_boot;

        if (fabs(row_mean) < DBL_EPSILON) {
            result[i] = 0.0;
            continue;
        }

        double sum_squared_diffs = 0.0;
        for (int j = 0; j < n_boot; j++) {
            double val = matrix[i + j * n_trans];
            double diff = val - row_mean;
            sum_squared_diffs += diff * diff;
        }

        result[i] = sum_squared_diffs / row_mean;
    }

    UNPROTECT(3);
    return result_sexp;
}
