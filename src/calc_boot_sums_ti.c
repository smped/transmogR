#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <R.h>  // Include R header for Rprintf

// After parsing the bootstraps.gz file, this calculates the first values in
// the sums described by Baldoni et al, which are then able to be summed across
// samples to create the final, moderated overdispersions.

void calc_boot_row_vals(char **filename, int *n_trans, int *n_boot, double *result) {

    // Open the gzipped file
    gzFile file = gzopen(*filename, "rb");
    if (file == NULL) {
        perror("Error opening gzipped file");
        return;
    }

    int total_values = (*n_trans) * (*n_boot);
    // Allocate memory to store matrix data
    double *matrix = (double *)malloc(total_values * sizeof(double));
    if (matrix == NULL) {
        perror("Memory allocation failed");
        gzclose(file);
        return;
    }

    // Parse the data
    size_t read_count = gzread(file, matrix, total_values * sizeof(double));
    if (read_count != total_values * sizeof(double)) {
        Rprintf("Expected %d values but read %zu bytes\n", total_values, read_count);
        free(matrix);
        gzclose(file);
        return;
    }

    gzclose(file); // Close the file

    // Compute row statistics as requested
    for (int i = 0; i < *n_trans; i++) {
        double row_sum = 0.0;

        // Calculate row mean
        for (int j = 0; j < *n_boot; j++) {
            row_sum += matrix[i + j * (*n_trans)];  // Fill by column (column-major order)
        }
        double row_mean = row_sum / (*n_boot);

        if (row_mean == 0.0) {
            result[i] = 0.0;  // Handle potential division by zero
            continue;
        }

        // Compute sum of squared differences divided by row mean
        double sum_squared_diffs = 0.0;
        for (int j = 0; j < *n_boot; j++) {
            double diff = matrix[i + j * (*n_trans)] - row_mean;  // Subtract row mean
            sum_squared_diffs += diff * diff;
        }
        result[i] = sum_squared_diffs / row_mean;  // Return result for this row
    }

    free(matrix);
}
