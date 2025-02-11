#include <R.h>
#include <Rinternals.h>
#include <zlib.h>
#include <string.h>
#include <stdlib.h>

// This function is to parse the names.tsv.gz files returned by salmon
// These are single-line files, but with tab-separated values making
// up that single line.
// In reality, there will be around 250000 values in the single line

SEXP parse_trans_names(SEXP r_filename) {
    const char *filename = CHAR(STRING_ELT(r_filename, 0));

    // Open gzipped file
    gzFile file = gzopen(filename, "rb");
    if (file == NULL) {
        Rf_error("Failed to open file: %s", filename);
    }

    // Read the file into a dynamically allocated buffer
    size_t buffer_size = 1024;
    size_t total_size = 0;
    char *buffer = (char *)malloc(buffer_size);
    if (buffer == NULL) {
        gzclose(file);
        Rf_error("Memory allocation failed");
    }

    int bytes_read;
    while ((bytes_read = gzread(file, buffer + total_size, buffer_size - total_size)) > 0) {
        total_size += bytes_read;
        if (total_size == buffer_size) {
            buffer_size *= 2;
            buffer = (char *)realloc(buffer, buffer_size);
            if (buffer == NULL) {
                gzclose(file);
                Rf_error("Memory reallocation failed");
            }
        }
    }
    gzclose(file);

    if (bytes_read < 0) {
        free(buffer);
        Rf_error("Error reading gzipped file");
    }

    // Null-terminate the buffer
    buffer[total_size] = '\0';

    // **Remove trailing newline if present**
    if (total_size > 0 && buffer[total_size - 1] == '\n') {
        buffer[total_size - 1] = '\0';
        total_size--;  // Adjust size to exclude the newline
    }

    // Split the buffer by tab characters
    size_t value_count = 0;
    char *token = strtok(buffer, "\t");
    char **values = NULL;

    while (token != NULL) {
        values = (char **)realloc(values, (value_count + 1) * sizeof(char *));
        if (values == NULL) {
            free(buffer);
            Rf_error("Memory reallocation failed");
        }
        values[value_count] = strdup(token);  // Duplicate the token string
        if (values[value_count] == NULL) {
            free(buffer);
            Rf_error("Memory allocation failed");
        }
        value_count++;
        token = strtok(NULL, "\t");
    }

    // Create an R character vector to store the result
    SEXP result = PROTECT(allocVector(STRSXP, value_count));
    for (size_t i = 0; i < value_count; i++) {
        SET_STRING_ELT(result, i, mkChar(values[i]));
        free(values[i]);
    }

    free(values);
    free(buffer);

    UNPROTECT(1);
    return result;
}
