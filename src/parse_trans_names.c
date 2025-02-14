#include <R.h>
#include <Rinternals.h>
#include <zlib.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

// This function is to parse the names.tsv.gz files returned by salmon
// These are single-line files, but with tab-separated values making
// up that single line.
// In reality, there will be around 250000 values in the single line

SEXP parse_trans_names(SEXP r_filename) {
    // Input validation - quick and early
    if (!isString(r_filename) || length(r_filename) != 1) {
        Rf_error("Filename must be a single character string");
    }

    const char *filename = CHAR(STRING_ELT(r_filename, 0));
    if (!filename || strlen(filename) == 0) {
        Rf_error("Empty filename provided");
    }

    gzFile file = gzopen(filename, "rb");
    if (file == NULL) {
        Rf_error("Failed to open file: %s", filename);
    }

    // Start with a larger buffer since we expect ~250k values
    // Assuming average string length of 20 chars + tab, initial 5MB should
    // cover many cases without reallocation
    size_t buffer_size = 5 * 1024 * 1024;
    size_t total_size = 0;
    size_t tab_count = 0;
    char *buffer = (char *)malloc(buffer_size);
    if (buffer == NULL) {
        gzclose(file);
        Rf_error("Memory allocation failed");
    }

    // Combined read and tab counting loop
    int bytes_read;
    while ((bytes_read = gzread(file, buffer + total_size, buffer_size - total_size)) > 0) {
        // Count tabs in this chunk
        for (int i = 0; i < bytes_read; i++) {
            if (buffer[total_size + i] == '\t') tab_count++;
        }

        total_size += bytes_read;

        // Only resize if really necessary
        if (total_size == buffer_size) {
            // Check for overflow before multiplying
            if (buffer_size > SIZE_MAX / 2) {
                free(buffer);
                gzclose(file);
                Rf_error("File too large, would cause buffer size overflow");
            }
            size_t new_size = buffer_size * 2;
            char *new_buffer = (char *)realloc(buffer, new_size);
            if (new_buffer == NULL) {
                free(buffer);
                gzclose(file);
                Rf_error("Memory reallocation failed");
            }
            buffer = new_buffer;
            buffer_size = new_size;
        }
    }

    // Check for gzip-specific errors
    int gz_err;
    const char *gz_err_msg = gzerror(file, &gz_err);
    gzclose(file);
    if (gz_err != Z_OK && gz_err != Z_STREAM_END) {
        free(buffer);
        Rf_error("Gzip error: %s", gz_err_msg);
    }

    if (total_size == 0) {
        free(buffer);
        Rf_error("Empty file");
    }

    // Null-terminate with bounds checking
    if (total_size >= buffer_size) {
        char *new_buffer = (char *)realloc(buffer, total_size + 1);
        if (new_buffer == NULL) {
            free(buffer);
            Rf_error("Memory reallocation failed");
        }
        buffer = new_buffer;
    }
    buffer[total_size] = '\0';

    // Handle trailing whitespace
    while (total_size > 0 && (buffer[total_size - 1] == '\n' ||
           buffer[total_size - 1] == '\r' ||
           buffer[total_size - 1] == '\t')) {
        buffer[--total_size] = '\0';
    }

    // Number of values is tab count + 1
    size_t value_count = tab_count + 1;

    // Pre-allocate array for all values
    char **values = (char **)malloc(value_count * sizeof(char *));
    if (values == NULL) {
        free(buffer);
        Rf_error("Memory allocation failed");
    }

    // Split the string - using a pointer-based approach instead of strtok
    // for better performance
    size_t current_value = 0;
    char *start = buffer;
    char *end = buffer;

    while (*end != '\0' && current_value < value_count) {
        // Find next tab or end of string
        while (*end != '\t' && *end != '\0') end++;

        // Temporarily null-terminate this segment
        char tmp = *end;
        *end = '\0';

        // Duplicate the string segment
        values[current_value] = strdup(start);
        if (values[current_value] == NULL) {
            for (size_t i = 0; i < current_value; i++) {
                free(values[i]);
            }
            free(values);
            free(buffer);
            Rf_error("Memory allocation failed");
        }

        // Restore the original character and move to next segment
        *end = tmp;
        current_value++;

        if (*end == '\t') {
            start = end + 1;
            end = start;
        }
    }

    // Create R vector
    SEXP result = PROTECT(allocVector(STRSXP, value_count));
    if (result == NULL) {
        for (size_t i = 0; i < value_count; i++) {
            free(values[i]);
        }
        free(values);
        free(buffer);
        UNPROTECT(1);
        Rf_error("Failed to allocate R vector");
    }

    // Convert to R strings
    for (size_t i = 0; i < value_count; i++) {
        SEXP str = mkChar(values[i]);
        if (str == NULL) {
            for (size_t j = 0; j < value_count; j++) {
                free(values[j]);
            }
            free(values);
            free(buffer);
            UNPROTECT(1);
            Rf_error("Failed to create R character element");
        }
        SET_STRING_ELT(result, i, str);
        free(values[i]);
    }

    free(values);
    free(buffer);

    UNPROTECT(1);
    return result;
}
