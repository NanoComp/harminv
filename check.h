#ifndef CHECK_H
#define CHECK_H

/* some useful error-checking macros: */

#define CHECK(condition, message) do { \
     if (!(condition))  { \
          fprintf(stderr, "CHECK failure on line %d of " __FILE__ ": " \
		  message "\n", __LINE__);  exit(EXIT_FAILURE); \
     } \
} while (0)

#define CHK_MALLOC(p, t, n) do {                              \
     size_t CHK_MALLOC_n_tmp = (n);                           \
     (p) = (t *) malloc(sizeof(t) * CHK_MALLOC_n_tmp);        \
     CHECK((p) || CHK_MALLOC_n_tmp == 0, "out of memory!");   \
} while (0)

#endif /* CHECK_H */
