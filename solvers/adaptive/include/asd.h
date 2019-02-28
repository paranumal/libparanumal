#ifndef ASD_H
#define ASD_H

#include <mpi.h>
#include <stdint.h>
#include <stdio.h>

#define ldouble long double

// {{{ compiler helper macros
// Define gcc version macro  __GNUC_PREREQ if not defined
#if defined __GNUC__ && !defined __GNUC_PREREQ
#ifndef __GNUC_MINOR__
#define __GNUC_PREREQ(maj, min) 0
#else
#define __GNUC_PREREQ(maj, min)                                                \
  ((__GNUC__ << 16) + __GNUC_MINOR__ >= ((maj) << 16) + (min))
#endif
#endif

#ifndef ASD_NORETURN
#if defined(__clang__)
#if __has_feature(attribute_analyzer_noreturn)
#define ASD_NORETURN __attribute__((analyzer_noreturn))
#else
#define ASD_NORETURN
#endif
#else
#define ASD_NORETURN
#endif
#endif
// }}}

// {{{ misc macros
#define ASD_BUFSIZ 8192

#define ASD_XSTR(a) ASD_STR(a)
#define ASD_STR(a) #a

#define ASD_APPEND(x, y) x##y
#define ASD_APPEND_EXPAND(x, y) ASD_APPEND(x, y)

#define ASD_MIN(a, b) (((a) < (b)) ? (a) : (b))
#define ASD_MAX(a, b) (((a) > (b)) ? (a) : (b))

#define ASD_SWAP_PTR(x, y)                                                     \
  do                                                                           \
  {                                                                            \
    void *swap_temp = y;                                                       \
    y = x;                                                                     \
    x = swap_temp;                                                             \
  } while (0)
// }}}

// {{{ error checking macros
#define ASD_SYS_ERROR_CHECK(cond, msg)                                         \
  do                                                                           \
  {                                                                            \
    if (cond)                                                                  \
    {                                                                          \
      perror(msg);                                                             \
      ASD_ABORT("perror");                                                     \
    }                                                                          \
  } while (0)

#define ASD_MPI_CHECK(c) ASD_ABORT_IF_NOT((c) == MPI_SUCCESS, "MPI Error")
// }}}

// {{{ abort and assert
#define ASD_NOOP()                                                             \
  do                                                                           \
  {                                                                            \
  } while (0)
#define ASD_ABORT(...) asd_abort_verbose(__FILE__, __LINE__, __VA_ARGS__)
#define ASD_ABORT_IF(q, ...) ((q) ? ASD_ABORT(__VA_ARGS__) : (void)0)
#define ASD_ABORT_IF_NOT(q, ...) ASD_ABORT_IF(!(q), __VA_ARGS__)
#ifdef ASD_DEBUG
#define ASD_ASSERT(expression)                                                 \
  ASD_ABORT_IF_NOT((expression), "Assert Failed: '" #expression "'")
#else
#define ASD_ASSERT(expression) ASD_NOOP()
#endif

/** Verbose abort function.
 *
 * Typically this function is not called directly but the \c ASD_ABORT,
 * \c ASD_ABORT_IF, and \c ASD_ABORT_IF_NOT macros are used.
 *
 * \param[in] file file name the abort occurred in (can be obtained from
 *                 \c __FILE__)
 * \param[in] line line number the abort occurred at (can be obtained from
 *                 \c __LINE__)
 * \param[in] note note to print before aborting
 *
 * \return This function called \c asd_abort() and will not return.
 */
void asd_abort_verbose(const char *file, int line, ...) ASD_NORETURN;
// }}}

// {{{ logging
#define ASD_LC_ALL 0
#define ASD_LC_ROOT 1

#define ASD_LL_DEFAULT -1
#define ASD_LL_ALWAYS 0
#define ASD_LL_TRACE 1
#define ASD_LL_DEBUG 2
#define ASD_LL_VERBOSE 3
#define ASD_LL_INFO 4
#define ASD_LL_WARNING 5
#define ASD_LL_ERROR 6
#define ASD_LL_SILENT 7

/*
 * Setting a hard log threshold is set at compile time
 */
#ifdef ASD_LOG_LEVEL
#define ASD_LL_THRESHOLD ASD_LOG_LEVEL
#else
#ifdef ASD_DEBUG
#define ASD_LL_THRESHOLD ASD_LL_TRACE
#else
#define ASD_LL_THRESHOLD ASD_LL_VERBOSE
#endif
#endif

#define ASD_ROOT_TRACE(...) ASD_LOG(ASD_LC_ROOT, ASD_LL_TRACE, __VA_ARGS__)
#define ASD_ROOT_LDEBUG(...) ASD_LOG(ASD_LC_ROOT, ASD_LL_DEBUG, __VA_ARGS__)
#define ASD_ROOT_VERBOSE(...) ASD_LOG(ASD_LC_ROOT, ASD_LL_VERBOSE, __VA_ARGS__)
#define ASD_ROOT_INFO(...) ASD_LOG(ASD_LC_ROOT, ASD_LL_INFO, __VA_ARGS__)
#define ASD_ROOT_WARNING(...) ASD_LOG(ASD_LC_ROOT, ASD_LL_WARNING, __VA_ARGS__)
#define ASD_ROOT_LERROR(...) ASD_LOG(ASD_LC_ROOT, ASD_LL_ERROR, __VA_ARGS__)

#define ASD_TRACE(...) ASD_LOG(ASD_LC_ALL, ASD_LL_TRACE, __VA_ARGS__)
#define ASD_LDEBUG(...) ASD_LOG(ASD_LC_ALL, ASD_LL_DEBUG, __VA_ARGS__)
#define ASD_VERBOSE(...) ASD_LOG(ASD_LC_ALL, ASD_LL_VERBOSE, __VA_ARGS__)
#define ASD_INFO(...) ASD_LOG(ASD_LC_ALL, ASD_LL_INFO, __VA_ARGS__)
#define ASD_WARNING(...) ASD_LOG(ASD_LC_ALL, ASD_LL_WARNING, __VA_ARGS__)
#define ASD_LERROR(...) ASD_LOG(ASD_LC_ALL, ASD_LL_ERROR, __VA_ARGS__)

#define ASD_LOG(category, level, ...)                                          \
  ((level) < ASD_LL_THRESHOLD                                                  \
       ? (void)0                                                               \
       : asd_log_printf(__FILE__, __LINE__, (category), (level), __VA_ARGS__))

/** Initialization function for the logging system.
 *
 * \param[in] rank      the MPI rank of the current process
 * \param[in] stream    the stream to output the logging information
 * \param[in] threshold the threshold for the logging system (use
 *                      \c ASD_LL_DEFAULT for the default value or
 *                      \c ASD_LL_ALWAYS to print all log messages)
 *
 */
void asd_log_init(int rank, FILE *stream, int threshold);

/** Logging function called by the logging macros.
 *
 * This function is typically not called directly in application code.  Instead
 * one of the logging macros (such as \c ASD_INFO or \c ASD_ROOT_INFO is
 * called) will call this function.
 *
 * Messages with \a category \c ASD_LC_ROOT are only printed by the root rank
 * and \a category \c ASD_LC_ALL are printed by all ranks.
 *
 * The parameter \a level is used to determine if the message will be printed
 * or not based on the threshold passed into \c asd_log_init.
 *
 * \param[in] file     file name the abort occurred in (can be obtained from
 *                     \c __FILE__)
 * \param[in] line     line number the abort occurred at (can be obtained from
 *                     \c __LINE__)
 * \param[in] category logging category (e.g., \c ASD_LC_ROOT or \c
 *                     ASD_LC_ALL)
 * \param[in] level    logging level (e.g., \c ASD_LL_TRACE, \c ASD_LL_ERROR,
 *                     ...)
 * \param[in] ...      It is assumed that the first of the variable arguments
 *                     is a \c printf style format and its required arguments
 *                     follow.
 */
void asd_log_printf(const char *file, int line, int category, int level, ...);
// }}}

// {{{ memory allocation and alignment
/** \c malloc wrapper.
 *
 * This wrapper calls \c malloc and aborts if there is a memory error.
 * The returned pointer needs to be freed by \c asd_free();
 *
 * \param[in] size allocation size
 *
 * \return pointer to allocated memory.
 */
void *asd_malloc(size_t size);

/** \c calloc wrapper.
 *
 * This wrapper calls \c calloc and aborts if there is a memory error.
 * The returned pointer needs to be freed by \c asd_free();
 *
 * \param[in] nmemb number of elements
 * \param[in] size  size of each element
 *
 * \return pointer to allocated memory.
 */
void *asd_calloc(size_t nmemb, size_t size);

/** \c realloc wrapper.
 *
 * This wrapper calls \c realloc and aborts if there is a memory error.
 * The returned pointer needs to be freed by \c asd_free();
 *
 * \param[in] ptr  pointer to memory to reallocate
 * \param[in] size allocation size
 *
 * \return pointer to reallocated memory.
 */
void *asd_realloc(void *ptr, size_t size);

/** \c free wrapper
 *
 * This function frees memory.
 *
 * \param[in,out] ptr pointer to memory to free.
 */
void asd_free(void *ptr);

#define ASD_IS_ALIGNED(p, a) (((intptr_t)(p) & ((a)-1)) == 0)

/*
 * Switch for different compilers taken from this web page:
 *
 *     http://nadeausoftware.com/articles/2012/10/c_c_tip_how_detect_compiler_name_and_version_using_compiler_predefined_macros
 *
 */
#if defined(__clang__)
/* Clang/LLVM. ---------------------------------------------- */
#define ASD_ASSUME_ALIGNED(lvalueptr, align) ASD_NOOP()
#elif defined(__ICC) || defined(__INTEL_COMPILER)
/* Intel ICC/ICPC. ------------------------------------------ */
#define ASD_ASSUME_ALIGNED(lvalueptr, align)                                   \
  __assume_aligned(lvalueptr, align);                                          \
  ASD_ASSERT(ASD_IS_ALIGNED(lvalueptr, align))

#elif defined(__GNUC__) || defined(__GNUG__)
/* GNU GCC/G++. --------------------------------------------- */
#if __GNUC_PREREQ(4, 7)
//      If  gcc_version >= 4.7
#define ASD_ASSUME_ALIGNED(lvalueptr, align)                                   \
  lvalueptr = __builtin_assume_aligned(lvalueptr, align);                      \
  ASD_ASSERT(ASD_IS_ALIGNED(lvalueptr, align))
#else
//       Else
#define ASD_ASSUME_ALIGNED(lvalueptr, align) ASD_NOOP()
#endif

#elif defined(__HP_cc) || defined(__HP_aCC)
/* Hewlett-Packard C/aC++. ---------------------------------- */
#define ASD_ASSUME_ALIGNED(lvalueptr, align) ASD_NOOP()

#elif defined(__IBMC__) || defined(__IBMCPP__)
/* IBM XL C/C++. -------------------------------------------- */
#define ASD_ASSUME_ALIGNED(lvalueptr, align) ASD_NOOP()

#elif defined(_MSC_VER)
/* Microsoft Visual Studio. --------------------------------- */
#define ASD_ASSUME_ALIGNED(lvalueptr, align) ASD_NOOP()

#elif defined(__PGI)
/* Portland Group PGCC/PGCPP. ------------------------------- */
#define ASD_ASSUME_ALIGNED(lvalueptr, align) ASD_NOOP()

#elif defined(__SUNPRO_C) || defined(__SUNPRO_CC)
/* Oracle Solaris Studio. ----------------------------------- */
#define ASD_ASSUME_ALIGNED(lvalueptr, align) ASD_NOOP()

#endif

#if (defined __GNUC__) || (defined __PGI) || (defined __IBMC__)
#define ASD_ALIGN(n) __attribute__((aligned(n)))
#elif (defined _MSC_VER)
#define ASD_ALIGN(n) __declspec(align(n))
#else
#error Need equilvent of __attribute__((aligned(n))) for this compiler
#endif

/** \c malloc wrapper for cache line aligned memory.
 *
 * This wrapper aligns memory to cache lines.  One needs to call \c
 * asd_free_aligned() to free the allocated memory.
 *
 * \param[in] size allocation size
 *
 * \return pointer to cache line aligned allocated memory.
 */
void *asd_malloc_aligned(size_t size);

/** \c free wrapper for cache line aligned memory.
 *
 * This function frees memory that was allocated with \c asd_malloc_aligned().
 *
 * \param[in,out] ptr pointer to cache line aligned memory to free.
 */
void asd_free_aligned(void *ptr);
// }}}

// {{{ signals
/** Set a signal handler which prints stack traces on terminating signals.
 */
void asd_signal_handler_set();
// }}}

// {{{ critbit
typedef struct
{
  void *root;
} asd_critbit0_tree_t;

/** Membership testing.
 *
 * The following function takes a tree, \a t, and a \c NULL terminated string,
 * \a u,  and returns non-zero iff \a u in \a t.
 *
 * \param [in] t tree
 * \param [in] u possible member
 * \returns non-zero iff \a u in \a t
 */
int asd_critbit0_contains(asd_critbit0_tree_t *t, const char *u);

/** Inserting into the tree.
 *
 * It takes a tree, \a t, and possibly mutates it such that a \c NULL
 * terminated string, \a u, is a member on exit.
 *
 * \param [in,out] t tree
 * \param [in] u possible member
 * \returns:
 *   $\cases{ 0 &if {\rm out of memory} \cr
 *            1 &if {\it u} {\rm was already a member} \cr
 *            2 &if {\it t} {\rm was mutated successfully}}$.
 */
int asd_critbit0_insert(asd_critbit0_tree_t *t, const char *u);

/** Deleting elements.
 *
 * This function takes a tree, \a t, and a \c NULL terminated string,
 * \a u, and possibly mutates the tree such that $u \notin t$.
 *
 * \param [in,out] t tree
 * \param [in] u possible member to remove from \a t
 * \returns It returns 1 if the tree was mutated, 0 otherwise.
 */
int asd_critbit0_delete(asd_critbit0_tree_t *t, const char *u);

/** Clearing a tree.
 *
 * Clearing a tree (freeing all members) brings us our first code for walking
 * the whole tree rather than just tracing a path through it.
 *
 * So, the \c critbit0_clear function takes a tree, \a t, and frees every
 * member of it, mutating the tree such that it is empty on exit.
 *
 * \param [in,out] t tree
 */
void asd_critbit0_clear(asd_critbit0_tree_t *t);

/** Fetching elements with a given prefix.
 *
 * One of the operations which crit-bit trees can perform efficiently that hash
 * tables cannot is the extraction of the subset of elements with a given
 * prefix.
 *
 * The following function takes a tree, \a t, and a \c NULL terminated string,
 * \a prefix. Let $S \subseteq t$ where $x \in S$ iff \a prefix is a prefix of
 * \c x, then $\forall x : S.$ \a handle is called with arguments \c x and
 * \c arg.
 * \returns:
 *   $\cases{ 0 &if {\it handle} {\rm returned 0} \cr
 *            1 &if {\rm successful} \cr
 *            2 &if {\it handle} {\rm returned a value} $\notin [0,1]$}$
 * \note (Note that, if |handle| returns 0, the iteration is aborted)
 */
int asd_critbit0_allprefixed(asd_critbit0_tree_t *t, const char *prefix,
                             int (*handle)(const char *, void *), void *arg);
// }}}

// {{{ dictionary
typedef struct
{
  size_t num_entries;
  asd_critbit0_tree_t t;
} asd_dictionary_t;

/** Dictionary to initialize.
 *
 * The following function takes a pointer to a dictionary, \a d, and initializes
 * it.
 *
 * \param [out] d pointer to dictionary
 */
void asd_dictionary_init(asd_dictionary_t *d);

/** Membership testing.
 *
 * The following function takes a dictionary, \a d, and a \c NULL terminated
 * string, \a u,  and returns non-zero iff \a u in \a d.
 *
 * \param [in] d dictionary
 * \param [in] u possible member
 * \returns non-zero iff \a u in \a d
 */
int asd_dictionary_contains(asd_dictionary_t *d, const char *u);

/** Inserting key and value pair into a dictionary.
 *
 * It takes a dictionary, \a d, and possibly mutates it such that a \c NULL
 * terminated strings, \a key and \a value, is a member on exit.
 *
 * \param [in,out] d dictionary
 * \param [in] key possible key
 * \param [in] val possible value
 * \returns:
 *   $\cases{ 0 &if {\rm out of memory} \cr
 *            1 &if {\it key} {\rm was already a member} \cr
 *            2 &if {\it d} {\rm was mutated successfully}}$.
 */
int asd_dictionary_insert(asd_dictionary_t *d, const char *key,
                          const char *val);

/** Inserting key and value pair where value is a pointer into a dictionary.
 *
 * It takes a dictionary, \a d, and possibly mutates it such that a \c NULL
 * terminated strings, \a key and \a value, is a member on exit.
 *
 * \param [in,out] d dictionary
 * \param [in] key possible key
 * \param [in] val possible pointer value
 * \returns:
 *   $\cases{ 0 &if {\rm out of memory} \cr
 *            1 &if {\it key} {\rm was already a member} \cr
 *            2 &if {\it d} {\rm was mutated successfully}}$.
 */
int asd_dictionary_insert_ptr(asd_dictionary_t *d, const char *key,
                              const void *val);

/** Inserting key and value pair where value is a \c int into a
 * dictionary.
 *
 * It takes a dictionary, \a d, and possibly mutates it such that a \c NULL
 * terminated strings, \a key and \a value \c snprintf'd, is a member on exit.
 *
 * \param [in,out] d dictionary
 * \param [in] key possible key
 * \param [in] val possible \c int
 * \returns:
 *   $\cases{ 0 &if {\rm out of memory} \cr
 *            1 &if {\it key} {\rm was already a member} \cr
 *            2 &if {\it d} {\rm was mutated successfully}}$.
 */
int asd_dictionary_insert_int(asd_dictionary_t *d, const char *key,
                              const int val);

/** Inserting key and value pair where value is a \c int32_t into a
 * dictionary.
 *
 * It takes a dictionary, \a d, and possibly mutates it such that a \c NULL
 * terminated strings, \a key and \a value \c snprintf'd, is a member on exit.
 *
 * \param [in,out] d dictionary
 * \param [in] key possible key
 * \param [in] val possible \c int32_t
 * \returns:
 *   $\cases{ 0 &if {\rm out of memory} \cr
 *            1 &if {\it key} {\rm was already a member} \cr
 *            2 &if {\it d} {\rm was mutated successfully}}$.
 */

int asd_dictionary_insert_int32(asd_dictionary_t *d, const char *key,
                                const int32_t val);

/** Inserting key and value pair where value is a \c size_t into a
 * dictionary.
 *
 * It takes a dictionary, \a d, and possibly mutates it such that a \c NULL
 * terminated strings, \a key and \a value \c snprintf'd, is a member on exit.
 *
 * \param [in,out] d dictionary
 * \param [in] key possible key
 * \param [in] val possible \c size_t
 * \returns:
 *   $\cases{ 0 &if {\rm out of memory} \cr
 *            1 &if {\it key} {\rm was already a member} \cr
 *            2 &if {\it d} {\rm was mutated successfully}}$.
 */

int asd_dictionary_insert_size_t(asd_dictionary_t *d, const char *key,
                                 const size_t val);

/** Return a value given a key.
 *
 * It takes a dictionary, \a d, returns a pointer to the value associated with
 * a \c NULL terminated \a key.
 *
 * \param [in] d dictionary
 * \param [in] key possible key
 * \returns:
 *   $\cases{ \c NULL & if {\it key} {\rm is not a member} \cr
 *          {\rm pointer to value} & if {\it key} {\rm is a member}}$
 */
char *asd_dictionary_get_value(asd_dictionary_t *d, const char *key);

/** Return a value given a key assuming value is a pointer.
 *
 * It takes a dictionary, \a d, returns a pointer that points to where the value
 * pointed associated with a \c NULL terminated \a key.
 *
 * \param [in] d dictionary
 * \param [in] key possible key
 * \returns:
 *   $\cases{ \c NULL & if {\it key} {\rm is not a member} \cr
 *          {\rm pointer to value} & if {\it key} {\rm is a member}}$
 */
void *asd_dictionary_get_value_ptr(asd_dictionary_t *d, const char *key);

/** Return a value given a key assuming value is a \c int.
 *
 * It takes a dictionary, \a d, returns an \c int associated with a \c NULL
 * terminated \a key.
 *
 * \param [in]  d dictionary
 * \param [in]  key possible key
 * \param [out] val value if function returned \c 1
 *
 * \returns:
 *   $\cases{ \c 0 & if {\it key} {\rm is not a member} \cr
 *            \c 1 & if {\it key} {\rm is a member}}$
 */
int asd_dictionary_get_value_int(asd_dictionary_t *d, const char *key,
                                 int *val);

/** Return a value given a key assuming value is a \c int32_t.
 *
 * It takes a dictionary, \a d, returns an \c int32_t associated with a \c NULL
 * terminated \a key.
 *
 * \param [in]  d dictionary
 * \param [in]  key possible key
 * \param [out] val value if function returned \c 1
 *
 * \returns:
 *   $\cases{ \c 0 & if {\it key} {\rm is not a member} \cr
 *            \c 1 & if {\it key} {\rm is a member}}$
 */
int asd_dictionary_get_value_int32(asd_dictionary_t *d, const char *key,
                                   int32_t *val);

/** Delete key and value pair into a dictionary.
 *
 * It takes a dictionary, \a d, and possibly mutates it such that a \c NULL
 * terminated string, \a key  is not member on exit.
 *
 * \param [in,out] d dictionary
 * \param [in] key possible key
 * \returns:
 *   $\cases{ 0 &if {\rm key not found} \cr
 *            1 &if {\it key deleted}}$.
 */
int asd_dictionary_delete(asd_dictionary_t *d, const char *key);

/** Clearing a dictionary.
 *
 * Clearing a dictionary (freeing all members) brings us our first code for
 * walking the whole dictionary rather than just tracing a path through it.
 *
 * So, the \c asd_dictionary_clear function takes a dictionary, \a d, and frees
 * every member of it, mutating the dictionary such that it is empty on exit.
 *
 * \param [in,out] d dictionary
 */
void asd_dictionary_clear(asd_dictionary_t *d);

/** Fetching values with a given prefix.
 *
 * The following function takes a dictionary, \a d, and a \c NULL terminated
 * string, \a prefix. Let $S \subseteq d$ where $x \in S$ iff \a prefix is a
 * prefix of \c x, then $\forall x : S.$ \a handle is called with arguments \c x
 * and \c arg.
 * \returns:
 *   $\cases{ 0 &if {\it handle} {\rm returned 0} \cr
 *            1 &if {\rm successful} \cr
 *            2 &if {\it handle} {\rm returned a value} $\notin [0,1]$}$
 * \note (Note that, if |handle| returns 0, the iteration is aborted)
 */
int asd_dictionary_allprefixed(asd_dictionary_t *t, const char *prefix,
                               int (*handle)(const char *, const char *,
                                             void *),
                               void *arg);

/** Fetching pointer values with a given prefix.
 *
 * The following function takes a dictionary, \a d, and a \c NULL terminated
 * string, \a prefix. Let $S \subseteq d$ where $x \in S$ iff \a prefix is a
 * prefix of \c x, then $\forall x : S.$ \a handle is called with arguments \c x
 * and \c arg.
 * \returns:
 *   $\cases{ 0 &if {\it handle} {\rm returned 0} \cr
 *            1 &if {\rm successful} \cr
 *            2 &if {\it handle} {\rm returned a value} $\notin [0,1]$}$
 * \note (Note that, if |handle| returns 0, the iteration is aborted)
 *
 * \note The void * input to the handle is the pointer stored in the value, not
 * the pointer to the pointer
 */
int asd_dictionary_allprefixed_ptr(asd_dictionary_t *t, const char *prefix,
                                   int (*handle)(const char *, void *, void *),
                                   void *arg);

/** Fetching int values with a given prefix.
 *
 * The following function takes a dictionary, \a d, and a \c NULL terminated
 * string, \a prefix. Let $S \subseteq d$ where $x \in S$ iff \a prefix is a
 * prefix of \c x, then $\forall x : S.$ \a handle is called with arguments \c x
 * and \c arg.
 * \returns:
 *   $\cases{ 0 &if {\it handle} {\rm returned 0} \cr
 *            1 &if {\rm successful} \cr
 *            2 &if {\it handle} {\rm returned a value} $\notin [0,1]$}$
 * \note (Note that, if |handle| returns 0, the iteration is aborted)
 */
int asd_dictionary_allprefixed_int(asd_dictionary_t *t, const char *prefix,
                                   int (*handle)(const char *, int, void *),
                                   void *arg);

/** Fetching size_t values with a given prefix.
 *
 * The following function takes a dictionary, \a d, and a \c NULL terminated
 * string, \a prefix. Let $S \subseteq d$ where $x \in S$ iff \a prefix is a
 * prefix of \c x, then $\forall x : S.$ \a handle is called with arguments \c x
 * and \c arg.
 * \returns:
 *   $\cases{ 0 &if {\it handle} {\rm returned 0} \cr
 *            1 &if {\rm successful} \cr
 *            2 &if {\it handle} {\rm returned a value} $\notin [0,1]$}$
 * \note (Note that, if |handle| returns 0, the iteration is aborted)
 */
int asd_dictionary_allprefixed_size_t(asd_dictionary_t *t, const char *prefix,
                                      int (*handle)(const char *, size_t,
                                                    void *),
                                      void *arg);

// }}}

// {{{ random number (pcg32)
// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
// Copied from: http://www.pcg-random.org/
// Copied from: https://github.com/imneme/pcg-c-basic

typedef struct
{
  uint64_t state;
  uint64_t inc;
} asd_pcg32_random_t;

void asd_pcg32_srandom_r(asd_pcg32_random_t *rng, uint64_t initstate,
                         uint64_t initseq);

uint32_t asd_pcg32_boundedrand_r(asd_pcg32_random_t *rng, uint32_t bound);
// }}}

// {{{ file
/** Read the contents of a file into a string.
 */
char *asd_read_file(const char *filename, size_t *len);
// }}}

// {{{ endian
#define ASD_LITTLE_ENDIAN 0
#define ASD_BIG_ENDIAN 1

/** Determine endian of integers
 *
 * \returns 0 if int is little endian and 1 if big endian
 */
int asd_endian();
// }}}

// {{{ MPI
/** Given a MPI Comm this returns the host rank.
 *
 * The host rank is the rank of the process in the list of processes that
 * have the same MPI processor name.
 */
int asd_get_host_rank(MPI_Comm comm);
// }}}

// {{{ Lua
#ifdef ASD_USE_LUA

#include <lauxlib.h>
#include <lua.h>
#include <lualib.h>

/*
 * Helper function for calling a lua function. Based on generic call function of
 * Listing 25.4-25.6 of
 * @book{Ierusalimschy2006Lua,
 *  author = {Ierusalimschy, Roberto},
 *  title = {Programming in Lua, Second Edition},
 *  year = {2006},
 *  isbn = {8590379825},
 *  publisher = {Lua.Org},
 * }
 */
int asd_lua_global_function_call(lua_State *L, const char *name,
                                 const char *sig, ...);

/*
 * Some of the following functions were modified from the code at:
 *
 *    http://windrealm.org/tutorials/reading-a-lua-configuration-file-from-c.php
 */

/** Evaluates a Lua expression and returns the boolean result.
 *
 * If an error occurs or the result is not a boolean, def is returned.
 */
int asd_lua_expr_boolean(lua_State *L, const char *expr, int def);

/** Evaluates a Lua expression and returns the integer result.
 *
 * If an error occurs or the result is not a integer, def is returned.
 */
lua_Integer asd_lua_expr_integer(lua_State *L, const char *expr,
                                 lua_Integer def);

/** Evaluates a Lua expression and returns the number result.
 *
 * If an error occurs or the result is not a number, def is returned.
 */
lua_Number asd_lua_expr_number(lua_State *L, const char *expr, lua_Number def);

/** Evaluates a Lua expression and returns the string result.
 *
 * If an error occurs or the result is not a string, a copy of def is returned.
 * The user of this function is responsible for cleaning up the memory.
 *
 */
char *asd_lua_expr_string(lua_State *L, const char *expr, const char *def);

#endif
// }}}

// {{{ Jacobi
/** Compute the (\c N + 1) node Jacobi--Gauss--Lobatto quadrature.
 *
 * This quadrature integrates $N$th order polynomials exactly with the weight
 * function $(1-x)^\alpha (1+x)^\beta$.
 *
 * \param[in]  alpha     $\alpha$ parameter from the weight function
 * \param[in]  beta      $\beta$ parameter from the weight function
 * \param[in]  N         number of quadrature nodes minus one
 * \param[out] x         location of the quadrature nodes in $[-1,1]$
 * \param[out] w         quadrature weights
 *
 */
void asd_jacobi_gauss_lobatto_quadrature(ldouble alpha, ldouble beta,
                                         int N, ldouble *x,
                                         ldouble *w);

/** Compute the (\c N + 1) column Jacobi Vandermonde matrix.
 *
 * The $n$th column of this matrix is the $n$th order Jacobi polynomial
 * $P^{(\alpha,\beta)}_n$ evaluated at $x$.
 *
 * \param[in]  alpha     $\alpha$ parameter from the Jacobi polynomials
 * \param[in]  beta      $\beta$ parameter from the Jacobi polynomials
 * \param[in]  N         The highest order Jacobi polynomial to be computed
 * \param[in]  nx        number of points in \c x
 * \param[in]  x         location of the evaluation of the Jacobi polynomials
 * \param[out] V         the Jacobi Vandermonde matrix
 *
 */
void asd_jacobi_p_vandermonde(ldouble alpha, ldouble beta, int N,
                              size_t nx, ldouble *x, ldouble *V);

/** Compute the $N$th order interpolation matrix.
 *
 * This builds the $N$th order nodal interpolation matrix associated with the
 * Jacobi Vandermonde matrix $V$ where the interpolant is evaluated at $x$.
 *
 * \param[in]  alpha     $\alpha$ parameter from the Jacobi polynomials
 * \param[in]  beta      $\beta$ parameter from the Jacobi polynomials
 * \param[in]  N         The order of \c V and \c D
 * \param[in]  nx        number of points in \c x
 * \param[in]  x         location of the evaluation of the interpolant
 * \param[in]  V         the Jacobi Vandermonde matrix; this should be the $N$
 *                       order Vandermonde matrix evaluated at the
 *                       interpolation nodes
 * \param[out] I         the nodal interpolation matrix
 *
 */
void asd_jacobi_p_interpolation(ldouble alpha, ldouble beta, int N,
                                size_t nx, ldouble *x, ldouble *V,
                                ldouble *I);

/** Compute the $N$th order differentiation matrix.
 *
 * This builds the $N$th order nodal differentiation matrix associated with the
 * Jacobi Vandermonde matrix $V$ where the derivative is evaluated at $x$.
 *
 * \param[in]  alpha     $\alpha$ parameter from the Jacobi polynomials
 * \param[in]  beta      $\beta$ parameter from the Jacobi polynomials
 * \param[in]  N         The order of \c V and \c D
 * \param[in]  nx        number of points in \c x
 * \param[in]  x         location of the evaluation of the Derivative
 * \param[in]  V         the Jacobi Vandermonde matrix; this should be the $N$
 *                       order Vandermonde matrix evaluated at the
 *                       interpolation nodes
 * \param[out] D         the nodal differentiation matrix
 *
 */
void asd_jacobi_p_differentiation(ldouble alpha, ldouble beta, int N,
                                  size_t nx, ldouble *x, ldouble *V,
                                  ldouble *D);

/** Compute the $N$th order mass matrix.
 *
 * This builds the $N$th order nodal mass matrix associated with the
 * Jacobi Vandermonde matrix $V$ evaluated at $x$.
 *
 * \param[in]  N         The order of \c V and \c M
 * \param[in]  V         the Jacobi Vandermonde matrix; this should be the $N$
 *                       order Vandermonde matrix evaluated at the
 *                       interpolation nodes
 * \param[out] M         the nodal mass matrix
 *
 */
void asd_jacobi_p_mass(int N, ldouble *V, ldouble *M);

/** Compute the $N$th order h L2 projection matrix.
 *
 * P = h V V^T I^T M
 *
 * \param[in]  N         The order of \c V and \c D
 * \param[in]  h         scalar projection ratio
 * \param[in]  V         the Jacobi Vandermonde matrix; this should be the $N$
 *                       order Vandermonde matrix evaluated at the
 *                       interpolation nodes
 * \param[in]  I         the nodal interpolation matrix
 * \param[in]  M         the nodal mass matrix
 * \param[out] P         the nodal L2 projection matrix
 *
 */
void asd_jacobi_p_h_project(int N, ldouble h, ldouble *V,
                            ldouble *I, ldouble *M, ldouble *P);
// }}}


#endif
