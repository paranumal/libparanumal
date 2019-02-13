#include "asd.h"

#include <err.h>
#include <execinfo.h>
#include <fcntl.h>
#include <float.h>
#include <inttypes.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sysexits.h>
#include <unistd.h>

// {{{ abort
static void asd_abort() { abort(); }

void asd_abort_verbose(const char *file, int line, ...)
{
  va_list ap;
  const char *fmt;
  char note[ASD_BUFSIZ];

  va_start(ap, line);
  fmt = va_arg(ap, const char *);
  vsnprintf(note, ASD_BUFSIZ, fmt, ap);
  va_end(ap);
  ASD_LERROR("Abort: [%s:%d] %s", file, line, note);
  asd_abort();
}
// }}}

// {{{ logging
static const char *const asd_log_level_names[] = {
    "ALWAYS", "TRACE", "DEBUG", " VERB", " INFO", " WARN", "ERROR", "SILENT",
};

static const int asd_log_root_rank = 0;
static int asd_log_rank = 0;
static FILE *asd_log_stream = NULL;
static int asd_log_threshold = ASD_LL_ALWAYS;

void asd_log_init(int rank, FILE *stream, int threshold)
{

  asd_log_rank = rank;
  if (stream == NULL)
  {
    asd_log_stream = stdout;
  }
  else
  {
    asd_log_stream = stream;
  }

  if (threshold == ASD_LL_DEFAULT)
  {
    asd_log_threshold = ASD_LL_INFO;
  }
  else
  {
    ASD_ABORT_IF_NOT(threshold <= ASD_LL_SILENT && threshold >= ASD_LL_ALWAYS,
                     "Invalid logging threshold");
    asd_log_threshold = threshold;
  }
}

void asd_log_printf(const char *file, int line, int category, int level, ...)
{
  FILE *stream = (asd_log_stream) ? asd_log_stream : stdout;
  va_list ap;
  const char *fmt;

  ASD_ASSERT(level <= ASD_LL_SILENT && level >= ASD_LL_ALWAYS);

  if (asd_log_threshold > level)
    return;

  if (category == ASD_LC_ROOT && asd_log_rank != asd_log_root_rank)
    return;

  if (category == ASD_LC_ALL)
  {
    fprintf(stream, "[%3d] ", asd_log_rank);
  }
  else
  {
    fprintf(stream, "[   ] ");
  }
  fprintf(stream, "%s: ", asd_log_level_names[level]);

  if (level == ASD_LL_TRACE)
  {

    fprintf(stream, "%s:%3d ", file, line);
  }

  va_start(ap, level);
  fmt = va_arg(ap, const char *);
  vfprintf(stream, fmt, ap);
  fprintf(stream, "\n");
  va_end(ap);
}
// }}}

// {{{ memory allocation and alignment
void *asd_malloc(size_t size)
{
  void *r;

  r = malloc(size);

  if (!r && !size)
  {
    r = malloc(1);
  }

  if (!r)
  {
    ASD_LERROR("Memory allocation of %lu bytes failed", (unsigned long)size);
    ASD_ABORT("Failed malloc");
  }

#ifdef ASD_DEBUG
  memset(r, 0xA3, size);
#endif

  return r;
}

void *asd_calloc(size_t nmemb, size_t size)
{
  void *r;

  r = calloc(nmemb, size);

  if (!r && (!nmemb || !size))
  {
    r = calloc(1, 1);
  }

  if (!r)
  {
    ASD_LERROR("Memory allocation of %lu elements of size %lu bytes failed",
               (unsigned long)nmemb, (unsigned long)size);
    ASD_ABORT("Failed calloc");
  }

  return r;
}

void *asd_realloc(void *ptr, size_t size)
{
  void *r;

  r = realloc(ptr, size);

  if (!r && !size)
  {
    r = realloc(ptr, 1);
  }

  if (!r)
  {
    ASD_LERROR("Memory reallocation of size %lu bytes failed",
               (unsigned long)size);
    ASD_ABORT("Failed realloc");
  }

  return r;
}

void asd_free(void *ptr) { free(ptr); }

static size_t asd_page_size()
{
  long page_size = sysconf(_SC_PAGE_SIZE);

  return (size_t)page_size;
}

#if defined(__APPLE__)
#include <sys/sysctl.h>
#include <sys/types.h>

static size_t asd_cache_line_size()
{
  size_t line_size = 0;
  size_t sizeof_line_size = sizeof(line_size);
  sysctlbyname("hw.cachelinesize", &line_size, &sizeof_line_size, 0, 0);
  return line_size;
}

#elif defined(__linux__)

static size_t asd_cache_line_size()
{
  long line_size = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);

  return (size_t)line_size;
}

#else
#error Unrecognized platform for cache line size
#endif

static void *asd_malloc_aligned_cache_line(size_t size, size_t line,
                                           size_t line_size, size_t page_size)
{
  void *r;
  intptr_t a;

  ASD_ASSERT(page_size >= 1);
  ASD_ASSERT(page_size >= line * line_size);

  r = asd_malloc(size + sizeof(intptr_t) + page_size);

  a = (((intptr_t)r + sizeof(intptr_t) + page_size - line * line_size - 1) /
       page_size) *
          page_size +
      line * line_size;

  ((intptr_t *)a)[-1] = (intptr_t)r;

  r = (void *)a;

  return r;
}

void *asd_malloc_aligned(size_t size)
{
  void *r;
  static size_t line_no = 0;

  const size_t line_size = asd_cache_line_size();
  const size_t page_size = asd_page_size();
  const size_t line_count = page_size / line_size;

  r = asd_malloc_aligned_cache_line(size, line_no, line_size, page_size);

  line_no = (line_no + 1) % line_count;

#ifdef ASD_DEBUG
  memset(r, 0xA3, size);
#endif

  return r;
}

void asd_free_aligned(void *ptr)
{
  ptr = (void *)((intptr_t *)ptr)[-1];
  ASD_ASSERT(ptr != NULL);
  free(ptr);
}
// }}}

// {{{ signals
#if defined __GNUC__
/*
 * This signal handler code is a modified version of the one presented in
 *
 *     http://spin.atomicobject.com/2013/01/13/exceptions-stack-traces-c/
 *
 * by Job Vranish.  The code can be found here
 *
 *     https://gist.github.com/4441299
 */

#define ASD_MAX_STACK_FRAMES 1024
static void *stack_traces[ASD_MAX_STACK_FRAMES];

static void asd_posix_print_stack_trace()
{
  int i, trace_size = 0;
  char **messages = (char **)NULL;

  trace_size = backtrace(stack_traces, ASD_MAX_STACK_FRAMES);
  messages = backtrace_symbols(stack_traces, trace_size);
  ASD_SYS_ERROR_CHECK(messages == NULL, "backtrace_symbols");

  for (i = 0; i < trace_size; ++i)
  {
    ASD_LERROR("%s", messages[i]);
  }

  if (messages)
  {
    free(messages);
  }
}

static void asd_posix_signal_handler(int sig, siginfo_t *siginfo, void *context)
{
  (void)context;
  switch (sig)
  {
  case SIGSEGV:
    ASD_LERROR("Caught SIGSEGV: Segmentation Fault");
    break;
  case SIGINT:
    ASD_LERROR("Caught SIGINT: Interactive attention signal, (usually ctrl+c)");
    break;
  case SIGFPE:
    switch (siginfo->si_code)
    {
    case FPE_INTDIV:
      ASD_LERROR("Caught SIGFPE: (integer divide by zero)");
      break;
    case FPE_INTOVF:
      ASD_LERROR("Caught SIGFPE: (integer overflow)");
      break;
    case FPE_FLTDIV:
      ASD_LERROR("Caught SIGFPE: (floating-point divide by zero)");
      break;
    case FPE_FLTOVF:
      ASD_LERROR("Caught SIGFPE: (floating-point overflow)");
      break;
    case FPE_FLTUND:
      ASD_LERROR("Caught SIGFPE: (floating-point underflow)");
      break;
    case FPE_FLTRES:
      ASD_LERROR("Caught SIGFPE: (floating-point inexact result)");
      break;
    case FPE_FLTINV:
      ASD_LERROR("Caught SIGFPE: (floating-point invalid operation)");
      break;
    case FPE_FLTSUB:
      ASD_LERROR("Caught SIGFPE: (subscript out of range)");
      break;
    default:
      ASD_LERROR("Caught SIGFPE: Arithmetic Exception");
      break;
    }
  case SIGILL:
    switch (siginfo->si_code)
    {
    case ILL_ILLOPC:
      ASD_LERROR("Caught SIGILL: (illegal opcode)");
      break;
    case ILL_ILLOPN:
      ASD_LERROR("Caught SIGILL: (illegal operand)");
      break;
    case ILL_ILLADR:
      ASD_LERROR("Caught SIGILL: (illegal addressing mode)");
      break;
    case ILL_ILLTRP:
      ASD_LERROR("Caught SIGILL: (illegal trap)");
      break;
    case ILL_PRVOPC:
      ASD_LERROR("Caught SIGILL: (privileged opcode)");
      break;
    case ILL_PRVREG:
      ASD_LERROR("Caught SIGILL: (privileged register)");
      break;
    case ILL_COPROC:
      ASD_LERROR("Caught SIGILL: (coprocessor error)");
      break;
    case ILL_BADSTK:
      ASD_LERROR("Caught SIGILL: (internal stack error)");
      break;
    default:
      ASD_LERROR("Caught SIGILL: Illegal Instruction");
      break;
    }
    break;
  case SIGTERM:
    ASD_LERROR("Caught SIGTERM: a termination request was sent to the program");
    break;
  case SIGABRT:
    ASD_LERROR("Caught SIGABRT: usually caused by an abort() or assert()");
    break;
  default:
    break;
  }
  asd_posix_print_stack_trace();
  _Exit(1);
}

static uint8_t alternate_stack[SIGSTKSZ];
void asd_signal_handler_set()
{
  /* setup alternate stack */
  {
    stack_t ss;
    /* malloc is usually used here, I'm not 100% sure my static allocation
       is valid but it seems to work just fine. */
    ss.ss_sp = (void *)alternate_stack;
    ss.ss_size = SIGSTKSZ;
    ss.ss_flags = 0;

    if (sigaltstack(&ss, NULL) != 0)
    {
      err(1, "sigaltstack");
    }
  }

  /* register our signal handlers */
  {
    struct sigaction sig_action;
    sig_action.sa_sigaction = asd_posix_signal_handler;
    sigemptyset(&sig_action.sa_mask);

#ifdef __APPLE__
    /* for some reason we backtrace() doesn't work on osx
       when we use an alternate stack */
    sig_action.sa_flags = SA_SIGINFO;
#else
    sig_action.sa_flags = SA_SIGINFO | SA_ONSTACK;
#endif

    if (sigaction(SIGSEGV, &sig_action, NULL) != 0)
    {
      err(1, "sigaction");
    }
    if (sigaction(SIGFPE, &sig_action, NULL) != 0)
    {
      err(1, "sigaction");
    }
    if (sigaction(SIGINT, &sig_action, NULL) != 0)
    {
      err(1, "sigaction");
    }
    if (sigaction(SIGILL, &sig_action, NULL) != 0)
    {
      err(1, "sigaction");
    }
    if (sigaction(SIGTERM, &sig_action, NULL) != 0)
    {
      err(1, "sigaction");
    }
    if (sigaction(SIGABRT, &sig_action, NULL) != 0)
    {
      err(1, "sigaction");
    }
  }
}

#else

void asd_almost_c99_signal_handler(int sig)
{
  switch (sig)
  {
  case SIGABRT:
    fputs("Caught SIGABRT: usually caused by an abort() or assert()\n", stderr);
    break;
  case SIGFPE:
    fputs("Caught SIGFPE: arithmetic exception, such as divide by zero\n",
          stderr);
    break;
  case SIGILL:
    fputs("Caught SIGILL: illegal instruction\n", stderr);
    break;
  case SIGINT:
    fputs("Caught SIGINT: interactive attention signal, probably a ctrl+c\n",
          stderr);
    break;
  case SIGSEGV:
    fputs("Caught SIGSEGV: segfault\n", stderr);
    break;
  case SIGTERM:
  default:
    fputs("Caught SIGTERM: a termination request was sent to the program\n",
          stderr);
    break;
  }
  _Exit(1);
}

void asd_signal_handler_set()
{
  signal(SIGABRT, asd_almost_c99_signal_handler);
  signal(SIGFPE, asd_almost_c99_signal_handler);
  signal(SIGILL, asd_almost_c99_signal_handler);
  signal(SIGINT, asd_almost_c99_signal_handler);
  signal(SIGSEGV, asd_almost_c99_signal_handler);
  signal(SIGTERM, asd_almost_c99_signal_handler);
}
#endif
// }}}


// {{{ critbit
//
// We use a slightly modified version of Adam Langley's implementation of the
// binary crit-bit from Dan Bernstein's qhasm.  The original file were
// obtained from [here](https://github.com/agl/critbit.git).

typedef struct
{
  void *child[2];
  uint32_t byte;
  uint8_t otherbits;
} asd_critbit0_node_t;

int asd_critbit0_contains(asd_critbit0_tree_t *t, const char *u)
{
  const uint8_t *ubytes = (uint8_t *)u;
  const size_t ulen = strlen(u);
  uint8_t *p = (uint8_t*)t->root;

  if (!p)
    return 0;

  while (1 & (intptr_t)p)
  {
    asd_critbit0_node_t *q = (asd_critbit0_node_t *)(p - 1);

    uint8_t c = 0;
    if (q->byte < ulen)
      c = ubytes[q->byte];
    const int direction = (1 + (q->otherbits | c)) >> 8;

    p = (uint8_t*) q->child[direction];
  }

  return 0 == strcmp(u, (const char *)p);
}

int asd_critbit0_insert(asd_critbit0_tree_t *t, const char *u)
{
  const uint8_t *const ubytes = (uint8_t *)u;
  const size_t ulen = strlen(u);
  uint8_t *p = (uint8_t *)t->root;

  if (!p)
  {
    char *x;
    int a = posix_memalign((void **)&x, sizeof(void *), ulen + 1);
    if (a)
      return 0;
    memcpy(x, u, ulen + 1);
    t->root = x;
    return 2;
  }

  while (1 & (intptr_t)p)
  {
    asd_critbit0_node_t *q = (asd_critbit0_node_t *)(p - 1);

    uint8_t c = 0;
    if (q->byte < ulen)
      c = ubytes[q->byte];
    const int direction = (1 + (q->otherbits | c)) >> 8;

    p = (uint8_t *)q->child[direction];
  }

  uint32_t newbyte;
  uint32_t newotherbits;

  for (newbyte = 0; newbyte < ulen; ++newbyte)
  {
    if (p[newbyte] != ubytes[newbyte])
    {
      newotherbits = p[newbyte] ^ ubytes[newbyte];
      goto different_byte_found;
    }
  }

  if (p[newbyte] != 0)
  {
    newotherbits = p[newbyte];
    goto different_byte_found;
  }
  return 1;

different_byte_found:

  newotherbits |= newotherbits >> 1;
  newotherbits |= newotherbits >> 2;
  newotherbits |= newotherbits >> 4;
  newotherbits = (newotherbits & ~(newotherbits >> 1)) ^ 255;
  uint8_t c = p[newbyte];
  int newdirection = (1 + (newotherbits | c)) >> 8;

  asd_critbit0_node_t *newnode;

  if (posix_memalign((void **)&newnode, sizeof(void *),
                     sizeof(asd_critbit0_node_t)))
  {
    return 0;
  }

  char *x;
  if (posix_memalign((void **)&x, sizeof(void *), ulen + 1))
  {
    free(newnode);
    return 0;
  }
  memcpy(x, ubytes, ulen + 1);

  newnode->byte = newbyte;
  newnode->otherbits = (int8_t)newotherbits;
  newnode->child[1 - newdirection] = x;

  void **wherep = &t->root;
  for (;;)
  {
    uint8_t *p = (uint8_t *)*wherep;
    if (!(1 & (intptr_t)p))
      break;
    asd_critbit0_node_t *q = (asd_critbit0_node_t *)(p - 1);
    if (q->byte > newbyte)
      break;
    if (q->byte == newbyte && q->otherbits > newotherbits)
      break;
    uint8_t c = 0;
    if (q->byte < ulen)
      c = ubytes[q->byte];
    const int direction = (1 + (q->otherbits | c)) >> 8;
    wherep = q->child + direction;
  }

  newnode->child[newdirection] = *wherep;
  *wherep = (void *)(1 + (char *)newnode);

  return 2;
}

int asd_critbit0_delete(asd_critbit0_tree_t *t, const char *u)
{
  const uint8_t *ubytes = (uint8_t *)u;
  const size_t ulen = strlen(u);
  uint8_t *p = (uint8_t *)t->root;
  void **wherep = &t->root;
  void **whereq = 0;
  asd_critbit0_node_t *q = 0;
  int direction = 0;

  if (!p)
    return 0;

  while (1 & (intptr_t)p)
  {
    whereq = wherep;
    q = (asd_critbit0_node_t *)(p - 1);
    uint8_t c = 0;
    if (q->byte < ulen)
      c = ubytes[q->byte];
    direction = (1 + (q->otherbits | c)) >> 8;
    wherep = q->child + direction;
    p = (uint8_t *)*wherep;
  }

  if (0 != strcmp(u, (const char *)p))
    return 0;
  free(p);

  if (!whereq)
  {
    t->root = 0;
    return 1;
  }

  *whereq = q->child[1 - direction];
  free(q);

  return 1;
}

static void traverse(void *top)
{
  uint8_t *p = (uint8_t *)top;

  if (1 & (intptr_t)p)
  {
    asd_critbit0_node_t *q = (asd_critbit0_node_t *)(p - 1);
    traverse(q->child[0]);
    traverse(q->child[1]);
    free(q);
  }
  else
  {
    free(p);
  }
}

void asd_critbit0_clear(asd_critbit0_tree_t *t)
{
  if (t->root)
    traverse(t->root);
  t->root = NULL;
}

static int allprefixed_traverse(uint8_t *top,
                                int (*handle)(const char *, void *), void *arg)
{

  if (1 & (intptr_t)top)
  {
    asd_critbit0_node_t *q = (asd_critbit0_node_t *)(top - 1);
    for (int direction = 0; direction < 2; ++direction)
      switch (allprefixed_traverse((uint8_t *)q->child[direction], handle, arg))
      {
      case 1:
        break;
      case 0:
        return 0;
      default:
        return -1;
      }
    return 1;
  }

  return handle((const char *)top, arg);
}

int asd_critbit0_allprefixed(asd_critbit0_tree_t *t, const char *prefix,
                             int (*handle)(const char *, void *), void *arg)
{
  const uint8_t *ubytes = (uint8_t *)prefix;
  const size_t ulen = strlen(prefix);
  uint8_t *p = (uint8_t *)t->root;
  uint8_t *top = p;

  if (!p)
    return 1;

  while (1 & (intptr_t)p)
  {
    asd_critbit0_node_t *q = (asd_critbit0_node_t *)(p - 1);
    uint8_t c = 0;
    if (q->byte < ulen)
      c = ubytes[q->byte];
    const int direction = (1 + (q->otherbits | c)) >> 8;
    p = (uint8_t *)q->child[direction];
    if (q->byte < ulen)
      top = p;
  }

  for (size_t i = 0; i < ulen; ++i)
  {
    if (p[i] != ubytes[i])
      return 1;
  }

  return allprefixed_traverse(top, handle, arg);
}
// }}}

// {{{ dictionary

#define ASD_KEYVALUE_SPLIT '\255'
#define ASD_PTR_STR_LEN ASD_BUFSIZ

void asd_dictionary_init(asd_dictionary_t *d)
{
  d->num_entries = 0;
  d->t.root = NULL;
}

void asd_dictionary_clear(asd_dictionary_t *d)
{
  asd_critbit0_clear(&(d->t));
  d->num_entries = 0;
}

static int asd_dictionary_contains_check(const char *value, void *arg)
{
  (*(int *)arg)++;
  return 1;
}

int asd_dictionary_contains(asd_dictionary_t *d, const char *key)
{
  const size_t keylen = strlen(key);

  char *u = (char *)asd_malloc(sizeof(char) * (keylen + 2));

  memcpy(u, key, keylen);
  u[keylen] = ASD_KEYVALUE_SPLIT;
  u[keylen + 1] = '\0';
  int found = 0;
  asd_critbit0_allprefixed(&(d->t), u, &asd_dictionary_contains_check, &found);
  asd_free(u);
  return found;
}

int asd_dictionary_insert(asd_dictionary_t *d, const char *key, const char *val)
{
  const size_t keylen = strlen(key);
  const size_t vallen = strlen(val);

  char *keyval = (char *)asd_malloc(sizeof(char) * (keylen + vallen + 2));

  memcpy(keyval, key, keylen);

  if (asd_dictionary_contains(d, key))
  {
    free(keyval);
    return 1;
  }

  keyval[keylen] = ASD_KEYVALUE_SPLIT;
  memcpy(&keyval[keylen + 1], val, vallen);
  keyval[keylen + vallen + 1] = '\0';

  int rval = asd_critbit0_insert(&(d->t), keyval);

  if (rval == 2)
    ++d->num_entries;

  asd_free(keyval);
  return rval;
}

int asd_dictionary_insert_ptr(asd_dictionary_t *d, const char *key,
                              const void *val_ptr)
{
  char val_str[ASD_PTR_STR_LEN + 1];
  snprintf(val_str, ASD_PTR_STR_LEN + 1, "%p", val_ptr);
  return asd_dictionary_insert(d, key, val_str);
}

int asd_dictionary_insert_int(asd_dictionary_t *d, const char *key,
                              const int val)
{
  char val_str[ASD_BUFSIZ];
  snprintf(val_str, ASD_BUFSIZ, "%d", val);
  return asd_dictionary_insert(d, key, val_str);
}

int asd_dictionary_insert_int32(asd_dictionary_t *d, const char *key,
                                const int32_t val)
{
  char val_str[ASD_BUFSIZ];
  snprintf(val_str, ASD_BUFSIZ, "%" PRId32, val);
  return asd_dictionary_insert(d, key, val_str);
}

int asd_dictionary_insert_size_t(asd_dictionary_t *d, const char *key,
                                 const size_t val)
{
  char val_str[ASD_BUFSIZ];
  snprintf(val_str, ASD_BUFSIZ, "%zu", val);
  return asd_dictionary_insert(d, key, val_str);
}

static int asd_dictionary_get_value_handle(const char *keyval, void *arg)
{
  char *key = (char *)((void **)arg)[0];
  char **val = (char **)((void **)arg)[1];
  const size_t keylen = strlen(key);

  *val = (char *)&keyval[keylen];

  return 1;
}

char *asd_dictionary_get_value(asd_dictionary_t *d, const char *key)
{
  if (!asd_dictionary_contains(d, key))
    return NULL;

  const size_t keylen = strlen(key);

  char *u = (char *)asd_malloc(sizeof(char) * (keylen + 2));

  memcpy(u, key, keylen);
  u[keylen] = ASD_KEYVALUE_SPLIT;
  u[keylen + 1] = '\0';

  char *value = NULL;
  void *arg[2];
  arg[0] = u;
  arg[1] = &value;

  asd_critbit0_allprefixed(&(d->t), u, &asd_dictionary_get_value_handle, arg);

  asd_free(u);

  return value;
}

void *asd_dictionary_get_value_ptr(asd_dictionary_t *d, const char *key)
{
  char *val_str = asd_dictionary_get_value(d, key);
  if (val_str == NULL)
    return NULL;
  void *val_ptr = NULL;
  sscanf(val_str, "%p", &val_ptr);
  return val_ptr;
}

int asd_dictionary_get_value_int(asd_dictionary_t *d, const char *key, int *val)
{
  char *val_str = asd_dictionary_get_value(d, key);
  if (val_str == NULL)
    return 0;
  int n = sscanf(val_str, "%d", val);
  return n;
}

int asd_dictionary_get_value_int32(asd_dictionary_t *d, const char *key,
                                   int32_t *val)
{
  char *val_str = asd_dictionary_get_value(d, key);
  if (val_str == NULL)
    return 0;
  int n = sscanf(val_str, "%" SCNd32, val);
  return n;
}

int asd_dictionary_delete(asd_dictionary_t *d, const char *key)
{
  const size_t keylen = strlen(key);
  char *val = asd_dictionary_get_value(d, key);
  if (val == NULL)
    return 0;
  const char *keyval = val - keylen - 1;
  return asd_critbit0_delete(&(d->t), keyval);
}

typedef struct
{
  int (*handle)(const char *, const char *, void *);
  void *arg;
} asd_dict_allprex;

static int asd_dictionary_allprefixed_usercall(const char *keyval, void *arg)
{
  char *key = (char *)keyval;
  char *split = strchr(key, ASD_KEYVALUE_SPLIT);
  *split = '\0';
  asd_dict_allprex *s_arg = (asd_dict_allprex *)arg;
  s_arg->handle(key, split + 1, s_arg->arg);
  *split = ASD_KEYVALUE_SPLIT;
  return 1;
}

int asd_dictionary_allprefixed(asd_dictionary_t *d, const char *prefix,
                               int (*handle)(const char *, const char *,
                                             void *),
                               void *arg)
{
  asd_dict_allprex args = {0, 0};
  args.handle = handle;
  args.arg = arg;
  return asd_critbit0_allprefixed(&(d->t), prefix,
                                  &asd_dictionary_allprefixed_usercall, &args);
}

typedef struct
{
  int (*handle)(const char *, void *, void *);
  void *arg;
} asd_dict_allprex_ptr;

static int asd_dictionary_allprefixed_usercall_ptr(const char *keyval,
                                                   void *arg)
{
  char *key = (char *)keyval;
  char *split = strchr(key, ASD_KEYVALUE_SPLIT);
  *split = '\0';

  void *val_ptr = NULL;
  sscanf(split + 1, "%p", &val_ptr);

  asd_dict_allprex_ptr *s_arg = (asd_dict_allprex_ptr *)arg;
  s_arg->handle(key, val_ptr, s_arg->arg);

  *split = ASD_KEYVALUE_SPLIT;
  return 1;
}

int asd_dictionary_allprefixed_ptr(asd_dictionary_t *d, const char *prefix,
                                   int (*handle)(const char *, void *, void *),
                                   void *arg)
{
  asd_dict_allprex_ptr args = {0, 0};
  args.handle = handle;
  args.arg = arg;
  return asd_critbit0_allprefixed(
      &(d->t), prefix, &asd_dictionary_allprefixed_usercall_ptr, &args);
}

typedef struct
{
  int (*handle)(const char *, int, void *);
  void *arg;
} asd_dict_allprex_int;

static int asd_dictionary_allprefixed_usercall_int(const char *keyval,
                                                   void *arg)
{
  char *key = (char *)keyval;
  char *split = strchr(key, ASD_KEYVALUE_SPLIT);
  *split = '\0';

  int val_int;
  sscanf(split + 1, "%d", &val_int);

  asd_dict_allprex_int *s_arg = (asd_dict_allprex_int *)arg;
  s_arg->handle(key, val_int, s_arg->arg);

  *split = ASD_KEYVALUE_SPLIT;
  return 1;
}

int asd_dictionary_allprefixed_int(asd_dictionary_t *d, const char *prefix,
                                   int (*handle)(const char *, int, void *),
                                   void *arg)
{
  asd_dict_allprex_int args = {0, 0};
  args.handle = handle;
  args.arg = arg;
  return asd_critbit0_allprefixed(
      &(d->t), prefix, &asd_dictionary_allprefixed_usercall_int, &args);
}

typedef struct
{
  int (*handle)(const char *, size_t, void *);
  void *arg;
} asd_dict_allprex_size_t;

static int asd_dictionary_allprefixed_usercall_size_t(const char *keyval,
                                                      void *arg)
{
  char *key = (char *)keyval;
  char *split = strchr(key, ASD_KEYVALUE_SPLIT);
  *split = '\0';

  size_t val_size_t;
  sscanf(split + 1, "%zu", &val_size_t);

  asd_dict_allprex_size_t *s_arg = (asd_dict_allprex_size_t *)arg;
  s_arg->handle(key, val_size_t, s_arg->arg);

  *split = ASD_KEYVALUE_SPLIT;
  return 1;
}

int asd_dictionary_allprefixed_size_t(asd_dictionary_t *d, const char *prefix,
                                      int (*handle)(const char *, size_t,
                                                    void *),
                                      void *arg)
{
  asd_dict_allprex_size_t args = {0, 0};
  args.handle = handle;
  args.arg = arg;
  return asd_critbit0_allprefixed(
      &(d->t), prefix, &asd_dictionary_allprefixed_usercall_size_t, &args);
}

// }}}

// {{{ random number (pcg32)
// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
// Copied from: http://www.pcg-random.org/
// Copied from: https://github.com/imneme/pcg-c-basic

static uint32_t asd_pcg32_random_r(asd_pcg32_random_t *rng)
{
  uint64_t oldstate = rng->state;
  // Advance internal state
  rng->state = oldstate * 6364136223846793005ULL + (rng->inc | 1);
  // Calculate output function (XSH RR), uses old state for max ILP
  uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
  uint32_t rot = (uint32_t)(oldstate >> 59u);
  return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

void asd_pcg32_srandom_r(asd_pcg32_random_t *rng, uint64_t initstate,
                         uint64_t initseq)
{
  rng->state = 0U;
  rng->inc = (initseq << 1u) | 1u;
  asd_pcg32_random_r(rng);
  rng->state += initstate;
  asd_pcg32_random_r(rng);
}

uint32_t asd_pcg32_boundedrand_r(asd_pcg32_random_t *rng, uint32_t bound)
{
  // To avoid bias, we need to make the range of the RNG a multiple of
  // bound, which we do by dropping output less than a threshold.
  // A naive scheme to calculate the threshold would be to do
  //
  //     uint32_t threshold = 0x100000000ull % bound;
  //
  // but 64-bit div/mod is slower than 32-bit div/mod (especially on
  // 32-bit platforms).  In essence, we do
  //
  //     uint32_t threshold = (0x100000000ull-bound) % bound;
  //
  // because this version will calculate the same modulus, but the LHS
  // value is less than 2^32.

  uint32_t threshold = -bound % bound;

  // Uniformity guarantees that this loop will terminate.  In practice, it
  // should usually terminate quickly; on average (assuming all bounds are
  // equally likely), 82.25% of the time, we can expect it to require just
  // one iteration.  In the worst case, someone passes a bound of 2^31 + 1
  // (i.e., 2147483649), which invalidates almost 50% of the range.  In
  // practice, bounds are typically small and only a tiny amount of the range
  // is eliminated.
  for (;;)
  {
    uint32_t r = asd_pcg32_random_r(rng);
    if (r >= threshold)
      return r % bound;
  }
}
// }}}

// {{{ file
/*
 * This uses a posix compilent solution from:
 *   https://www.securecoding.cert.org/confluence/display/c/FIO19-C.+Do+not+use+fseek%28%29+and+ftell%28%29+to+compute+the+size+of+a+regular+file
 */
static size_t asd_file_size(const char *filename)
{
  struct stat stbuf;
  int fd, err;

  fd = open(filename, O_RDONLY);
  if (fd == -1)
  {
    ASD_LERROR("Error opening file: %s", filename);
    ASD_ABORT("Failed open");
  }

  if ((fstat(fd, &stbuf) != 0) || (!S_ISREG(stbuf.st_mode)))
  {
    ASD_LERROR("Error determining size of the file: %s", filename);
    ASD_ABORT("Failed fstat");
  }

  ASD_ASSERT(stbuf.st_size >= 0);

  err = close(fd);
  if (err)
  {
    ASD_LERROR("Error closing file: %s", filename);
    ASD_ABORT("Failed close");
  }

  return (size_t)stbuf.st_size;
}

char *asd_read_file(const char *filename, size_t *len)
{
  size_t readsize;
  char *buffer;
  size_t filesize = asd_file_size(filename);
  FILE *stream = fopen(filename, "r");
  if (!stream)
  {
    ASD_LERROR("Error opening the file: %s", filename);
    ASD_ABORT("Failed fopen");
  }

  buffer = (char *)asd_malloc(filesize + 1);
  readsize = fread(buffer, sizeof(char), filesize, stream);
  buffer[filesize] = '\0';

  if (readsize != filesize)
  {
    ASD_LERROR("Error determining reading the file: %s", filename);
    ASD_ABORT("Failed fread");
  }

  if (fclose(stream))
  {
    ASD_LERROR("Error closing the file: %s", filename);
    ASD_ABORT("Failed fclose");
  }

  if (len)
    *len = filesize;

  return buffer;
}
// }}}

// {{{ endian
/*
 *  The follow code snippet is from:
 *
 *    http://www.ibm.com/developerworks/aix/library/au-endianc/
 */
int asd_endian()
{
  int i = 1;
  char *p = (char *)&i;

  if (p[0] == 1)
    return ASD_LITTLE_ENDIAN;
  else
    return ASD_BIG_ENDIAN;
}
// }}}

// {{{ MPI
int asd_get_host_rank(MPI_Comm comm)
{
  int host_rank, rank, size, length;
  char name[MPI_MAX_PROCESSOR_NAME];

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Get_processor_name(name, &length);

  int *host_ranks = (int *)asd_malloc(size * sizeof(int));
  char *allnames = (char*)asd_malloc(size * MPI_MAX_PROCESSOR_NAME * sizeof(char));

  ASD_MPI_CHECK(MPI_Allgather(name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, allnames,
                              MPI_MAX_PROCESSOR_NAME, MPI_CHAR, comm));

  asd_dictionary_t hosts;
  asd_dictionary_init(&hosts);

  for (int h = 0; h < size; ++h)
  {
    const char *key = allnames + h * MPI_MAX_PROCESSOR_NAME;
    int hr = 0;
    if (asd_dictionary_get_value_int(&hosts, key, &hr))
    {
      ++hr;
      if (1 != asd_dictionary_delete(&hosts, key))
        ASD_ABORT("get host rank dictionary delete fail");
    }
    if (2 != asd_dictionary_insert_int(&hosts, key, hr))
      ASD_ABORT("get host rank dictionary insert fail");

    host_ranks[h] = hr;
  }

  ASD_ROOT_VERBOSE("-----------Host Ranks----------");
  for (int h = 0; h < size; ++h)
  {
    ASD_ROOT_VERBOSE("  host_rank[%3d] = %2d", h, host_ranks[h]);
  }
  ASD_ROOT_VERBOSE("-------------------------------");

  host_rank = host_ranks[rank];

  asd_dictionary_clear(&hosts);
  asd_free(allnames);
  asd_free(host_ranks);

  return host_rank;
}
// }}}

// {{{ Lua
#ifdef ASD_USE_LUA

#define ASD_LUA_MAX_COMMAND_LEN 4096
#define ASD_LUA_EVALEXP_VAR "XXX_asd_evalexpr_XXX"

int asd_lua_global_function_call(lua_State *L, const char *name,
                                 const char *sig, ...)
{
  va_list vl;
  int num_arg = 0;
  int num_res = 0;

  va_start(vl, sig);

  char buf[ASD_LUA_MAX_COMMAND_LEN];
  /* Assign the Lua expression to a Lua global variable. */
  snprintf(buf, ASD_LUA_MAX_COMMAND_LEN, ASD_LUA_EVALEXP_VAR "=%s", name);
  if (!luaL_dostring(L, buf))
  {
    /* Get the value of the global varibable */
    lua_getglobal(L, ASD_LUA_EVALEXP_VAR);
  }
  else
  {
    ASD_ROOT_WARNING("function `%s' not found in lua file", name);
    return 1;
  }

  if (!lua_isfunction(L, -1))
  {
    ASD_ROOT_WARNING("function `%s' not found in lua file", name);
    lua_pop(L, 1);
    return 1;
  }

  for (num_arg = 0; sig[num_arg] && sig[num_arg] != '>'; num_arg++)
  {
    luaL_checkstack(L, 1, "too many arguments");

    switch (sig[num_arg])
    {
    case 'd':
      lua_pushnumber(L, va_arg(vl, double));
      break;
    case 'l':
      lua_pushnumber(L, (double)va_arg(vl, long double));
      break;
    case 'i':
      lua_pushinteger(L, va_arg(vl, int));
      break;
    case 's':
      lua_pushstring(L, va_arg(vl, char *));
      break;
    case '>':
      break;
    default:
      ASD_ABORT("function '%s' invalid input argument (%c)", name,
                sig[num_arg]);
    }
  }

  ASD_ABORT_IF_NOT(sig[num_arg] == '>', "arguments for '%s' does not contain "
                                        " a '>' character",
                   name);

  num_res = (int)strlen(sig) - num_arg - 1;

  ASD_ABORT_IF_NOT(lua_pcall(L, num_arg, num_res, 0) == 0,
                   "error running function %s: %s", name, lua_tostring(L, -1));

  for (int n = 0; n < num_res; n++)
  {
    switch (sig[num_arg + 1 + n])
    {
    case 'd':
      ASD_ABORT_IF_NOT(lua_isnumber(L, n - num_res),
                       "for '%s' return %d expected number got '%s'", name, n,
                       lua_tostring(L, n - num_res));
      *va_arg(vl, double *) = (double)lua_tonumber(L, n - num_res);
      break;
    case 'l':
      ASD_ABORT_IF_NOT(lua_isnumber(L, n - num_res),
                       "for '%s' return %d expected number got '%s'", name, n,
                       lua_tostring(L, n - num_res));
      *va_arg(vl, long double *) = (long double)lua_tonumber(L, n - num_res);
      break;
    case 'i':
      ASD_ABORT_IF_NOT(lua_isnumber(L, n - num_res),
                       "for '%s' return %d expected number got '%s'", name, n,
                       lua_tostring(L, n - num_res));
      *va_arg(vl, int *) = (int)lua_tointeger(L, n - num_res);
      break;
    case 's':
      ASD_ABORT_IF_NOT(lua_isstring(L, n - num_res),
                       "for '%s' return %d expected string got '%s'", name, n,
                       lua_tostring(L, n - num_res));
      *va_arg(vl, const char **) = lua_tostring(L, n - num_res);
      break;
    default:
      ASD_ABORT("function '%s' invalid output argument (%c)", name,
                sig[num_arg]);
    }
  }

  lua_pop(L, num_res);

  va_end(vl);
  return 0;
}

int asd_lua_expr_boolean(lua_State *L, const char *expr, int def)
{
  int r = def;
  char buf[ASD_LUA_MAX_COMMAND_LEN];
  snprintf(buf, ASD_LUA_MAX_COMMAND_LEN, ASD_LUA_EVALEXP_VAR "=%s", expr);
  if (!luaL_dostring(L, buf))
  {
    /* Get the value of the global varibable */
    lua_getglobal(L, ASD_LUA_EVALEXP_VAR);
    if (lua_isboolean(L, -1))
      r = lua_toboolean(L, -1);
  }
  lua_pop(L, 1);

  return r;
}

lua_Integer asd_lua_expr_integer(lua_State *L, const char *expr,
                                 lua_Integer def)
{
  lua_Integer r = def;
  char buf[ASD_LUA_MAX_COMMAND_LEN];
  /* Assign the Lua expression to a Lua global variable. */
  snprintf(buf, ASD_LUA_MAX_COMMAND_LEN, ASD_LUA_EVALEXP_VAR "=%s", expr);
  if (!luaL_dostring(L, buf))
  {
    /* Get the value of the global varibable */
    lua_getglobal(L, ASD_LUA_EVALEXP_VAR);
    if (lua_isnumber(L, -1))
      r = lua_tointeger(L, -1);
  }
  lua_pop(L, 1);

  return r;
}

lua_Number asd_lua_expr_number(lua_State *L, const char *expr, lua_Number def)
{
  lua_Number r = def;
  char buf[ASD_LUA_MAX_COMMAND_LEN];
  /* Assign the Lua expression to a Lua global variable. */
  snprintf(buf, ASD_LUA_MAX_COMMAND_LEN, ASD_LUA_EVALEXP_VAR "=%s", expr);
  if (!luaL_dostring(L, buf))
  {
    /* Get the value of the global varibable */
    lua_getglobal(L, ASD_LUA_EVALEXP_VAR);
    if (lua_isnumber(L, -1))
      r = lua_tonumber(L, -1);
  }
  lua_pop(L, 1);

  return r;
}

char *asd_lua_expr_string(lua_State *L, const char *expr, const char *def)
{
  const char *r = def;
  size_t len = strlen(r);
  char buf[ASD_LUA_MAX_COMMAND_LEN];
  /* Assign the Lua expression to a Lua global variable. */
  snprintf(buf, ASD_LUA_MAX_COMMAND_LEN, ASD_LUA_EVALEXP_VAR "=%s", expr);
  if (!luaL_dostring(L, buf))
  {
    /* Get the value of the global varibable */
    lua_getglobal(L, ASD_LUA_EVALEXP_VAR);
    if (lua_isstring(L, -1))
      r = lua_tolstring(L, -1, &len);
  }

  lua_pop(L, 1);

  char *s = asd_malloc((len + 1) * sizeof(char));
  strncpy(s, r, len + 1);
  s[len] = '\0';

  return s;
}

#endif
// }}}

// {{{ linear algebra
static void asd_la_mtranspose(size_t m, size_t n, long double *A,
                              size_t lda, long double *B, size_t ldb)
{
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < n; ++j)
      B[i * ldb + j] = A[j * lda + i];
}

static void asd_la_mTmmult(size_t m, size_t n, size_t k,
                           long double *A, size_t lda,
                           long double *B, size_t ldb,
                           long double *C, size_t ldc)
{
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < n; ++j)
      for (size_t p = 0; p < k; ++p)
        C[j * ldc + i] += A[i * lda + p] * B[j * ldb + p];
}

static void asd_la_mmmult(size_t m, size_t n, size_t k, long double *A,
                          size_t lda, long double *B, size_t ldb,
                          long double *C, size_t ldc)
{
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < n; ++j)
      for (size_t p = 0; p < k; ++p)
        C[j * ldc + i] += A[p * lda + i] * B[j * ldb + p];
}

static void asd_la_lu_factor(size_t n, long double *A,
                             size_t *p, size_t *q)
{
  /*
   * Algorithm 3.4.2 (Gaussian Elimination with Complete Pivoting) from p. 118
   * of Matrix Computations, Third Edition, by Golub and van Loan.
   */

  /*
   * Initialize pivots
   */
  for (size_t k = 0; k < n; ++k)
  {
    p[k] = q[k] = k;
  }

  for (size_t k = 0; k < n - 1; ++k)
  {
    size_t mu = k, lambda = k;
    long double a_max = 0.0;

    for (size_t j = k; j < n; ++j)
    {
      for (size_t i = k; i < n; ++i)
      {
        const long double a_abs = fabsl(A[i + n * j]);
        if (a_abs > a_max)
        {
          a_max = a_abs;
          mu = i;
          lambda = j;
        }
      }
    }

    /*
     * Swap rows
     */
    const size_t ptmp = p[k];
    p[k] = p[mu];
    p[mu] = ptmp;
    for (size_t j = 0; j < n; ++j)
    {
      const long double rtmp = A[k + n * j];
      A[k + n * j] = A[mu + n * j];
      A[mu + n * j] = rtmp;
    }

    /*
     * Swap columns
     */
    const size_t qtmp = q[k];
    q[k] = q[lambda];
    q[lambda] = qtmp;
    for (size_t i = 0; i < n; ++i)
    {
      const long double rtmp = A[i + n * k];
      A[i + n * k] = A[i + n * lambda];
      A[i + n * lambda] = rtmp;
    }

    if (fabsl(A[k + n * k]) > (LDBL_EPSILON * LDBL_EPSILON))
    {
      for (size_t i = k + 1; i < n; ++i)
      {
        A[i + n * k] = A[i + n * k] / A[k + n * k];
      }

      for (size_t i = k + 1; i < n; ++i)
      {
        for (size_t j = k + 1; j < n; ++j)
        {
          A[i + n * j] = A[i + n * j] - A[i + n * k] * A[k + n * j];
        }
      }
    }
  }
}

static void asd_la_lu_solve(size_t n, long double *LU,
                            long double *x, size_t *p,
                            size_t *q, long double *work)
{
  /*
   * Compute $Pb$
   */
  for (size_t k = 0; k < n; ++k)
  {
    work[k] = x[p[k]];
  }

  /*
   * Algorithm 3.1.3 (Forward Substitution: Column Version) from p. 90 of
   * Matrix Computations, Third Edition, by Golub and van Loan.
   *
   * Note: here we have L is unit lower triangular.
   *
   * Solve $Ly=Pb$.
   */
  for (size_t j = 0; j < n - 1; ++j)
  {
    /*
     * work[j] = work[j] / LU[j + n * j];
     */
    for (size_t i = j + 1; i < n; ++i)
    {
      work[i] = work[i] - work[j] * LU[i + n * j];
    }
  }
  /*
   * work[n - 1] = work[n - 1] / LU[n - 1 + n * (n - 1)];
   */

  /*
   * Algorithm 3.1.4 (Back Substitution: Column Version) from p. 90 of
   * Matrix Computations, Third Edition, by Golub and van Loan.
   *
   * Solve $Uw=y$.
   */
  for (size_t j = n - 1; j > 0; --j)
  {
    work[j] = work[j] / LU[j + n * j];
    for (size_t i = 0; i < j; ++i)
    {
      work[i] = work[i] - work[j] * LU[i + n * j];
    }
  }
  work[0] = work[0] / LU[0 + n * 0];

  /*
   * Compute $Qw$
   */
  for (size_t k = 0; k < n; ++k)
  {
    x[q[k]] = work[k];
  }
}

static void asd_la_backslash(size_t m, size_t n, long double *A,
                             long double *B, long double *C)
{
  long double *LU = (long double *)asd_malloc_aligned(m * m * sizeof(long double));
  long double *work = (long double *)asd_malloc_aligned(m * sizeof(long double));

  size_t *p = (size_t *)asd_malloc_aligned(m * sizeof(size_t));
  size_t *q = (size_t *)asd_malloc_aligned(m * sizeof(size_t));

  memcpy(LU, A, m * m * sizeof(long double));
  memcpy(C, B, m * n * sizeof(long double));

  asd_la_lu_factor(m, LU, p, q);

  for (size_t j = 0; j < n; ++j)
    asd_la_lu_solve(m, LU, C + j * m, p, q, work);

  asd_free_aligned(q);
  asd_free_aligned(p);
  asd_free_aligned(work);
  asd_free_aligned(LU);
}

static void asd_la_forwardslash(size_t m, size_t n, long double *A,
                                long double *B,
                                long double *C)
{
  long double *AT = (long double *)asd_malloc_aligned(n * m * sizeof(long double));
  long double *BT = (long double *)asd_malloc_aligned(n * n * sizeof(long double));
  long double *CT = (long double *)asd_malloc_aligned(n * m * sizeof(long double));

  asd_la_mtranspose(m, n, A, m, AT, n);
  asd_la_mtranspose(n, n, B, n, BT, n);

  asd_la_backslash(n, m, BT, AT, CT);

  asd_la_mtranspose(n, m, CT, n, C, m);

  asd_free_aligned(CT);
  asd_free_aligned(BT);
  asd_free_aligned(AT);
}
// }}}

// {{{ Jacobi
#define ASD_APPROX_EQ(x, y, K, abs, eps, min)                                  \
  ((abs)((x) - (y)) < (min) + (K) * (eps)*ASD_MAX((abs)((x)), (abs)((y))))

#define ASD_LONG_DOUBLE_APPROX_EQ(x, y, K)                                     \
  ASD_APPROX_EQ((x), (y), (K), fabsl, LDBL_EPSILON, LDBL_EPSILON *LDBL_EPSILON)

/*
 * This function computes the normalization of the Jacobi polynomial
 * $\left\{h_N^{(\alpha,\beta)}\right\}^{-\frac12}$ where [see @Szego39
 * (4.3.3)]
 *
 * $$
 *  \begin{aligned}
 *    h_N^{(\alpha,\beta)} &=
 *      \int_{-1}^{+1} (1-x)^{\alpha} (1+x)^{\beta}
 *    \left\{P_N^{(\alpha,\beta)} (x)\right\}^2 \, dx \\
 *    &=
 *    \frac{2^{\alpha+\beta+1}}{2N+\alpha+\beta+1}
 *    \frac{\Gamma(N+\alpha+1)\Gamma(N+\beta+1)}
 *         {\Gamma(N+\alpha+\beta+1)\Gamma(N+1)}.
 *  \end{aligned}
 * $$
 */
static long double asd_jacobi_h_inv_sqrt(long double alpha, long double beta,
                                         int N)
{
  ASD_ASSERT(N >= 0);
  ASD_ASSERT(alpha >= -1.0L);
  ASD_ASSERT(beta >= -1.0L);
  ASD_ASSERT(!(ASD_LONG_DOUBLE_APPROX_EQ(alpha, -0.5L, 10) &&
               ASD_LONG_DOUBLE_APPROX_EQ(beta, -0.5L, 10)));

  long double lgn = -(alpha + beta + 1) * logl(2) - lgammal(N + alpha + 1) -
                    lgammal(N + beta + 1) + logl(2 * N + alpha + beta + 1) +
                    lgammal(N + 1) + lgammal(N + alpha + beta + 1);
  return sqrtl(expl(lgn));
}

/*
 * This function evaluates the orthonormal polynomial $p_N(x)$ associated with
 * the Jacobi polynomial where $p_N(x) =
 * \left\{h_N^{(\alpha,\beta)}\right\}^{-\frac12} P_N^{(\alpha,\beta)} (x)$.
 *
 * The Jacobi polynomials are a set of polynomial functions for $\alpha > -1$,
 * $\beta > -1$, and $N$ a non-negative integer.  The functions are defined
 * on $[-1, +1]$ and orthogonal with respect to the weight function
 * $$w(x)=(1-x)^\alpha(1+x)^\beta.$$  Here we use the same normalization as
 * @Szego39, i.e., $P_N^{(\alpha,\beta)}(1) = {n+\alpha\choose n}$.  Thus we
 * have
 * $$
 *   \int_{-1}^{+1} p_N(x) p_M(x) w(x) \, dx = \delta_{NM}.
 * $$
 *
 * The three term recurrence relation arrived at by rearranging @Szego39
 * [(4.5.1)] is
 * $$
 *   P_n^{(\alpha,\beta)}(x) = (ax-b) P_{n-1}^{(\alpha,\beta)}(x)
 *                             - c P_{n-2}^{(\alpha,\beta)}
 * $$
 * where
 * $$
 * \begin{aligned}
 *   a &= \frac{(2n + \alpha + \beta -1)(2n + \alpha + \beta)}
 *             {2n (n + \alpha + \beta)} \\
 *   b &= \frac{(\beta^2 - \alpha^2)(2n + \alpha + \beta - 1)}
 *             {2n(n + \alpha + \beta)(2n + \alpha + \beta - 2)} \\
 *   c &= \frac{(n + \alpha - 1)(n + \beta - 1)(2n + \alpha + \beta)}
 *             {n(n + \alpha + \beta)(2n + \alpha + \beta - 2)}
 * \end{aligned}
 * $$
 * with $P_0^{(\alpha,\beta)}(x) = 1$ and
 * $P_1^{(\alpha,\beta)}(x) =  \frac12(\alpha + \beta + 2)x
 *                           + \frac12(\alpha - \beta)$.
 */
static void asd_jacobi_p(long double alpha, long double beta, int N, size_t nx,
                         long double *x, long double *P)
{
  ASD_ASSERT(N >= 0);
  ASD_ASSERT(alpha >= -1.0L);
  ASD_ASSERT(beta >= -1.0L);
  ASD_ASSERT(!(ASD_LONG_DOUBLE_APPROX_EQ(alpha, -0.5L, 10) &&
               ASD_LONG_DOUBLE_APPROX_EQ(beta, -0.5L, 10)));

  for (size_t i = 0; i < nx; ++i)
  {
    long double P_n_2;
    long double P_n_1 = 1;
    long double P_n_0 = ((alpha + beta + 2) / 2) * x[i] + (alpha - beta) / 2;
    if (N == 0)
    {
      P[i] = P_n_1;
    }
    else if (N == 1)
    {
      P[i] = P_n_0;
    }
    else
    {
      for (int n = 2; n < N + 1; ++n)
      {
        long double a = (2 * n + alpha + beta - 1) * (2 * n + alpha + beta) /
                        (2 * n * (n + alpha + beta));
        long double b =
            (beta * beta - alpha * alpha) * (2 * n + alpha + beta - 1) /
            (2 * n * (n + alpha + beta) * (2 * n + alpha + beta - 2));
        long double c = (n + alpha - 1) * (n + beta - 1) *
                        (2 * n + alpha + beta) /
                        (n * (n + alpha + beta) * (2 * n + alpha + beta - 2));

        P_n_2 = P_n_1;
        P_n_1 = P_n_0;
        P_n_0 = (a * x[i] - b) * P_n_1 - c * P_n_2;
      }
      P[i] = P_n_0;
    }
  }

  /*
   * Normalize the Jacobi polynomials
   */
  long double h_inv_sqrt = asd_jacobi_h_inv_sqrt(alpha, beta, N);
  for (size_t i = 0; i < nx; ++i)
    P[i] *= h_inv_sqrt;

  return;
}

/*
 * This function evaluates the derivative of the orthonormal polynomial
 * $p_N(x)$ associated with the Jacobi polynomial where
 * $p_N(x) = \left\{h_N^{(\alpha,\beta)}\right\}^{-\frac12}
 * P_N^{(\alpha,\beta)} (x)$.
 *
 * For the evaluation of the derivative we use the identity
 * $$
 *   \frac{d}{dx} P_N^{(\alpha,\beta)} (x) =
 *    \frac{N+\alpha+\beta+1}{2} P_{N-1}^{(\alpha+1,\beta+1)} (x)
 * $$
 * along with
 * $$
 *   h_N^{(\alpha,\beta)} =
 *     \frac{N+\alpha+\beta+1}{4N} h_{N-1}^{(\alpha+1,\beta+1)}
 * $$
 * to get
 * $$
 * \begin{aligned}
 *   \frac{d}{dx} p_N^{(\alpha,\beta)} (x)
 *     &= \left\{h_N^{(\alpha,\beta)}\right\}^{-\frac12}
 *        \frac{d}{dx} P_N^{(\alpha,\beta)} (x) \\
 *     &= \left\{h_N^{(\alpha,\beta)}\right\}^{-\frac12}
 *        \frac{N+\alpha+\beta+1}{2}
 *        P_{N-1}^{(\alpha+1,\beta+1)} (x) \\
 *     &= \left(\frac{4N}{N+\alpha+\beta+1}\right)^{\frac12}
 *        \frac{N+\alpha+\beta+1}{2}
 *        \left\{h_{N-1}^{(\alpha+1,\beta+1)}\right\}^{-\frac12}
 *        P_{N-1}^{(\alpha+1,\beta+1)} (x) \\
 *     &= \left(N(N+\alpha+\beta+1)\right)^{\frac12}
 *        p_{N-1}^{(\alpha+1,\beta+1)} (x).
 * \end{aligned}
 * $$
 */
static void asd_grad_jacobi_p(long double alpha, long double beta, int N,
                              size_t nx, long double *x, long double *dP)
{
  ASD_ASSERT(N >= 0);
  ASD_ASSERT(alpha >= -1.0L);
  ASD_ASSERT(beta >= -1.0L);

  if (N == 0)
  {
    for (size_t i = 0; i < nx; ++i)
    {
      dP[i] = 0.0L;
    }
  }
  else
  {
    asd_jacobi_p(alpha + 1, beta + 1, N - 1, nx, x, dP);
    long double scale = sqrtl(N * (N + alpha + beta + 1));
    for (size_t i = 0; i < nx; ++i)
    {
      dP[i] *= scale;
    }
  }

  return;
}

static void asd_jacobi_gauss_quadrature_half(long double alpha,
                                             long double beta, int N, int half,
                                             long double *x,
                                             long double *w)
{
  ASD_ASSERT(N >= 0);
  ASD_ASSERT(alpha >= -1.0L);
  ASD_ASSERT(beta >= -1.0L);
  ASD_ASSERT(!(ASD_LONG_DOUBLE_APPROX_EQ(alpha, -0.5L, 10) &&
               ASD_LONG_DOUBLE_APPROX_EQ(beta, -0.5L, 10)));

  const int MAX_ITERATIONS = 200;

  int nk = (half) ? (((N + 1) % 2) ? (N + 1) / 2 + 1
                                   : (N + 1) / 2) /* ceil((N + 1)/2) */
                  : (N + 1) / 2;                  /* floor((N + 1)/2) */

  if (nk == 0)
    return;

  long double tworho = 2 * (N + 1) + alpha + beta + 1;
  long double *tmp;

  long double *theta0 = (long double *)asd_malloc_aligned(nk * sizeof(long double));
  long double *theta1 = (long double *)asd_malloc_aligned(nk * sizeof(long double));
  long double *p0 = (long double *)asd_malloc_aligned(nk * sizeof(long double));
  long double *dp0 = (long double *)asd_malloc_aligned(nk * sizeof(long double));

  /*
   * Use Gatteschi and Pittaluga's approximation for the roots of the Jacobi
   * polynomials as an initial guess.  See equation (3.19) of
   * Nicholas Hale and Alex Townsend ``Fast and Accurate Computation of
   * Gauss–Legendre and Gauss–Jacobi Quadrature Nodes and Weights'' SIAM J.
   * SCI. COMPUT. Vol. 35, No. 2, pp. A652–A674.
   */
  for (int k = nk; k > 0; --k)
  {
    int khat = (half) ? nk - k : k - 1;

    long double pi = 4 * atanl(1);

    long double phik = (2 * k + alpha - 0.5L) * pi / tworho;

    theta1[khat] = phik +
                   1 / (tworho * tworho) *
                       ((0.25L - alpha * alpha) * 1 / tanl(0.5L * phik) -
                        (0.25L - beta * beta) * tanl(0.5L * phik));
  }

  /*
   * Use Newton's method for finding the roots of the Jacobi polynomial.
   */
  int converged = 0;
  for (int i = 0; i < MAX_ITERATIONS; ++i)
  {
    tmp = theta0;
    theta0 = theta1;
    theta1 = tmp;

    for (int k = 0; k < nk; ++k)
    {
      x[k] = cosl(theta0[k]);
    }

    asd_jacobi_p(alpha, beta, N + 1, nk, x, p0);
    asd_grad_jacobi_p(alpha, beta, N + 1, nk, x, dp0);

    for (int k = 0; k < nk; ++k)
    {
      theta1[k] = theta0[k] - p0[k] / (-sinl(theta0[k]) * dp0[k]);
    }

    int diff = 0;
    for (int k = 0; k < nk; ++k)
    {
      diff += !ASD_LONG_DOUBLE_APPROX_EQ(theta0[k], theta1[k], 10);
    }
    if (!diff)
    {
      converged = 1;
      break;
    }
  }

  ASD_ABORT_IF(
      !converged,
      "Newton's method does not converge when computing Jacobi Gauss points");
  /*
   * Nodes
   */
  for (int k = 0; k < nk; ++k)
  {
    x[k] = cosl(theta1[k]);
  }

  /*
   * Weights
   */
  asd_grad_jacobi_p(alpha, beta, N + 1, nk, x, dp0);

  for (int k = 0; k < nk; ++k)
  {
    long double sint = sinl(theta1[k]);
    w[k] = tworho / (sint * sint * dp0[k] * dp0[k]);
  }

  asd_free_aligned(theta0);
  asd_free_aligned(theta1);
  asd_free_aligned(p0);
  asd_free_aligned(dp0);

  return;
}

static void asd_jacobi_gauss_quadrature(long double alpha, long double beta,
                                        int N, long double *x,
                                        long double *w)
{
  int nk_floor = (N + 1) / 2; /* floor((N + 1)/2) */

  asd_jacobi_gauss_quadrature_half(alpha, beta, N, 1, x + nk_floor,
                                   w + nk_floor);
  asd_jacobi_gauss_quadrature_half(beta, alpha, N, 0, x, w);

  for (int k = 0; k < nk_floor; ++k)
    x[k] *= -1;

  return;
}

void asd_jacobi_gauss_lobatto_quadrature(long double alpha, long double beta,
                                         int N, long double *x,
                                         long double *w)
{
  ASD_ASSERT(N >= 1);

  x[0] = -1;
  x[N] = 1;

  if (N > 1)
  {
    asd_jacobi_gauss_quadrature(alpha + 1, beta + 1, N - 2, x + 1, w + 1);
  }

  asd_jacobi_p(alpha, beta, N, N + 1, x, w);
  long double fac = (2 * N + alpha + beta + 1) / (N * (N + alpha + beta + 1));
  for (int k = 0; k < N + 1; ++k)
  {
    w[k] = fac / (w[k] * w[k]);
  }

  w[0] *= (1 + beta);
  w[N] *= (1 + alpha);

  return;
}

void asd_jacobi_p_vandermonde(long double alpha, long double beta, int N,
                              size_t nx, long double *x, long double *V)
{
  for (int j = 0; j <= N; ++j)
    asd_jacobi_p(alpha, beta, j, nx, x, V + j * nx);

  return;
}

static void asd_grad_jacobi_p_vandermonde(long double alpha, long double beta,
                                          int N, size_t nx, long double *x,
                                          long double *V)
{
  for (int j = 0; j <= N; ++j)
    asd_grad_jacobi_p(alpha, beta, j, nx, x, V + j * nx);

  return;
}

void asd_jacobi_p_interpolation(long double alpha, long double beta, int N,
                                size_t nx, long double *x, long double *V,
                                long double *I)
{

  long double *Vx = (long double *)asd_malloc_aligned((N + 1) * nx * sizeof(long double));

  asd_jacobi_p_vandermonde(alpha, beta, N, nx, x, Vx);

  asd_la_forwardslash(nx, N + 1, Vx, V, I);

  asd_free_aligned(Vx);

  return;
}

void asd_jacobi_p_differentiation(long double alpha, long double beta, int N,
                                  size_t nx, long double *x, long double *V,
                                  long double *D)
{

  long double *Vx = (long double *)asd_malloc_aligned((N + 1) * nx * sizeof(long double));

  asd_grad_jacobi_p_vandermonde(alpha, beta, N, nx, x, Vx);

  asd_la_forwardslash(nx, N + 1, Vx, V, D);

  asd_free_aligned(Vx);

  return;
}

void asd_jacobi_p_mass(int N, long double *V, long double *M)
{
  long double *I = (long double *)asd_malloc_aligned((N + 1) * (N + 1) * sizeof(long double));

  long double *invV =
      (long double *)asd_malloc_aligned((N + 1) * (N + 1) * sizeof(long double));

  for (int i = 0; i < (N + 1) * (N + 1); ++i)
    I[i] = 0;

  for (int i = 0; i < (N + 1) * (N + 1); ++i)
    M[i] = 0;

  for (int i = 0; i <= N; ++i)
    I[(N + 1) * i + i] = 1;

  asd_la_backslash(N + 1, N + 1, V, I, invV);

  asd_la_mTmmult(N + 1, N + 1, N + 1, invV, N + 1, invV, N + 1, M, N + 1);

  asd_free_aligned(I);
  asd_free_aligned(invV);

  return;
}

void asd_jacobi_p_h_project(int N, long double h, long double *V,
                            long double *I, long double *M, long double *P)
{

  long double *ITM =
      (long double *)asd_malloc_aligned((N + 1) * (N + 1) * sizeof(long double));
  long double *VTITM =
      (long double *)asd_malloc_aligned((N + 1) * (N + 1) * sizeof(long double));

  for (int i = 0; i < (N + 1) * (N + 1); ++i)
    P[i] = 0;

  for (int i = 0; i < (N + 1) * (N + 1); ++i)
    ITM[i] = 0;

  for (int i = 0; i < (N + 1) * (N + 1); ++i)
    VTITM[i] = 0;

  asd_la_mTmmult(N + 1, N + 1, N + 1, I, N + 1, M, N + 1, ITM, N + 1);
  asd_la_mTmmult(N + 1, N + 1, N + 1, V, N + 1, ITM, N + 1, VTITM, N + 1);

  asd_la_mmmult(N + 1, N + 1, N + 1, V, N + 1, VTITM, N + 1, P, N + 1);

  for (int i = 0; i < (N + 1) * (N + 1); ++i)
    P[i] *= h;

  asd_free_aligned(ITM);
  asd_free_aligned(VTITM);

  return;
}
// }}}
