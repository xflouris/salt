#include "salt.h"

long gcd(long a, long b)
{
  if (b == 0)
  {
    return a;
  }
  else
  {
    return gcd(b, a % b);
  }
}

long roundup(long offset, long step)
{
  return ((offset - 1) | (step - 1)) + 1;
}

void fatal(const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf(stderr, format, argptr);
  va_end(argptr);
  fprintf(stderr, "\n");
  exit(1);
}

void * xmalloc(size_t size, const size_t alignment)
{
  void * t;
  posix_memalign(& t, alignment, size);

  if (t==NULL)
    fatal("Unable to allocate enough memory.");

  return t;
}

void * xrealloc(void *ptr, size_t size)
{
  void * t = realloc(ptr, size);
  if (!t)
    fatal("Unable to allocate enough memory.");
  return t;
}

void xfree (void* ptr)
{
    free(ptr);
}

char * xstrchrnul(char *s, int c)
{
  char * r = strchr(s, c);

  if (r)
    return r;
  else
    return (char *)s + strlen(s);
}

void * xstrdup(char * s)
{
  int len = strlen(s);
  char * x = xmalloc(len+1, 8);
  strcpy(x,s);
  x[len]=0;
  return x;
}
