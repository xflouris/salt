/*
    Copyright (C) 2014 Tomas Flouri & Lucas Czech

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

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

/* TODO: separate it into xstrdup and xstrdup_aligned */
void * xstrdup(char * s)
{
  int len = strlen(s);
  //char * x = xmalloc(len+1, 8);
  char * x = xmalloc(roundup(len,16)+1, 16);
  strcpy(x,s);

  /* zero out remaining bytes */
  for (long i = len; i < roundup(len,16); ++i)
    x[i] = 0;

  return x;
}

char* strdup(const char *str)
{
    int n = strlen(str) + 1;
    char *dup = malloc(n * sizeof(char));
    if(dup)
    {
        strcpy(dup, str);
    }
    return dup;
}
