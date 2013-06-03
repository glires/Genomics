/*                                                                           */
/* NAME                                                                      */
/*   jf0007 - count high-scored bases in FASTQ                               */
/*                                                                           */
/* SYNOPSIS                                                                  */
/*   $ gcc -O -ansi -o jf0007 jf0007.c                                       */
/*   $ jf0007 < mouse.fastq                                                  */
/*   $ jf0007 30 < mouse.fastq                                               */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   This programme counts high-scored bases in a FASTQ format.              */
/*   The minimal score is specified as the first command argument.           */
/*   Default minimal score is 20.                                            */
/*   Bases that have no less than the minimal socre are counted up.          */
/*   Number of high-quality bases with number of entries are finally         */
/*   printed onto the standard output.                                       */
/*                                                                           */
/* OPTIONS                                                                   */
/*   No options are supported.                                               */
/*   The first command option can be the minimal quality score.              */
/*                                                                           */
/* AUTHOR                                                                    */
/*   Coded by Kohji OKAMURA, Ph.D.                                           */
/*                                                                           */
/* HISTORY                                                                   */
/*   Oct 01, 2012  Copied from count_bases.c coded on Aug 06, 2012           */
/*   Apr 19, 2013  Use the return value of fgets to avoid warnings           */
/*                                                                           */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_CHAR 2048
#define MIN_QUAL_SCORE 20


void fgets_error(int error_code)
{
  fprintf(stderr, "Unexpected error in fgets: %d\n", error_code);
  exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
  int i, l;
  unsigned short minimal_score;
  unsigned long long int count_base = 0;
  unsigned long long int count_line = 0;
  char line[MAX_CHAR];

  if (argc > 1)
  { minimal_score = (unsigned short)strtol(argv[1], (char **)NULL, 10); }
  else { minimal_score = MIN_QUAL_SCORE; }

  while (fgets(line, MAX_CHAR, stdin) != NULL)  /* 1st line */
  {
    if (fgets(line, MAX_CHAR, stdin) == NULL)   /* 2nd line */
    { fgets_error(2); }
    if (fgets(line, MAX_CHAR, stdin) == NULL)   /* 3rd line */
    { fgets_error(3); }
    if (fgets(line, MAX_CHAR, stdin) == NULL)   /* 4th line */
    { fgets_error(4); }
    l = strlen(line) - 1;       /* minus one for the new line */

    for (i = 0; i < l; i++)
    {
      if (((unsigned)line[i] - 33) >= minimal_score)
      { count_base++; }
    }
    count_line++;
  }

 fprintf(stdout, "%llu bases in %llu sequences\n", count_base, count_line);

 return EXIT_SUCCESS;
}
