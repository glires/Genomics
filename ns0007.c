/*                                                                           */
/* NAME                                                                      */
/*   ns0007 - count reads and bases in a FASTQ file                          */
/*                                                                           */
/* SYNOPSIS                                                                  */
/*   ns0007 [min_qscore or -q min_qscore] [-h] [-v] input.fastq              */
/*                                                                           */
/* USAGE                                                                     */
/*   $ gcc -W -Wall -O -std=c99 -pedantic -Werror -o ns0007 ns0007.c         */
/*   $ ns0007 -h                                                             */
/*   $ ns0007 -v                                                             */
/*   $ ns0007 mouse_R1.fastq                                                 */
/*   $ ns0007 < mouse_R1.fastq                                               */
/*   $ ns0007 20 mouse_R1.fastq                                              */
/*   $ ns0007 -q 21 mouse_R1.fastq    # recommended                          */
/*   $ ns0007 -q 22 mouse_R1.fastq | awk '{ print $3 }'                      */
/*   $ ns0007 -q 23 < mouse_R1.fastq                                         */
/*   $ head -n 40000 mouse_R1.fastq | ns0007 -q24                            */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   This programme counts numbers of reads and bases in a FASTQ file.       */
/*   The minimal score can be specified as the first command argument.       */
/*   The default value for the  minimal score is 16.                         */
/*   Bases that have no less than the minimal socre are counted up.          */
/*   As a result, three numbers are printed onto the standard output.        */
/*   They are numbers of reads, higher-scored bases, and all bases,          */
/*   respectively.                                                           */
/*   Since this programme uses the unsigned long long type which is not      */
/*   defined in ANSI C, option -ansi cannot be specified in gcc.             */
/*   Instead, option -std=c99 can be specified.                              */
/*                                                                           */
/* OPTIONS                                                                   */
/*   -h  print simple usage                                                  */
/*   -q  minimal quality score (default: 16)                                 */
/*         However, the first command argument can be the minimal quality    */
/*         score without using -q                                            */
/*   -v  print programme name and its version                                */
/*                                                                           */
/* AUTHOR                                                                    */
/*   Coded by Kohji OKAMURA, Ph.D.                                           */
/*                                                                           */
/* HISTORY                                                                   */
/*   Oct 01, 2012  Copied from count_bases.c coded on Aug 06, 2012           */
/*   Apr 19, 2013  Use the return value of fgets to avoid warnings           */
/*   Nov 07, 2013  Rename from jf0007.c to ns0007.c                          */
/*   Jun 14, 2016  Process integer in the second command argument            */
/*   Jun 15, 2016  Employ lseek() to check pipe                              */
/*   Jun 16, 2016  Execute fclose() only when it is opened                   */
/*   Jun 17, 2016  Release as ns0007 ver. 1.3                                */
/*                                                                           */


#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define VERSION "1.3"
#define MAX_CHAR 4096
#define MIN_QUAL_SCORE 16

int getopt(int, char * const [], const char *);

extern char *optarg;


void fgets_error(int error_code)
{
  fprintf(stderr, "Unexpected error in fgets(): %d\n", error_code);
  exit(EXIT_FAILURE);
}


int is_small_int(char *first_arg)
{ /* return 1 if 0-99 */
  int length;

  length = (int)strlen(first_arg);
  if (length == 1)
  { if (isdigit(first_arg[0])) return 1; else return 0; }
  else if (length == 2)
  {
    if (isdigit(first_arg[0]) && isdigit(first_arg[1])) return 1;
    else return 0;
  }
  return 0;	/* more than 2 */
}


int main(int argc, char *argv[])
{
  int i, l, opt;
  unsigned short minimal_score = MIN_QUAL_SCORE;
  unsigned long long int count_reads = 0;
  unsigned long long int count_bases = 0;
  unsigned long long int count_all   = 0;
  char line[MAX_CHAR];
  FILE *check;

  if (argc == 1) goto read_data;
  if (is_small_int(argv[1]))
  { minimal_score = (unsigned short)strtol(argv[1], (char **)NULL, 10); }
  else
  {
    while ((opt = getopt(argc, argv, "hq:v")) != -1)
    {
      switch (opt)
      {
        case 'h': fprintf(stdout, "Usage: ns0007 [min_qscore or -q min_qscore "
                    "(default: %d)] input.fastq\n"
                    "Output: number of reads, higher-scored bases, "
                    "and all bases\n", MIN_QUAL_SCORE);
                  return EXIT_SUCCESS;
        case 'q': minimal_score
                    = (unsigned short)strtol(optarg, (char **)NULL, 10);
                  break;
        case 'v': fprintf(stdout, "ns0007 ver. %s\n", VERSION);
                  return EXIT_SUCCESS;
        default:  fprintf(stderr, "Unknown option: %c\n", opt);
                    /* does not exit */
      }
    }
  }

  check = fopen(argv[argc - 1], "r");
  if (check != NULL) { freopen(argv[argc - 1], "r", stdin); fclose(check); }
	/* Segmentation fault occurs when fclose an unopened file pointer */

  read_data: while (fgets(line, MAX_CHAR, stdin) != NULL)	/* 1st line */
  {
    if (fgets(line, MAX_CHAR, stdin) == NULL)	/* 2nd line */
    { fgets_error(2); }
    if (fgets(line, MAX_CHAR, stdin) == NULL)	/* 3rd line */
    { fgets_error(3); }
    if (fgets(line, MAX_CHAR, stdin) == NULL)	/* 4th line */
    { fgets_error(4); }
    l = strlen(line) - 1;	/* minus one for the new line */
    count_all += l;
    for (i = 0; i < l; i++)
    { if (((unsigned)line[i] - 33) >= minimal_score) count_bases++; }
    count_reads++;
  }

 fprintf(stdout, "%llu\t%llu\t%llu\n", count_reads, count_bases, count_all);

 return EXIT_SUCCESS;
}
