#include <stdio.h>
#include <stdlib.h>

int fastq2fasta(char argv[], char outfile[]);

int fastq2fasta(char argv[], char outfile[])
{
  /* 1 second to read without print and 6 seconds in total if we print */

  FILE *fp;
  char *line = NULL;
  FILE *fout;
  size_t len = 0;
  ssize_t read;

  /* int correspond to 4.9 billions. should be enough*/
  unsigned int count = 0;

  /* use argv*/

  fp = fopen(argv, "r");
  if (fp == NULL){
        printf("No such file\n");
        return 1;
  }

  fout = fopen(outfile, "w");
  if (fp == NULL){
        printf("No such file\n");
        return 1;
   }

  while ((read = getline(&line, &len, fp)) != -1) {
      /* skip empty lines (last one for example)*/
      if (read == 1UL){
        continue;
      }

      if ((count % 4UL == 0UL )){
          fprintf(fout, ">%s", &line[1]);
      }
      else if (count % 4UL == 1UL ){
          fprintf(fout, "%s", &line[0]);
      }
      count++;
  }

  fclose(fp);
  fclose(fout);

  if (line){
     free(line);
  }
  return 0;
}




/* use stdout and this can be an executable
 *
int fastq2fasta(argc, *argv[]) 
{
  FILE *fp;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;

  unsigned int count = 0;

  fp = fopen(argv[0], "r");
  if (fp == NULL){
        printf("No such file\n");
        return 1;
  }

  while ((read = getline(&line, &len, fp)) != -1) {
      if (read == 1UL){
        continue;
      }

      if ((count % 4UL == 0UL )){
          fputs(">", stdout);
          fputs(&line[1], stdout);
      }
      else if (count % 4UL == 1UL ){
          fputs(line, stdout);
      }
      count++;
  }

  fclose(fp);

  if (line){
     free(line);
  }
  return 0;
}
*/
