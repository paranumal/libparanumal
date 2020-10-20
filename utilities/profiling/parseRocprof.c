//Index,KernelName,gpu-id,queue-id,queue-index,pid,tid,grd,wgr,lds,scr,vgpr,sgpr,fbar,sig,WriteSize,DispatchNs,BeginNs,EndNs,CompleteNs
//0,"_occa_weightedNorm2_0.kd",0,1,0,16433,16439,262144,512,4096,0,16,32,0,0x7fa90163e700,4,4556288436792544,4556288440165491,4556288440211571,4556288440226587

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define mymin(a,b) ((a)<(b))?(a):(b)
#define mymax(a,b) ((a)>(b))?(a):(b)

typedef struct{
  double totalT;
  double aveT;
  char *name;
  double aveBW;
}entry_t;

typedef struct{
  char *str;
  int Nentries;
  char **entries;
} row_t;

typedef struct{

  int Nrows;

  row_t *rows;

}db_t;

void printRow(row_t *row){
  printf("test: %s, kernel: %s\n", row->entries[0], row->entries[1]);
}

row_t readRow(char *bufin){

  char buf[BUFSIZ];
  strcpy(buf, bufin);
  
  row_t *row = (row_t*) calloc(1, sizeof(row_t));

  row->str = strdup(bufin);
  
  int cntLabels = 0;
  char *dummy = strtok(buf, ",");
  ++cntLabels;
  
  while(strtok(NULL, ",")){
    ++cntLabels;
  }
  
  row->Nentries = cntLabels;
  
  // create database
  row->entries = (char**) calloc(cntLabels, sizeof(char*));

  // revert to original labels string
  strcpy(buf, bufin);
  
  int cnt = 0;
  row->entries[cnt] = strdup(strtok(buf, ","));
  ++cnt;
  
  for(int n=1;n<cntLabels;++n){
    row->entries[n] = strdup(strtok(NULL, ","));
  }

  return *row;
}

int compareRows(const void *a, const void *b){

  row_t *r1 = (row_t*) a;
  row_t *r2 = (row_t*) b;

  return strcmp(r1->entries[1], r2->entries[1]);
}

int compareTimes(const void *a, const void *b){

  entry_t *r1 = (entry_t*) a;
  entry_t *r2 = (entry_t*) b;

  return (r1->totalT < r2->totalT);
}


int main(int argc, char **argv){

  FILE *fp = fopen(argv[1], "r");

  char buf[BUFSIZ];
  
  // count lines
  int Nrows = 0;
  char *str = NULL;
  while(fgets(buf, BUFSIZ, fp))
    ++Nrows;

  --Nrows;
  
  fclose(fp);

  printf("Nrows=%d\n", Nrows);
  
  // reopen to read rows
  fp = fopen(argv[1], "r");
  
  db_t *db = (db_t*) calloc(1, sizeof(db_t));

  db->Nrows =  Nrows;
  db->rows = (row_t*) calloc(Nrows, sizeof(row_t));

  for(int r=0;r<Nrows;++r){
    fgets(buf, BUFSIZ, fp);
    db->rows[r] = readRow(buf);
  }

  qsort(db->rows+1, db->Nrows-1, sizeof(row_t), compareRows);

  int maxNentries = 1000;
  entry_t *entries = (entry_t*) calloc(maxNentries, sizeof(entry_t));
  
  //  printf("%40.40s, aveBW, minBW, maxBW, launches\n", "Kernel Name");
  //  printf("%40.40s, aveBW\n", "Kernel Name");
  int r = 1;
  int entry = 0;
  while(r<Nrows){
    long long int bytes = 0, T = 0, startT = 0, endT = 0;
    double sumBW = 0;
    int rbase = r, cnt = 0;
    double minBW = 1e20, maxBW = 0;
    while(strcmp(db->rows[rbase].entries[1],db->rows[r].entries[1])==0){
      long long int newbytes = 0;
      sscanf(db->rows[r].entries[15], "%lld", &newbytes);
      sscanf(db->rows[r].entries[17], "%lld", &startT);
      sscanf(db->rows[r].entries[18], "%lld", &endT);
      
      long long int  newT = endT-startT;

      T+=newT;
      
      // why the factor of 1024 ? (KB)
      double newBW = 1024.*newbytes/(double)newT;

      sumBW += newBW;

      minBW = mymin(minBW, newBW);
      maxBW = mymax(maxBW, newBW);
      
      //      printf("newBW=%lf\n", newBW);
      ++r;
      ++cnt;
      if(r>=Nrows) break;
    }

    double aveBW = sumBW/cnt;
    char name[BUFSIZ];
    sscanf(db->rows[rbase].entries[1], "\"_occa_%s_", name);
    char *basename = strtok(name, "_");
    
    //    printf("%40.40s, %5.4f, %5.4f, %5.4f, %d\n", name, aveBW, minBW, maxBW, cnt);
    //    printf("%40.40s, %5.4f, %5.4e\n", name, aveBW, T/1.e9);
    entries[entry].name = strdup(name);
    entries[entry].aveT = T/cnt;
    entries[entry].totalT = T;
    entries[entry].aveBW = aveBW;
    ++entry;
    
  }

  qsort(entries, entry, sizeof(entry_t), compareTimes);

  printf("%40.40s, aveBW, average T, total T\n", "Kernel Name");
  for(int e=0;e<entry;++e){
    
    printf("%40.40s, %5.4f, %5.4e, %5.4e\n", entries[e].name, entries[e].aveBW, entries[e].aveT/1.e9,
	   entries[e].totalT/1.e9);


  }
  
  fclose(fp);
  
}
