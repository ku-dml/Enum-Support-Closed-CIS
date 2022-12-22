/********** mylib.cpp **********/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "mylib.h"


/***** split a string into strings *****/
char **split( char *s, char c ){
  typedef struct list{
    int i;
    struct list *p;
  } List;
  char **str;
  List head,*tail,*ip;
  int countnum=0,counts=0;
  int i=0,j=0,i0=0;
  head.i=-1;
  head.p=NULL;
  tail=&head;
  while(s[i]!='\0'){
    if(s[i]==c){
      ip=(List*)mallocE(sizeof(List));
      ip->i=countnum;
      ip->p=NULL;
      tail->p=ip;
      tail=ip;
      countnum=-1;
      counts++;
    }
    countnum++;
    i++;
  }
 ip=(List*)mallocE(sizeof(List));
  ip->i=countnum;
  ip->p=NULL;
  tail->p=ip;
  tail=ip;
  str=(char**)mallocE((counts+1)*sizeof(char*));
  ip=&head;
  do{
    int k=0;
    ip=ip->p;
    str[j]=(char*)mallocE(((ip->i)+1)*sizeof(char));
    for(i=i0;i<i0+ip->i;i++){
      str[j][k]=s[i];
      k++;
    }
    str[j][k]='\0';
    i0+=ip->i+1;
    j++;
  }while(ip->p!=NULL);
  // free memory
  tail = head.p;
  while( tail != NULL ){
    ip = tail->p;
    free( tail );
    tail = ip;
  }
  return str;

}

/***** number of char c in string s *****/
int getNumChar( char *s, char c ){
  int numChar = 0;
  int k = 0;
  while( s[k] != '\0' ){
    if( s[k] == c )
      numChar++;
    k++;
  }
  return numChar;
}


/***** functions for memory allocation *****/
void *mallocE( size_t size ){
  void *s;
  if ( (s=malloc(size)) == NULL ) {
    fprintf( stderr, "malloc: not enough memory.\n" );
    exit( EXIT_FAILURE );
  }
  return s;
}

FILE *open_file(char *filename, char *mode){
  FILE *fp;
  fp = fopen(filename, mode);
  if(fp==NULL){
    fprintf(stderr,"failed to open %s\n",filename);
    exit(EXIT_FAILURE);
  }
  return fp;
}
