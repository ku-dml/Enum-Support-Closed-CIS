/********** mylib.h **********/

char **split(char *s, char c);
int getNumChar(char *s, char c);
void *mallocE(size_t size);
FILE *open_file(char *filename, char *mode);
int getVirtualMem();
int getPyisicalMem();
