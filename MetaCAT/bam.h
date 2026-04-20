#ifndef BAM_H
#define BAM_H

#include <stdint.h>

typedef struct
{
    uint32_t n;
    uint64_t headerFileOffset;
    uint64_t headerDataOffset;
    char **sequences;
    uint32_t *lengths;
    uint64_t *fileOffsets;
    uint64_t *dataOffsets;
    uint64_t *dataSizes;
} BAM;

typedef struct
{
    int32_t n;
    uint64_t headerFileOffset;
    uint64_t headerDataOffset;
    char **sequences;
    int32_t *lengths;
    uint64_t *fileOffsets;
    uint64_t *dataOffsets;
    uint64_t *dataSizes;
} INDEX;

int indexBam(char *, char *);

INDEX *readIndex(char *);

int freeIndex(INDEX *);

#endif