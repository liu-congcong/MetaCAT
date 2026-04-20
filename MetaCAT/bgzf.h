#ifndef BGZF_H
#define BGZF_H

#include <stdint.h>
#include <stdio.h>

typedef struct
{
    FILE *filePointer;
    uint8_t *compressedData;
    uint8_t *decompressedData;
    uint8_t *x;
    uint64_t xSize;
    uint64_t fileSize;
    uint64_t fileOffset;
    uint64_t FileOffset;
    uint64_t dataSize;
    uint64_t dataOffset;
} BGZF;

BGZF *bgzfOpen(char *, uint64_t);
int bgzfRead(BGZF *, uint64_t);
int bgzfClose(BGZF *);

#endif