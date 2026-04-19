#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "bgzf.h"

#ifdef _WIN32
#define ftell64 _ftelli64
#define fseek64 _fseeki64
#else
#define ftell64 ftello
#define fseek64 fseeko
#endif

static int getFileSize(BGZF *bgzf)
{
    fseek64(bgzf->filePointer, 0, SEEK_END);
    bgzf->fileSize = ftell64(bgzf->filePointer);
    return 0;
}

static int readBGZFBlock(BGZF *bgzf)
{
    /*
    ID1: uint8, ID2: uint8, CM: uint8, FLG: uint8, MTIME: uint32, XFL: uint8, OS: uint8, [10B]
    XLEN: uint16 [12B]
    SI1: uint8, SI2: uint8, SLEN: uint16, [16B]
    BSIZE: uint16, [18B]
    X: XLEN-6 bytes, [XLEN-6 B]
    CDATA: uint8[BSIZE-XLEN-19]
    CRC32: uint32, ISIZE: uint32 [8B]
    */
    fread(bgzf->compressedData, 1, 18, bgzf->filePointer);
    uint16_t xLen = bgzf->compressedData[10] | (bgzf->compressedData[11] << 8);
    uint16_t bSize = bgzf->compressedData[16] | (bgzf->compressedData[17] << 8);
    fseek64(bgzf->filePointer, xLen - 6, SEEK_CUR); // x
    uint64_t cDataSize = bSize - xLen - 19;
    fread(bgzf->compressedData, 1, cDataSize + 8, bgzf->filePointer);
    bgzf->dataSize = bgzf->compressedData[cDataSize + 4] | (bgzf->compressedData[cDataSize + 5] << 8) | (bgzf->compressedData[cDataSize + 6] << 16) | (bgzf->compressedData[cDataSize + 7] << 24); // uint32_t
    bgzf->dataOffset = 0;
    bgzf->fileOffset = bgzf->FileOffset;
    bgzf->FileOffset += bSize + 1;
    z_stream stream;
    memset(&stream, 0, sizeof(z_stream));
    stream.next_in = bgzf->compressedData;
    stream.avail_in = cDataSize;
    stream.next_out = bgzf->decompressedData;
    stream.avail_out = bgzf->dataSize;
    inflateInit2(&stream, -15);
    int flag = inflate(&stream, Z_FINISH);
    inflateEnd(&stream);
    if (flag != Z_STREAM_END) fputs("An error has occurred while decompressing data.", stderr);
    return 0;
}

BGZF *bgzfOpen(char *file, uint64_t fileOffset)
{
    BGZF *bgzf = malloc(sizeof(BGZF));
    bgzf->filePointer = fopen(file, "rb");
    getFileSize(bgzf);
    bgzf->fileOffset = fileOffset;
    fseek64(bgzf->filePointer, bgzf->fileOffset, SEEK_SET);
    bgzf->FileOffset = bgzf->fileOffset;
    bgzf->compressedData = malloc(65536);
    bgzf->decompressedData = malloc(65536);
    bgzf->x = malloc(65536);
    bgzf->xSize = 65536;
    bgzf->dataSize = 0;
    bgzf->dataOffset = 0;
    return bgzf;
}

int bgzfRead(BGZF *bgzf, uint64_t n)
{
    if (n > bgzf->xSize)
    {
        bgzf->x = realloc(bgzf->x, n);
        bgzf->xSize = n;
    }
    int eof = 0;
    uint64_t N = 0;
    while (N < n && !eof)
    {
        if (bgzf->dataSize)
        {
            uint64_t x = n - N < bgzf->dataSize ? n - N : bgzf->dataSize;
            memcpy(bgzf->x + N, bgzf->decompressedData + bgzf->dataOffset, x);
            bgzf->dataOffset += x;
            bgzf->dataSize -= x;
            N += x;
        }
        else if (bgzf->FileOffset < bgzf->fileSize) readBGZFBlock(bgzf);
        else eof = 1;
    }
    return eof;
}

int bgzfClose(BGZF *bgzf)
{
    fclose(bgzf->filePointer);
    free(bgzf->compressedData);
    free(bgzf->decompressedData);
    free(bgzf->x);
    free(bgzf);
    return 0;
}
