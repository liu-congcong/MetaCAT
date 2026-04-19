#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "bam.h"
#include "bgzf.h"

#ifdef _WIN32
#define fseek64 _fseeki64
#else
#define fseek64 fseeko
#endif

static int isValidBam(char *file)
{
    uint8_t x[30] = {
        0x1f, 0x8b,
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00,
        0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
        0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
    };
    uint8_t buffer[30];
    FILE *fp = fopen(file, "rb");
    fread(buffer, 1, 2, fp);
    fseek64(fp, -28, SEEK_END);
    fread(buffer + 2, 1, 28, fp);
    fclose(fp);
    return memcmp(buffer, x, 30);
}

static int isSortedBam(char *x)
{
    int flag = 1;
    char *p = strchr(x, '\n');
    if (p) flag = (strstr(x, "coordinate") == NULL);
    return flag;
}

static int readHeader(BGZF *bgzf, BAM *bam)
{
    uint8_t bamFlag[4] = {0x42, 0x41, 0x4D, 0x01};
    bgzfRead(bgzf, 8);
    if (memcmp(bgzf->x, bamFlag, 4)) return 1; /* magic code */
    uint32_t headerSize = bgzf->x[4] | (bgzf->x[5] << 8) | (bgzf->x[6] << 16) | (bgzf->x[7] << 24);
    bgzfRead(bgzf, headerSize + 4);
    if (isSortedBam((char *)bgzf->x)) return 2; /* sorted bam */
    bam->n = bgzf->x[headerSize] | (bgzf->x[headerSize + 1] << 8) | (bgzf->x[headerSize + 2] << 16) | (bgzf->x[headerSize + 3] << 24); // n + 1 //
    if (!bam->n) return 3; /* number of sequences > 0 */
    bam->sequences = malloc((bam->n + 1) * sizeof(char *));
    bam->lengths = malloc((bam->n + 1) * sizeof(uint32_t));
    for (uint32_t i = 0; i < bam->n; i++)
    {
        bgzfRead(bgzf, 4);
        uint32_t sequenceNameSize = bgzf->x[0] | (bgzf->x[1] << 8) | (bgzf->x[2] << 16) | (bgzf->x[3] << 24);
        bgzfRead(bgzf, sequenceNameSize + 4);
        bam->sequences[i] = malloc(sequenceNameSize * sizeof(char));
        memcpy(bam->sequences[i], bgzf->x, sequenceNameSize);
        bam->lengths[i] = bgzf->x[sequenceNameSize] | (bgzf->x[sequenceNameSize + 1] << 8) | (bgzf->x[sequenceNameSize + 2] << 16) | (bgzf->x[sequenceNameSize + 3] << 24);
    }
    bam->headerFileOffset = bgzf->fileOffset;
    bam->headerDataOffset = bgzf->dataOffset;
    bam->sequences[bam->n] = malloc(sizeof(char));
    memset(bam->sequences[bam->n], 0, 1);
    bam->lengths[bam->n] = 0;
    return 0;
}

static int readBody(BGZF *bgzf, BAM *bam)
{
    uint32_t alignmentSize;
    int32_t sequence;
    uint64_t fileOffset = bgzf->fileOffset;
    uint64_t dataOffset = bgzf->dataOffset;
    int32_t temp = -1;

    while (!bgzfRead(bgzf, 8))
    {
        alignmentSize = bgzf->x[0] | (bgzf->x[1] << 8) | (bgzf->x[2] << 16) | (bgzf->x[3] << 24);
        sequence = bgzf->x[4] | (bgzf->x[5] << 8) | (bgzf->x[6] << 16) | (bgzf->x[7] << 24);
        if (sequence == -1) sequence = bam->n;
        bgzfRead(bgzf, alignmentSize - 4);
        if (sequence != temp)
        {
            bam->fileOffsets[sequence] = fileOffset;
            bam->dataOffsets[sequence] = dataOffset;
            bam->dataSizes[sequence] = 0;
            temp = sequence;
        }
        fileOffset = bgzf->fileOffset;
        dataOffset = bgzf->dataOffset;
        bam->dataSizes[sequence] += 4 + alignmentSize;
    }
    return 0;
}

static int writeIndex(BAM *bam, char *file)
{
    bam->n++;
    gzFile fp = gzopen(file, "wb9");
    gzwrite(fp, &bam->n, sizeof(uint32_t));
    gzwrite(fp, &bam->headerFileOffset, sizeof(uint64_t));
    gzwrite(fp, &bam->headerDataOffset, sizeof(uint64_t));
    for (uint32_t i = 0; i < bam->n; i++)
    {
        //printf("%s\t%u\t%lu\t%lu\t%lu\n", bam->sequences[i], bam->lengths[i], bam->fileOffsets[i], bam->dataOffsets[i], bam->dataSizes[i]);
        uint32_t sequenceNameSize = strlen(bam->sequences[i]);
        gzwrite(fp, &sequenceNameSize, sizeof(uint32_t));
        gzwrite(fp, bam->sequences[i], sequenceNameSize);
        gzwrite(fp, bam->lengths + i, sizeof(uint32_t));
        gzwrite(fp, bam->fileOffsets + i, sizeof(uint64_t));
        gzwrite(fp, bam->dataOffsets + i, sizeof(uint64_t));
        gzwrite(fp, bam->dataSizes + i, sizeof(uint64_t));
        free(bam->sequences[i]);
    }
    gzclose(fp);

    free(bam->sequences);
    free(bam->lengths);
    free(bam->fileOffsets);
    free(bam->dataOffsets);
    free(bam->dataSizes);
    free(bam);
    return 0;
}

int indexBam(char *input, char *output)
{
    if (access(input, R_OK))
    {
        fprintf(stderr, "\"%s\" is not accessible.\n", input);
        return 1;
    }
    if (isValidBam(input))
    {
        fprintf(stderr, "\"%s\" is not a BGZF compressed BAM file.\n", input);
        return 2;
    }
    BAM *bam = malloc(sizeof(BAM));
    memset(bam, 0, sizeof(BAM));
    BGZF *bgzf = bgzfOpen(input, 0);
    int flag = readHeader(bgzf, bam);
    if (!flag)
    {
        bam->fileOffsets = malloc((bam->n + 1) * sizeof(uint64_t));
        memset(bam->fileOffsets, 0, (bam->n + 1) * sizeof(uint64_t));
        bam->dataOffsets = malloc((bam->n + 1) * sizeof(uint64_t));
        memset(bam->dataOffsets, 0, (bam->n + 1) * sizeof(uint64_t));
        bam->dataSizes = malloc((bam->n + 1) * sizeof(uint64_t));
        memset(bam->dataSizes, 0, (bam->n + 1) * sizeof(uint64_t));
        readBody(bgzf, bam);
        writeIndex(bam, output);
        flag -= 2;
    }
    else if (flag == 1) fprintf(stderr, "\"%s\" is not a BAM file.\n", input);
    else if (flag == 2) fprintf(stderr, "\"%s\" is not a sorted BAM file.\n", input);
    else if (flag == 3) fprintf(stderr, "\"%s\" contains no sequences.\n", input);
    bgzfClose(bgzf);
    flag += 2;
    return flag;
}

INDEX *readIndex(char *file)
{
    gzFile fp = gzopen(file, "rb");
    INDEX *index = malloc(sizeof(INDEX));
    gzread(fp, &index->n, 4);
    gzread(fp, &index->headerFileOffset, 8);
    gzread(fp, &index->headerDataOffset, 8);
    index->sequences = malloc(sizeof(char *) * index->n);
    index->lengths = malloc(sizeof(int32_t) * index->n);
    index->fileOffsets = malloc(sizeof(uint64_t) * index->n);
    index->dataOffsets = malloc(sizeof(uint64_t) * index->n);
    index->dataSizes = malloc(sizeof(uint64_t) * index->n);

    int32_t sequenceIDLength;
    for (int32_t i = 0; i < index->n; i++)
    {
        gzread(fp, &sequenceIDLength, 4);
        index->sequences[i] = malloc(sequenceIDLength + 1);
        gzread(fp, index->sequences[i], sequenceIDLength);
        index->sequences[i][sequenceIDLength] = 0;
        gzread(fp, index->lengths + i, 4);
        gzread(fp, index->fileOffsets + i, 8);
        gzread(fp, index->dataOffsets + i, 8);
        gzread(fp, index->dataSizes + i, 8);
    }
    gzclose(fp);
    return index;
}

int freeIndex(INDEX *index)
{
    for (int32_t i = 0; i < index->n; i++) free(index->sequences[i]);
    free(index->sequences);
    free(index->lengths);
    free(index->fileOffsets);
    free(index->dataOffsets);
    free(index->dataSizes);
    free(index);
    return 0;
}

// static int printIndex(char *file)
// {
//     if (access(file, R_OK))
//     {
//         fprintf(stderr, "\"%s\" is not accessible.\n", file);
//         return 1;
//     }
//     if (isValidBam(file))
//     {
//         fprintf(stderr, "\"%s\" is not a BGZF compressed BAM file.\n", file);
//         return 2;
//     }
//     BAM *bam = malloc(sizeof(BAM));
//     memset(bam, 0, sizeof(BAM));
//     BGZF *bgzf = bgzfOpen(file, 0);
//     int flag = readHeader(bgzf, bam);
//     if (!flag)
//     {
//         bam->fileOffsets = malloc((bam->n + 1) * sizeof(uint64_t));
//         memset(bam->fileOffsets, 0, (bam->n + 1) * sizeof(uint64_t));
//         bam->dataOffsets = malloc((bam->n + 1) * sizeof(uint64_t));
//         memset(bam->dataOffsets, 0, (bam->n + 1) * sizeof(uint64_t));
//         bam->dataSizes = malloc((bam->n + 1) * sizeof(uint64_t));
//         memset(bam->dataSizes, 0, (bam->n + 1) * sizeof(uint64_t));
//         readBody(bgzf, bam);
//         bam->n++;
//         puts("Sequence\tLength\tFileOffset\tDataOffset\tDataSize");
//         for (uint32_t i = 0; i < bam->n; i++)
//         {
//             printf("%s\t%u\t%llu\t%llu\t%llu\n", bam->sequences[i], bam->lengths[i], bam->fileOffsets[i], bam->dataOffsets[i], bam->dataSizes[i]);
//             free(bam->sequences[i]);
//         }
//         free(bam->sequences);
//         free(bam->lengths);
//         free(bam->fileOffsets);
//         free(bam->dataOffsets);
//         free(bam->dataSizes);
//         free(bam);
//     }
//     else if (flag == 1) fprintf(stderr, "\"%s\" is not a BAM file.\n", file);
//     else if (flag == 2) fprintf(stderr, "\"%s\" is not a sorted BAM file.\n", file);
//     else if (flag == 3) fprintf(stderr, "\"%s\" contains no sequences.\n", file);
//     bgzfClose(bgzf);
//     return flag;
// }

// static int printHelp()
// {
//     puts("indexBam v1.0.0");
//     puts("Create index of the sorted bam file.");
//     puts("https://github.com/liu-congcong/MetaCAT");
//     puts("\nUsage:");
//     puts("  indexBam [options]");
//     puts("\nOptions:");
//     puts("  -i    Path to the input sorted bam file");
//     puts("  -o    Path to the output index file");
//     puts("        If not specified, the index will be shown on stdout in a human-readable format");
//     return 0;
// }

// int main(int argc, char *argv[])
// {
//     char *inputFile = NULL;
//     char *outputFile = NULL;
//     for (int i = 1; i < argc; i++)
//     {
//         if (!strncmp(argv[i], "-i", 2) && i + 1 < argc) inputFile = argv[i + 1];
//         else if (!strncmp(argv[i], "-o", 2) && i + 1 < argc) outputFile = argv[i + 1];
//     }
//     if (!inputFile)
//     {
//         printHelp();
//         exit(EXIT_FAILURE);
//     }
//     if (outputFile) indexBam(inputFile, outputFile);
//     else printIndex(inputFile);
//     return 0;
// }