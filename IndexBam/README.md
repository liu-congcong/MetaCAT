# IndexBam

Create index of the sorted bam file.

indexBam is written in C and used by MetaCAT to build BAM indexes.

## Usage

```TEXT
indexBam v1.0.0
Create index of the sorted bam file.
https://github.com/liu-congcong/MetaCAT

Usage:
  indexBam [options]

Options:
  -i    Path to the input sorted bam file
  -o    Path to the output index file
        If not specified, the index will be shown on stdout in a human-readable format
```

## Examples

```TEXT
indexBam -i bam -o bam.index (binary)
indexBam -i bam (non-binary)
```
