# bamIndexCreate

Create index of the sorted bam file.

bamIndexCreate is written in C and used by MetaCAT to build BAM indexes.

## Usage

```text
bamIndexCreate v1.1.0
Create index of the sorted bam file.
https://github.com/liu-congcong/MetaCAT

Usage:
  bamIndexCreate [options]

Options:
  -i    Path to the input sorted bam file
  -o    Path to the output index file
        If not specified, the index will be shown on stdout in a human-readable format
```

## Examples

```text
bamIndexCreate -i bam -o bam.index (binary)
bamIndexCreate -i bam (non-binary)
```
