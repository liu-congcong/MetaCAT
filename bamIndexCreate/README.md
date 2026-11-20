# bamIndexCreate

Create index of the sorted bam file.

bamIndexCreate is written in C and used by MetaCAT to build BAM indexes.

## Usage

```text
bamIndexCreate v1.5.0
Create index of the sorted bam file.
https://github.com/liu-congcong/MetaCAT

Usage:
  bamIndexCreate [options]

Options:
  -i    Path to the input sorted bam files
  -t    Number of threads (default: 1)
```

## Examples

```text
bamIndexCreate -t 10 -i *.bam
```
