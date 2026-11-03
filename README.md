# vcf_to_snp

Converts variant call format (VCF) data into a SNP-focused representation for downstream genetic analysis.

## VCF file format

The VCF files we are working with have a few sections:
 - Header section (lines starting with `##`)
 - Column header line (starting with `#CHROM`)
 - Data lines (one per SNP/variant)

### Header section

The header section contains metadata, which is not of major concern if we are only interested in extracting the SNPs.

### Column header line

The column header line defines the columns in the VCF file. The standard columns are:

| Column | Description |
|--------|-------------|
| CHROM  | Chromosome name |
| POS    | Position of the variant on the chromosome |
| ID     | Identifier for the variant (e.g., rs number) |
| REF    | Reference allele |
| ALT    | Alternative allele(s) |
| QUAL   | Quality score for the variant |
| FILTER | Filter status (e.g., PASS, LOWQUAL) |
| INFO   | Additional information about the variant |
| FORMAT | Format string for the genotype fields |
| \*samples\* | Variance information for each genome sample |

For the 1000 Genome Project dataset, there is 1 sample column for each individual genome, with the column header name being the genome ID (e.g. `HG00096`).

### Data lines

The variance information is formatted as `A/B`, where `A` and `B` are booleans representing whether the genome has a variant allele at that position. Since humans have two copies of each chromosome, each individual will have two alleles for each variant, represented by the two values in the `A/B` format. For example, `0/0` represents an individual with no variant alleles at that position, `0/1` represents a heterozygous individual with one variant allele, and `1/1` represents a homozygous individual with two variant alleles. It can also be `./.` representing missing genotype information.

The separator is not always `/` (sometimes `|` is used instead).