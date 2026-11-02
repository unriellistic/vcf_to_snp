import gzip

import numpy as np


# SNP class adapted from https://github.com/diazale/1KGP_dimred/blob/master/Genotype_dimred_demo.ipynb
# with permission from Alex Diaz-Papkovich
class SNP:

    # memory optimization with __slots__
    __slots__ = ("failed", "failure_cause", "chrom", "position", "rid", "ref_allele", "alt_allele", "genotype_array")

    def __init__(self, line: str, autosome_only = True):
        """The initialization method takes in a line from the vcf file, as a string, 
        and records the relevant information. 
        line: a string from a vcf file
        """ 
        
        split_line = line.split()  #  First break down the line into a list of each field
        
        self.failed = False  # A label that we will set to True if something goes wrong.
        self.failure_cause = ""
        
        if line.startswith('#'):
            self.failed = True
            self.failure_cause = "line was a header line, not a snp"
            return
        
        if len(split_line)<=5:
            self.failed = True
            self.failure_cause = "incorrectly formatted line, should have at least 5 fields " + line
            return
          
        self.chrom = split_line[0]
        if autosome_only:
            if self.chrom not in ["%d" % (i,) for i in range(1,23)]:
                self.failed = True
                self.failure_cause = "not recognized as an autosome while autosome_only set to True"
                return
        
        self.chrom = int(split_line[0]) # Chromosome (numbered)
        self.position = int(split_line[1])  # The coordinates of the snp
        self.rid = split_line[2] # Name/Record ID
        self.ref_allele = split_line[3]
        self.alt_allele = split_line[4] # The alterate allele according to the vcf; also a string 
        # Only accept snps in ACGT. 
        if self.ref_allele not in ["A","C","G","T"] or self.alt_allele not in ["A","C","G","T"]:
            self.failed = True
            self.failure_cause = "ref or alt not in ACGT"
            return
        filter = split_line[6]  # See vcf format specifications for the interpretation of 
                                    # the filter field
        if filter not in ['PASS', '.'] :  # PASS indicates a SNP that passed all QC filters.
            self.failed = True
            self.failure_cause = filter
            return
              
        genotype_strings = split_line[9:]

        # Prepare a list that will contain the transformed genotypes. 
        # Since we already know how long the list will be, it makes sense 
        # to create an array of zeros of the same length as self.gtypes, 
        
        self.genotype_array = np.zeros(len(genotype_strings), dtype = np.int8)             

        # Count the number of each genotype. 
        # There may be different strings giving the same genotype so we increment the 
        # counts found so far for the genotype by the number of times the  
        # For example, "0/0" and "0\0" give homref, and "0|1" and "1|0" give het
        
        n_missing = 0
        for index,genotype_string in enumerate(genotype_strings):
            if genotype_string == './.':
                n_missing +=1 
                self.genotype_array[index]=-1
                continue # missing data will be left as 0
            allele_0 = genotype_string[0] # Get the first allele (as a string)
            allele_1 = genotype_string[2]
            if (allele_0=='1' and allele_1=='1'): # Use rstrip because windows machines will occasionally have extra \n
                self.genotype_array[index]=2
            elif ((allele_0=='0' and allele_1=='1') or (allele_0=='1' and allele_1=='0')):
                self.genotype_array[index]=1   
            elif (allele_0=='0' and allele_1=='0'):
                # The array was initialized to zero, so nothing to do here!
                continue
            else:
                self.failed=True
                self.failure_cause="unknown genotype " + genotype_string
                return

    def __repr__(self):
        if self.failed:
            return f"SNP(failed={self.failed}, failure_cause={self.failure_cause})"
        return f"SNP(rid={self.rid}, ref={self.ref_allele}, alt={self.alt_allele})"


def read_vcf(path: str, verbose = False) -> list[SNP]:
    snps = []

    # handle both gzipped and uncompressed VCF files
    with (gzip.open if path.endswith(".gz") else open)(path, 'r') as f:
        for line in f:
            snp = SNP(line)
            if snp.failed:
                if verbose:
                    print(f"Reading SNP failed: {snp.failure_cause}")
                continue
            snps.append(snp)

    return snps