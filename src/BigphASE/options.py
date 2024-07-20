# a interface to options for other tools   
class __Options__:
 
    def __init__(self):
        self._opts = []

    def get_opts(self) -> list:
        return self._opts
    
    def get_opts_str(self) -> str:
        opts =  ' '.join(self._opts)
        return opts

class NucmerOptions(__Options__):
    """
    Class representing options for the NUCmer alignment tool.
    
    Args:
        ref (str): Reference sequence file path.
        query (str): Query sequence file path.
        prefix (str, optional): Prefix or path/Prefix for output files. If not provided, it will be generated based on the reference and query file names.
    """    
    def __init__(self):
        super().__init__()
    
    def set_mum(self):
        """
        Use anchor matches that are unique in both the reference and query (false)
        """
        
        if self._opts.__contains__('--maxmatch'):
            self._opts.remove('--maxmatch')
            self._opts.append('--mum')
        else:
            self._opts.append('--mum')
        return self
        
    
    def set_maxmatch(self):
        """
        Use all anchor matches regardless of their uniqueness (false)
        """
        if self._opts.__contains__('--mum'):
            self._opts.remove('--mum')
            self._opts.append('--maxmatch')
        else:
            self._opts.append('--maxmatch')
        return self
            
    def set_b(self, value:str):
        """
        Set the distance an alignment extension will attempt to extend poor scoring regions before giving up (200)
        """
        self._opts.append('-b')
        self._opts.append(value)
        return self

    def set_c(self, value:str):
        """
        Sets the minimum length of a cluster of matches (65)
        """
        self._opts.append('-c')
        self._opts.append(value)
        return self
    
    def set_d(self, value:float):
        """
        Set the maximum diagonal difference between two adjacent anchors \n
        in a cluster as a differential fraction of the gap length (0.12)
        """
        self._opts.append('-d')
        self._opts.append(value)
        return self
    
    def set_D(self, value:str):
        """
        Set the maximum diagonal difference between two adjacent anchors in a cluster (5)
        """
        self._opts.append('-D')
        self._opts.append(value)
        return self
    
    def set_g(self, value:str):
        """
        Set the maximum gap between two adjacent matches in a cluster (90)
        """
        self._opts.append('-g')
        self._opts.append(value)
        return self
    
    def set_noextend(self):
        """
        Do not perform cluster extension step (false)
        """
        self._opts.append('--noextend')
        return self
    
    def set_f(self):
        """
        Use only the forward strand of the Query sequences (false)
        """
        self._opts.append('-f')
        return self
    
    def set_l(self, value:str):
        """
        Set the minimum length of a single exact match (20)
        """
        self._opts.append('-l')
        self._opts.append(value)
        return self
    
    def set_L(self, value:str):
        """
        Minimum length of an alignment, after clustering and extension (0)
        """
        self._opts.append('-L')
        self._opts.append(value)
        return self
     
    def set_nooptimize(self):
        """
        No alignment score optimization, i.e. if an alignment extension reaches the end of a sequence, \n
        it will not backtrack to optimize the alignment score and instead terminate the alignment at the end of the sequence (false)
        """
        self._opts.append('--nooptimize')
        return self
 
    def set_r(self):
        """
        Use only the reverse complement of the Query sequences (false)
        """
        self._opts.append('-r')
        return self
    
    # def set_p(self, value:str):
    #     """
    #     Write output to PREFIX.delta (out)
    #     """
    #     if self._opts.__contains__('-p'):
    #         self._opts[self._opts.index('-p')+1] = value
    #     else:
    #         self._opts.append('-p')
    #         self._opts.append(value)
    #     return self
        
    def set_delta(self, value:str):
        """
        Output delta file to PATH (instead of PREFIX.delta)
        """
        self._opts.append('--delta')
        self._opts.append(value)
        return self
    
    def set_sam_short(self, value:str):
        """
        Output SAM file to PATH, short format
        """
        self._opts.append('--sam-short')
        self._opts.append(value)
        return self
    
    def set_sam_long(self, value:str):
        """
        Output SAM file to PATH, long format
        """
        self._opts.append('--sam-long')
        self._opts.append(value)
        return self
    
    def set_save(self):
        """
        Save suffix array to files starting with PREFIX
        """
        self._opts.append('--save')
        return self
    
    def set_load(self):
        """
        Load suffix array from file starting with PREFIX
        """
        self._opts.append('--load')
        return self
    
    def set_batch(self):
        """
        Proceed by batch of chunks of BASES from the reference
        """
        self._opts.append('--batch')
        return self
     
    def set_t(self, value:int):
        """
        Use NUM threads (# of cores).
        """
        self._opts.append('-t')
        self._opts.append(str(value))
        return self
    


class Delta2hapOptions(__Options__):
    def __init__(self):
        super().__init__()
        
    def set_C(self):
        """
        out all SNPs (including SNPs found in ambiguous mapped alignments).
        """
        self._opts.append('-C')
        return self
    
    def set_I(self):
        """
        do not out indels
        """
        self._opts.append('-I')
        return self
        
    def set_l(self,gap:int):
        """
        maximum gap between variations in a single haplotype(default 30).
        """
        self._opts.append('-l')
        self._opts.append(str(gap))
        return self
        
    def set_L(self,length:int):
        """
        maximum number of variations in a single haplotype(default 50).
        """
        self._opts.append('-L')
        self._opts.append(str(length))
        return self

    def set_q(self):
        """
        sort by query.
        """
        self._opts.append('-q')
        return self
        
    def set_r(self):
        """
        sort by reference(default).
        """
        self._opts.append('-r')
        return self
        
    def set_filter(self,maximum:int,window_size:int):
        """
        SNPs will be discarded if more than maximum number in window_size bp in a sequence.
        """
        self._opts.append('-f')
        self._opts.append(str(maximum))
        self._opts.append('-F')
        self._opts.append(str(window_size))
        return self
        
        
class Hisat2BuildOptions(__Options__):
    def __init__(self):
        super().__init__()

    def set_a(self):
        """
        disable automatic -p/--bmax/--dcv memory-fitting
        """
        self._opts.append('-a')
        return self
    
    def set_p(self,value:int):
        """
        number of threads
        """
        self._opts.append('-p')
        self._opts.append(str(value))
        return self
    
    def set_r(self):
        """
        don't build .3/.4.ht2 (packed reference) portion
        """
        self._opts.append('-r')
        return self
    
    def set_3(self):
        """
        just build .3/.4.ht2 (packed reference) portion
        """
        self._opts.append('-3')
        return self
    
    def set_o(self,value:int):
        """
        SA is sampled every 2^offRate BWT chars (default: 5)
        """
        self._opts.append('-o')
        self._opts.append(str(value))
        return self
    
    def set_t(self,value:int):
        """
        of chars consumed in initial lookup (default: 10)
        """
        self._opts.append('-t')
        self._opts.append(str(value))
        return self
    
    def set_localoffrate(self,value:int):
        """
        SA (local) is sampled every 2^offRate BWT chars (default: 3)
        """
        self._opts.append('--localoffrate')
        self._opts.append(str(value))
        return self
    
    def set_localftabchars(self,value:int):
        """
        of chars consumed in initial lookup in a local index (default: 6)
        """
        self._opts.append('--localftabchars')
        self._opts.append(str(value))
        return self
    
# class GraphBuild will take care of this two options
    # def set_snp(self,value:str):
    #     """
    #     SNP file name
    #     """
    #     self._opts.append('--snp')
    #     self._opts.append(value)
    #     return self   
    
    # def set_haplotype(self,value:str):
    #     """
    #     haplotype file name
    #     """
    #     self._opts.append('--haplotype')
    #     self._opts.append(value)
    #     return self
    
    def set_ss(self,value:str):
        """
        Splice site file name
        """
        self._opts.append('--ss')
        self._opts.append(value)
        return self
    
    def set_exon(self,value:str):
        """
        Exon file name
        """
        self._opts.append('--exon')
        self._opts.append(value)
        return self
    
    def set_repeat_ref(self,value:str):
        """
        Repeat reference file name
        """
        self._opts.append('--repeat-ref')
        self._opts.append(value)
        return self
    
    def set_repeat_info(self,value:str):
        """
        Repeat information file name
        """
        self._opts.append('--repeat-info')
        self._opts.append(value)
        return self
    
    def set_repeat_snp(self,value:str):
        """
        Repeat snp file name
        """
        self._opts.append('--repeat-snp')
        self._opts.append(value)
        return self
    
    def set_repeat_haplotype(self,value:str):
        """
        Repeat haplotype file name
        """
        self._opts.append('--repeat-haplotype')
        self._opts.append(value)
        return self

class Hisat2Options(__Options__):
    def __init__(self):
        super().__init__()
        
    def set_q(self):
        """
        Reads are FASTQ files. FASTQ files usually have extension .fq or .fastq. \n
        FASTQ is the default format.
        """
        self._opts.append('-q')
        
        return self

    def set_qseq(self):
        """
        Reads are QSEQ files. QSEQ files usually end in _qseq.txt.
        """
        self._opts.append('--qseq')
        
        return self
    
    def set_f(self):
        """
        Reads are FASTA files. FASTA files usually have extension .fa, .fasta, .mfa, .fna or similar.
        """
        self._opts.append('-f')
        
        return self

    def set_r(self):
        """
        Reads are files with one input sequence per line,\n
        without any other information (no read names, no qualities).
        """
        self._opts.append('-r')
        
        return self
    
    def set_skip(self,value:int):
        """
        Skip (i.e. do not align) the first <int> reads or pairs in the input.
        """
        self._opts.append('--skip')
        self._opts.append(str(value))
        
        return self
    
    def set_p(self,value:int):
        """
        Launch NTHREADS parallel search threads (default: 1).
        """
        self._opts.append('-p')
        self._opts.append(str(value))
        
        return self
    
              