import subprocess
import os


from typing import overload

from .options import NucmerOptions,Delta2hapOptions,Hisat2BuildOptions,Hisat2Options
from .utils import ExonTree, target_exist, make_dir, Record

__all__ = ["GraphBuilder","Mapper","Analyzer"]

# ----------------------- PRIVATE ----------------------- #
delta2hap_dir = os.path.join(os.path.dirname(__file__),'delta2hap/bin/delta2hap')

# ----------------------- PUBLIC ----------------------- #

# GraphBuilder 类负责hisat2索引的建立并生成snp和haplotype文件
class GraphBuilder:
    _out_dir = ''
    _out_prefix = ''
    _ref = ''
    _qry = ''
    _delta = ''
    _snp = ''
    _haplotype = ''
    
    @overload
    def __init__(self, prefix:str, ref:str, qry:str): ...
    
    @overload
    def __init__(self, prefix:str, delta:str): ...
    
    @overload
    def __init__(self, prefix:str,snp:str, haplotype:str): ...
    
    
    def __init__(self,prefix: str, ref: str = None, qry: str = None, 
                 delta: str = None, 
                 snp:str = None, haplotype:str = None, 
                 run_single:bool=False) :
        """
        ## parames:
        `prefix`: prefix of out files(required).
        
        - combination1: \n
        `ref`: reference file.
        `qry`: query file.
        
        - combination2: \n
        `delta`: delta file from nucmer.
        
        - ccombination3: \n
        `snp`:snp file.
        `haplotype`:haplotype file.
        """
        if (ref and qry) or delta or (ref and snp and haplotype):
            self.OPT_nucmer = NucmerOptions()
            self.OPT_delta2hap = Delta2hapOptions()
            self.OPT_hisatBuild = Hisat2BuildOptions()
            
            self._out_dir = os.path.abspath(os.path.dirname(prefix))
            make_dir(self._out_dir)
            self._out_prefix = os.path.basename(prefix)
            
            self._run_single = run_single
            if snp:
                if target_exist(snp,haplotype,ref):
                    self._snp = os.path.abspath(snp)
                    self._haplotype = os.path.abspath(haplotype)
                    self._ref = os.path.abspath(ref)
            elif delta:
                with open(delta) as delta_f:
                    line = delta_f.readline().strip()
                    ref_path, qry_path = line.split(' ')[0], line.split(' ')[1]
                    if target_exist(ref_path,qry_path):
                        self._ref = os.path.abspath(ref_path)
                        self._qry = os.path.abspath(qry_path)
                        self._delta = os.path.abspath(delta)  
            else:
                if target_exist(ref,qry):
                    self._ref = os.path.abspath(ref)
                    self._qry = os.path.abspath(qry)
        else:
            raise ValueError("Must provide either (ref, qry, prefix) or (delta, prefix) or (ref,snp,haplotype) combination.")

    # nucmer [options] <reference> <query>
    @Record
    def _run_nucmer(self):
        out_f_prefix = os.path.join(self._out_dir, self._out_prefix)
        self.OPT_nucmer._opts.extend(['-p', self._out_prefix]) #set out put file path
        
        cmd_list = ["nucmer"] + self.OPT_nucmer.get_opts()
        cmd_list.append(self._ref)
        cmd_list.append(self._qry)
        
        result = subprocess.run(cmd_list, stdout = subprocess.PIPE, stderr = subprocess.PIPE, check = True, cwd=self._out_dir, encoding="utf-8")
        if result.returncode == 0:
            self._delta = out_f_prefix + '.delta'
            
        return result
    
    # delta2hap [options] <delta>
    @Record
    def _run_delta2hap(self):
        out_f_prefix = os.path.join(self._out_dir, self._out_prefix)
        self.OPT_delta2hap._opts.extend(['-p',out_f_prefix])
        
        cmd_list = [delta2hap_dir] + self.OPT_delta2hap.get_opts()
        cmd_list.append(self._delta)
        
        result = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, encoding="utf-8")
        if result.returncode == 0:
            self._snp = out_f_prefix + '.snp'
            self._haplotype = out_f_prefix + '.haplotype'
            
        return result
    
    # hisat2-build -p <int> --snp <snp> --haplotype <haplotype> <reference> <prefix>
    @Record
    def _run_hisat_build(self):
        self.OPT_hisatBuild._opts.extend(['--snp',self._snp,'--haplotype',self._haplotype])

        cmd_list = ["hisat2-build"] + self.OPT_hisatBuild.get_opts() + [self._ref, self._out_prefix]
        
        result = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, encoding="utf-8",cwd=self._out_dir)
        
        return result
    
    def run(self):
        if self._run_single:  # self._snp and self._haplotype and self._ref
            if self._snp:
                self._run_hisat_build()
                return
            if self._delta:
                self._run_delta2hap()
                return
            
            self._run_nucmer()
        else:
            if self._snp: # self._snp and self._haplotype and self._ref
                self._run_hisat_build()
            if self._delta:
                self._run_delta2hap()
                self._run_hisat_build()
            else:
                self._run_nucmer()
                self._run_delta2hap()
                self._run_hisat_build()
            
            
class Mapper:
    _out_dir = ''
    _out_prefix = ''
    _ht2_idx_dir = ''
    _ht2_idx_prefix = ''
    _ref = ''
    _read1 =''
    _read2 = ''
    
    def __init__(self,prefix:str,ht2_index:str,read1:str,read2:str):
        
        self._out_dir = os.path.abspath(os.path.dirname(prefix))
        make_dir(self._out_dir)
        self._out_prefix = os.path.basename(prefix)
        
        self._ht2_idx_dir = os.path.abspath(os.path.dirname(ht2_index)) if os.path.dirname(ht2_index) != '' else '.'
        
        if target_exist(self._ht2_idx_dir,read1,read2):
            self._ht2_idx_prefix = os.path.basename(ht2_index)
            self._read1 = os.path.abspath(read1)
            self._read2 = os.path.abspath(read2)
        
        self.OPT_hisat2 = Hisat2Options()
    
    def _clean_seq(self):
        pass
    

    @Record
    def _run_hisat2(self):
        cmd_hisat = "hisat2 " + self.OPT_hisat2.get_opts_str() + " -t -x " + os.path.join(self._ht2_idx_dir, self._ht2_idx_prefix) + " -1 " + self._read1 + " -2 " + self._read2
        cmd_samtools = "samtools " + "view " + "-Sbq " + "30 " + "-o " + os.path.join(self._out_dir, self._out_prefix + ".bam")
        cmd = cmd_hisat + " | " + cmd_samtools
        
        
        result = subprocess.run(cmd,shell=True,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        self._bam = os.path.join(self._out_dir, self._out_prefix + ".bam")
        if result.returncode != 0:
            raise ValueError('reads file type (default FASTQ) may not correct')
        return result
        
    @Record
    def _run_sort_and_index(self):
        cmd_sam_sort = ["samtools","sort",self._bam,"-o",os.path.join(self._out_dir, self._out_prefix + ".sorted.bam")]
        cmd_sam_index = ["samtools","index",os.path.join(self._out_dir, self._out_prefix + ".sorted.bam")]
        
        result = subprocess.run(cmd_sam_sort,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
        subprocess.run(cmd_sam_index,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
        
        return result

    
    def run(self):
        self._run_hisat2()
        self._run_sort_and_index()

class Analyzer:
    _out_dir = ''
    _out_prefix = ''
    _ref = ''
    _sorted_bam = ''
    _snp = ''
    _gff = ''
    
    def __init__(self,prefix:str,ref:str,sorted_bam:str,snp:str,gff:str):
        self._out_dir = os.path.abspath(os.path.dirname(prefix))
        make_dir(self._out_dir)
        if target_exist(ref,sorted_bam,snp,gff):
            self._ref =  os.path.abspath(ref)
            self._sorted_bam = os.path.abspath(sorted_bam)
            self._snp = os.path.abspath(snp)
            self._gff = os.path.abspath(gff)
            self._out_prefix = os.path.basename(prefix)
        

    def _extract_snp_and_gene_list (self):
        snp_in_exon_file = os.path.join(self._out_dir, self._out_prefix + ".inexon.snp.2l")
        gene_list_file =  os.path.join(self._out_dir, self._out_prefix + ".genes.bed")
        
        gff_info = {}
        gene_info =[]
        snps_in_exon = []
        
        with open(self._gff, 'r') as gff:
            for line in gff:
                fields = line.strip().split('\t')
                if len(fields) == 9 and fields[2] == 'exon':
                    if fields[0] not in gff_info:
                        gff_info[fields[0]] = ExonTree()
                    gff_info[fields[0]].insert(int(fields[3]),int(fields[4]))
                    
                if len(fields) == 9 and fields[2] == 'gene':
                    chrom, start, end, gene_id  = fields[0].split("chr")[-1], fields[3], fields[4], fields[8].split(';')[0].split("=")[1]
                    bed_line ='\t'.join([chrom,start,end,gene_id]) + '\n'
                    gene_info.append(bed_line)
                    
        with open(self._snp, 'r') as snp:
            for line in snp:
                fields = line.strip().split('\t')
                chr_id = fields[2]
                snp_pos = int(fields[3])
                if chr_id in gff_info: 
                    if gff_info[chr_id].search(snp_pos):
                        snps_in_exon.append(fields[2] + '\t' + str(snp_pos) + '\n')
            if len(snps_in_exon) == 0:
                raise ValueError('chr name in gff and snp file dos not match')
            
        with open(snp_in_exon_file, 'w') as f1:
            for snp in snps_in_exon:
                f1.write(snp)
        self._snp_in_exon = snp_in_exon_file
        
        with open(gene_list_file,'w') as f2:
            for gene in gene_info:
                f2.write(gene)
        self._gene_list = gene_list_file


    @Record
    def _run_mpileup(self):
        out_pileup = os.path.join(self._out_dir,self._out_prefix + ".sorted.inexon.pileup")
        cmd = ["samtools","mpileup","--no-output-ins","--no-output-del", "--no-output-del", "--no-output-ends","--output-QNAME"]
        cmd += ["-l",self._snp_in_exon] + ["-f",self._ref] + [self._sorted_bam]
        cmd += ["-o",out_pileup]
        result = subprocess.run(cmd,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        self._pileup = out_pileup
        return result
   
    
    def _pileup2csv(self):
        #pileup to bedlike format
        top_bedlike_pileup = os.path.join(self._out_dir,self._out_prefix + ".sorted.inexon.bedlike.pileup")
        cmd = "awk '{print $1,$2,$2,$3,$4,$5,$6,$7}' " + \
            self._pileup + " | " + "tr ' ' '\t' " + "> " + \
            top_bedlike_pileup
        subprocess.run(cmd,shell=True,check=True)

        #intersect bedlike pileup whit gene list
        top_wao = os.path.join(self._out_dir,self._out_prefix + ".sorted.inexon.wao")
        cmd = "bedtools intersect -a " + \
            self._gene_list + " -b " + \
            top_bedlike_pileup + " -wao" + " | " + r"grep -v -w '\-1' " + "> " + \
            top_wao
        subprocess.run(cmd,shell=True,check=True)
        
        # find reads match or mismatch to reference, and remove reads are on both sides
        genes = []
        with open(top_wao,'r') as wao:
            chr, start_pos, end_pos, gene_id = '','','',''
            match_ref, mismatch_ref = set(), set()
            for line in wao:
                fields = line.strip().split('\t')
                if gene_id != fields[3]:
                    tmp_set = match_ref & mismatch_ref
                    match_ref -= tmp_set
                    mismatch_ref -= tmp_set
                    line = '\t'.join([chr, start_pos, end_pos, gene_id, ','.join(match_ref), ','.join(mismatch_ref), str(len(match_ref)), str(len(mismatch_ref))]) + '\n'
                    genes.append(line)
                    match_ref.clear()#dont forget to clear the set
                    mismatch_ref.clear()
                    chr, start_pos, end_pos, gene_id = fields[0],fields[1],fields[2], fields[3]
                
                simbols, reads = fields[9], fields[11].split(',')
                index = 0
                for sim in simbols:
                    if sim in '.,':
                        match_ref.append(reads[index])
                        index +=1
                    elif sim in 'ACGTacgt':
                        mismatch_ref.append(reads[index])
                        index +=1
                    else:
                        continue
            # add lastest gene info 
            tmp_set = match_ref & mismatch_ref
            match_ref -= tmp_set
            mismatch_ref -= tmp_set
            line = '\t'.join([chr, start_pos, end_pos, gene_id, ','.join(match_ref), ','.join(mismatch_ref), str(len(match_ref)), str(len(mismatch_ref))]) + '\n'
            genes.append(line)
            del match_ref, mismatch_ref 
            
            op_csv = os.os.path.join(self._out_dir,self._out_prefix + "SNPs.inexon.csv")
            with open(op_csv,'w') as csv:
                for line in genes:
                    csv.write(line)
        
       
    def run(self):
        self._extract_snp_and_gene_list()
        self._run_mpileup()
        self._pileup2csv()