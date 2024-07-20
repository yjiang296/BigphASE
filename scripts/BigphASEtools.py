import os
import sys
import subprocess
import shutil
import multiprocessing
import re
import argparse
import BigphASE


BigphASE_url = 'https://github.com/yjiang296/BigphASE'
envs = [
    {
        "nucmer":{'version_cmd':'nucmer --version','min_version':[4,0]},
    },
    {
        "hisat2":{'version_cmd':'hisat2 --version','min_version':[2,0]},
    },
    {
        "samtools":{'version_cmd':'samtools --version','min_version':[1,0]},
    },
    {
        "bedtools":{'version_cmd':'bedtools --version','min_version':[2,0]},
    },
    
]

# global vars
global t_num_default

def check_os():
    if os.name != 'posix':
        sys.exit("This script can only be run on Linux systems.")
        
# check if the version is begger than minimum requirement.
def version_ok(cmd,min_ver)->bool:
    reg_pattern = r'\bv?(\d+\.\d+\.\d+\w*)\b'
    stdout = subprocess.check_output(cmd, shell=True, text=True)
    match = re.search(reg_pattern, stdout)
    assert match, f"{match},none"
    big = int(match.group(1).split('.')[0])
    mid = int(match.group(1).split('.')[1])
    min_big = min_ver[0]
    min_mid = min_ver[1]
    if big >= min_big and mid >= min_mid:
        return True
    return False

# Check if valid softwares are installed
def check_software_installed(softwares:list):
    bad_softwares = []
    requirements = ''
    for software in softwares:
        for name,info in software.items():
            if shutil.which(name) is None:
                bad_softwares.append(name)
            else:
                if not version_ok(info.get('version_cmd'), info.get('min_version')):
                    bad_softwares.append(name)
            requirements += name +f'(>= {str(info.get('min_version')[0])}.{str(info.get('min_version')[1])})\n'    
    if len(bad_softwares) != 0:        
        sys.exit(f"This script depends on:\n{requirements}\nPlease check whether you have installed or added the correct version to your PATH of these:\n{'\n'.join(bad_softwares)}")


def get_available_threads():
    return multiprocessing.cpu_count()
    
def init():
    global t_num_default
    check_os()
    check_software_installed(envs)
    t_num_suggested = int(get_available_threads() * 0.7)
    t_num_default = t_num_suggested if t_num_suggested > 1 else 1
    
    # del PREFIX.inexon.snp.2l
    # del PREFIX.sorted.inexon.pileup
    # del PREFIX.sorted.inexon.wao
    # del PREFIX.sorted.inexon.bedlike.pileup
    # del PREFIX.genes.bed
    # del PREFIX.bam
    # del PREFIX.sorted.bam
    # del PREFIX.sorted.bam.bai
    # del PREFIX.n.ht2
    # del PREFIX.delta
def clear_output(prefix):
    targets = [f'.{n}.ht2' for n in range(1, 9)] + ['.delta']
    targets += ['.bam', '.sorted.bam', ".sorted.bam.bai"]
    targets += ['.inexon.snp.2l', '.sorted.inexon.pileup', '.sorted.inexon.wao', '.sorted.inexon.bedlike.pileup', '.genes.bed']
    
    for target in targets:
        try:
            os.remove(prefix + target)
        except:
            continue
    

if __name__ == "__main__":
    init()
    
    parser = argparse.ArgumentParser()
    # Required arguments
    parser.add_argument("ref", help="parent_A genome sequence file",
                    type=str)
    parser.add_argument("qry", help="parent_B genome sequence file",
                    type=str)
    parser.add_argument("genome_info", help="parent_A genome gff3 file",
                    type=str)
    parser.add_argument("r1", help="clean RNAseq R1 file(.fq, .fastq, .fa, .fasta, .mfa, .fna, _qseq.txt and their .gz file).",
                    type=str)
    parser.add_argument("r2", help="clean RNAseq R2 file(.fq, .fastq, .fa, .fasta, .mfa, .fna, _qseq.txt and their .gz file).",
                    type=str)
    parser.add_argument("prefix", help="out put file prefix.",
                    type=str)
    
    # Optional arguments
    parser.add_argument("-t", default=0,help="number of threads allow to use(default = 0.7 of available threads or 1).",
                    type=int, metavar="int")
    # For nucmer
    parser.add_argument("-nb", default=200, help="Set the distance an alignment extension will attempt to extend poor scoring regions before giving up (200)."
                        , type=int, metavar="int")
    parser.add_argument("-nc", default=65, help="Sets the minimum length of a cluster of matches (65)."
                        , type=int, metavar="int")
    parser.add_argument("-nl", default=20, help="Set the minimum length of a single exact match (20)."
                        , type=int, metavar="int")
    # For delta2hap
    parser.add_argument("-C", action='store_true', help="out all SNPs (including SNPs found in ambiguous mapped alignments).")
    parser.add_argument("-I", action='store_true', help="do not out indels as snp.")
    parser.add_argument("-dl", default=30, help="maximum gap between variations in a single haplotype(default 30)."
                        , type=int, metavar="int")
    parser.add_argument("-dL", default=50, help="maximum number of variations in a single haplotype(default 50)."
                        , type=int, metavar="int")
    parser.add_argument("-df",nargs=2, help="SNPs will be discarded if more than maximum number in window_size bp in a sequence. Use like -df 3 10"
                        , type=int, metavar="int")
    
    # For hisat2
    parser.add_argument("-r", action='store_true', help="Reads are files with one input sequence per line,without any other information (no read names, no qualities).")  
    args = parser.parse_args()
   
   # Global vars
    REF = args.ref
    QRY = args.qry
    GFF = args.genome_info
    R1:str = args.r1
    R2:str = args.r2
    PREFIX = args.prefix
    
    # thread num
    t_num = args.t if args.t > 0 else t_num_default
    
    
    gb = BigphASE.GraphBuilder(PREFIX, ref=REF, qry=QRY)
    gb.OPT_nucmer.set_t(t_num)\
                .set_b(args.nb)\
                .set_c(args.nc)\
                .set_l(args.nl)

    gb.OPT_delta2hap.set_l(args.dl)\
                    .set_L(args.dL)
    if args.df:
        maxnum, window_size = args.df
        gb.OPT_delta2hap.set_filter(maxnum,window_size)
    if args.C:
        gb.OPT_delta2hap.set_C()
    if args.I:
        gb.OPT_delta2hap.set_I()
        
    gb.OPT_hisatBuild.set_p(t_num)

    
    hisat2_index = PREFIX
    mp = BigphASE.Mapper(PREFIX, hisat2_index, R1, R2)
    mp.OPT_hisat2.set_p(t_num)
    
    if R1.endswith(('.fa', '.fasta', '.mfa', '.fna','.fa.gz', '.fasta.gz', '.mfa.gz', '.fna.gz')):
        mp.OPT_hisat2.set_f()
    elif R1.endswith(('.fq', '.fastq','.fq.gz', '.fastq.gz')):
        mp.OPT_hisat2.set_q()
    elif R1.endswith('_qseq.txt') and R2.endswith('_qseq.txt'):
        mp.OPT_hisat2.set_qseq()
    
    else:
        raise ValueError("Read file format is not supported")
        

    sorted_bam = PREFIX + '.sorted.bam'
    snp = PREFIX + '.snp'
    an = BigphASE.Analyzer(PREFIX, REF, sorted_bam, snp, GFF)

    # Run all model
    gb.run()
    mp.run()
    an.run()
    
    # remove non importent files
    clear_output(PREFIX)
    
        
# Usage:
# python BigphASEtools.py [-options] <ref> <qry> <genome_info> <r1> <r2> <prefix>