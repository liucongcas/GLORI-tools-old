import sys
import os
import re
import argparse
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import reverse_complement
import time
import numpy as np
import pysam


parser = argparse.ArgumentParser(description = "building index for reference with three bases")
parser.add_argument("-r", "--reads", nargs="?", type=str, default='read2', help="read1,read2,paired")
parser.add_argument("-f", "--reference", nargs="?", type=str, default=sys.stdin, help="reference file")
parser.add_argument("-p", "--Threads", nargs="?", type=str, default='1', help = "number of alignment threads to launch")
parser.add_argument("-t", "--tools", nargs="?", type=str, default=sys.stdin, help="bowtie,bowtie2,bwa,hisat2,STAR")
parser.add_argument("-mate_length", "--mate_length", nargs="?", type=int, default=100, help="for STAR building index")
parser.add_argument("-pre", "--outname_prefix", nargs="?", type=str, default='default',help = "--outname_prefix")
parser.add_argument("-o", "--outputdir", nargs="?", type=str, default='./', help="outputdir")
args = parser.parse_args()


reads = args.reads
reference = args.reference
Threads = args.Threads
tools = args.tools
mate_length = args.mate_length
outputdir = args.outputdir
outname_prx = args.outname_prefix


def change_reference(reads,reference,outputdir,outname_prx):
    if outname_prx !='default':
        refer_name = outname_prx
    else:
        refer_name = "_".join(os.path.basename(reference).split(".")[:-1])
    if os.path.exists(outputdir):
        pass
    else:
        os.mkdir(outputdir)
    if reads == "read1":
        changed_refer = outputdir + "/" + refer_name + ".TC_conversion.fa"
        subprocess.call("rm -f " + changed_refer + " 2>/dev/null", shell=True)
        changed_refer2 = open(changed_refer, 'a')
        for record in SeqIO.parse(reference, "fasta"):
            record.seq = Seq(re.sub('T', 'C', str(record.seq).upper()))
            record.id = record.id + "_TC_converted"
            SeqIO.write(record, changed_refer2, "fasta")
    elif reads == "read2":
        changed_refer = outputdir+"/" + refer_name + ".AG_conversion.fa"
        subprocess.call("rm -f " + changed_refer + " 2>/dev/null", shell=True)
        changed_refer2 = open(changed_refer,'a')
        for record in SeqIO.parse(reference, "fasta"):
            record.seq = Seq(re.sub('A','G',str(record.seq).upper()))
            # print(record.id)
            record.id = record.id + "_AG_converted"
            # print(record)
            # time.sleep(1)
            SeqIO.write(record, changed_refer2, "fasta")
    return changed_refer


def build_index(changed_refer,tool,Threads):
    if tool == "bowtie":
        subprocess.call('bowtie-build -q '+ changed_refer + ' ' + changed_refer,shell=True)
    elif tool == "bowtie2":
        subprocess.call('bowtie2-build -q '+ changed_refer + ' ' + changed_refer,shell=True)
    elif tool == "bwa":
        subprocess.call('bwa index ' + changed_refer + ' ' + changed_refer,shell=True)
    elif tool == "hisat2":
        subprocess.call('hisat2-build -q ' + changed_refer + ' ' + changed_refer,shell=True)
    elif tool == "STAR":
        filedir_STAR=changed_refer[:-3]
        # subprocess.call('rm -rf ' + filedir_STAR + " 2>/dev/null", shell=True)
        # subprocess.call('mkdir -p ' + filedir_STAR, shell=True)
        fh = pysam.FastaFile(changed_refer)
        genomesize = sum(fh.lengths)
        Nbases = int(round(min(14, np.log2(genomesize)/2 - 1)))
        command2 = 'STAR --runMode genomeGenerate -runThreadN ' + str(Threads) + ' --genomeDir ' + filedir_STAR + \
                 ' --genomeFastaFiles ' + changed_refer + ' --genomeSAindexNbases '+ str(Nbases) + ' --limitGenomeGenerateRAM 84807429045'
        print(command2)
        subprocess.call(command2,shell=True)


if __name__ == "__main__":
    refername2 = os.path.basename(reference)
    print("**********Building index for " + refername2 + " with " + tools+"************")
    refer_name = reference
    changed_refer = outputdir+"/" + refer_name + ".AG_conversion.fa"
    # changed_refer = change_reference(reads,reference,outputdir,outname_prx)
    build_index(changed_refer,tools,Threads)
    print("**********Results will be found in "+outputdir + "************")
