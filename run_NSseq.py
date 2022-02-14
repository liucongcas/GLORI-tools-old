import os, sys
import argparse
from collections import defaultdict,Counter
import time
import subprocess
import os
import re
import multiprocessing


def run_command(file,combine,untreated,local):
    file3 = outputdir+"/"+prx+"_un_rs.bam"
    file4 = outputdir+"/"+prx+".trans2Genome.bam"
    file5 = outputdir+"/"+prx+"_rs.bam"
    file6 = outputdir+"/"+prx+".trans2Genome.sorted.bam"
    file7_1 = outputdir+"/"+prx+"_merged.bam"
    file7_2 = outputdir+"/"+prx+"_merged.sorted.bam"
    file8 = outputdir+"/"+prx+".pileup"
    chr_file = outputdir + "/" + prx + "_chrlist"
    global file9
    file9 = outputdir+"/"+prx+".referbase.mpi"
    mapping_1 = "python "+NStoolsdir+"mapping_reads.py -q "+ file +" -p "+ Threads + " -f "+ genome
    mapping_2 = " -mulMax " + mulMax + " -t " + tool + " -m " + mismatch +" -pre "+ prx+ " -o "+ outputdir
    print(untreated,combine)
    if untreated:
        file3 = outputdir + "/" + prx + "_un_s.bam"
        file5 = outputdir + "/" + prx + "_s.bam"
        mapping_1 = mapping_1 + " --flag 4 "
        if combine and local:
            mapping_command = mapping_1 + " -Tf " + transgenome + mapping_2 + " --untreated --combine --local "
        elif combine and not local:
            mapping_command = mapping_1 + " -Tf " + transgenome + mapping_2 + " --untreated --combine "
        elif not combine and not local:
            mapping_command = mapping_1 + mapping_2 + " --untreated "
        print(mapping_command)
        if combine:
            print("**************combine,untreated")
            subprocess.call(mapping_command, shell=True)
            print("python " + NStoolsdir + "Transcrip2genome.py --input " + file3 + " --output " + file4 + " --anno " + anno + " --fasta " + genome + " --untreated --sort --index")
            subprocess.call("python " + NStoolsdir + "Transcrip2genome.py --input " + file3 + " --output " + file4 \
                            + " --anno " + anno + " --fasta " + genome + " --untreated --sort --index",shell=True)
            subprocess.call("python " + NStoolsdir + "concat_bam.py -i " + file5 + " " + file6 + " -o " + file7_1 + " -t " + Threads + " --sort --index ",
                shell=True)
            finalbam = file7_2
        else:
            print("**************uncombine,untreated")
            subprocess.call(mapping_command, shell=True)
            finalbam = outputdir + "/" + prx + "_s_sorted.bam"
            print("samtools view -F 4 -@ " + Threads + " -bS -h " + file5 + " | samtools sort -@ " + Threads + " > " + finalbam)
            subprocess.call("samtools view -F 4 -@ " + Threads + " -bS -h " + file5 + " | samtools sort -@ " + Threads + " > " + finalbam,
                shell=True)
            subprocess.call("samtools index " + finalbam, shell=True)
    else:
        if combine and local:
            mapping_command = mapping_1 + " -Tf " + transgenome + mapping_2 + " --combine " + " --local "
        elif combine and not local:
            mapping_command = mapping_1 + " -Tf "+ transgenome + mapping_2 + " --combine "
        elif not combine and not local:
            mapping_command = mapping_1 + mapping_2
        print(mapping_command)
        if combine:
            print("**************combine,treated")
            subprocess.call(mapping_command, shell=True)
            print("python "+NStoolsdir+"Transcrip2genome.py --input " + file3 + " --output "+file4 + " --anno "+anno+" --fasta "+genome+" --sort --index")
            subprocess.call("python "+NStoolsdir+"Transcrip2genome.py --input " + file3 + " --output "+file4 + " --anno "+anno+" --fasta "+genome+" --sort --index",shell=True)
            subprocess.call("python "+NStoolsdir+"concat_bam.py -i " + file5 + " " + file6 + " -o " + file7_1 + " -t " + Threads + " --sort --index ",shell=True)
            finalbam = file7_2
        else:
            print("**************uncombine,treated")
            subprocess.call(mapping_command, shell=True)
            finalbam = outputdir+"/"+prx+"_rs_sorted.bam"
            print("samtools view -F 20 -@ " + Threads + " -bS -h " + file5 + " | samtools sort -@ " + Threads + " > " + finalbam)
            subprocess.call("samtools view -F 20 -@ " + Threads + " -bS -h " + file5 + " | samtools sort -@ " + Threads + " > " + finalbam,shell=True)
            subprocess.call("samtools index " + finalbam, shell=True)
    print("python "+NStoolsdir+"pileup_genome_multiprocessing.py -P " + Threads + " -f " + genome + " -i " + finalbam + " -o " + file8)
    subprocess.call("python "+NStoolsdir+"pileup_genome_multiprocessing.py -P "+Threads+" -f "+genome+" -i "+ finalbam +" -o "+file8,shell=True)
    print("python "+NStoolsdir+"get_referbase.py -input " + file8 + " -referFa "+ genome.split(".AG_conversion.fa")[0] + " -outname_prx "\
                    + outputdir +'/'+prx)
    subprocess.call("python "+NStoolsdir+"get_referbase.py -input " + file8 + " -referFa "+ genome.split(".AG_conversion.fa")[0] + " -outname_prx "\
                    + outputdir +'/'+prx, shell=True)
    subprocess.call("rm -f " + file8, shell=True)
    print("cut -f 1 "+file9 + " | sort -u > " + chr_file)
    subprocess.call("cut -f 1 "+file9 + " | sort -u > " + chr_file, shell=True)
    chr_list = [r1.strip().split("\t")[0] for r1 in open(chr_file).readlines()]
    multiprocessing.freeze_support()
    pool = multiprocessing.Pool(processes=int(Threads))  # """Create a process pool with p processes"""
    for chr in chr_list:
        pool.apply_async(func=get_sites, args=(chr,))
    print('Waiting for all subprocesses done...')
    pool.close()
    pool.join()
    print("Sub-processes done.")

    """optain totalsites"""
    # subprocess.call("rm -f " + file9, shell=True)
    final_sites1 = outputdir + "/" + prx + ".totalm6A.txt"
    final_format = outputdir + "/" + prx + ".totalformat.txt"
    final_CR = outputdir + "/" + prx + ".totalCR.txt"
    subprocess.call("cat " + outputdir + "/" + prx + ".referbase.mpi.formatted.txt" + ".* | sed '/#/d' > " + final_format,shell=True)
    subprocess.call("cat " + outputdir + "/" + prx + ".CR.txt" + ".* > " + final_CR, shell=True)
    print("cat " + outputdir + "/" + prx + ".callsites" + "*.3.txt > " + final_sites1)
    subprocess.call("cat " + outputdir + "/" + prx + ".callsites" + "*.3.txt > " + final_sites1,shell=True)
    # subprocess.call("rm -f " + outputdir + "/" + prx + ".CR.txt" + ".* ", shell=True)
    # subprocess.call("rm -f " + outputdir + "/" + prx + ".referbase.mpi.formatted.txt" + ".* ", shell=True)
    # subprocess.call("rm -f " + outputdir + "/" + prx + ".referbase.mpi" + ".* ", shell=True)
    # subprocess.call("rm -f " + outputdir + "/" + prx + ".callsites" + ".* ", shell=True)
    """FDR filter"""
    final_sites2 = outputdir + "/" + prx + ".totalm6A.FDR"
    print("python " + NStoolsdir + "m6A_caller_FDRfilter.py -i " + final_sites1 + ' -o ' + final_sites2\
                    + ' -adp ' + adjP)
    subprocess.call("python " + NStoolsdir + "m6A_caller_FDRfilter.py -i " + final_sites1 + ' -o ' + final_sites2\
                    + ' -adp ' + adjP, shell=True)


def get_sites(chr):
    chr2 = chr.split("_AG_converted")[0]
    baseanno_chr = baseanno + "." + chr2
    file_mpi = file9 + "." + chr2
    file_format = outputdir + "/" + prx + ".referbase.mpi.formatted.txt" + "." + chr2
    file_CR = outputdir + "/" + prx + ".CR.txt" + "." + chr2
    file_sites = outputdir + "/" + prx + ".callsites" + "." + chr2
    subprocess.call("awk \'$1==\"\'\"" + chr + "\"\'\"\' " + file9 + " > " + file_mpi,shell = True)

    if anno != 'None':
        if not os.path.exists(baseanno_chr):
            print("awk \'$1==\"\'\"" + chr2 + "\"\'\"\' " + baseanno + " > " + baseanno_chr)
            subprocess.call("awk \'$1==\"\'\"" + chr2 + "\"\'\"\' " + baseanno + " > " + baseanno_chr,shell = True)
        print("python "+NStoolsdir+"m6A_pileup_formatter.py --db "+ baseanno_chr + " -i " + file_mpi + " -o "+ file_format +" --CR "+ file_CR)
        subprocess.call("python "+NStoolsdir+"m6A_pileup_formatter.py --db "+ baseanno_chr + " -i " + file_mpi + " -o "+file_format +" --CR "+ file_CR,shell=True)
    else:
        print("python " + NStoolsdir + "m6A_pileup_formatter.py " + " -i " + file_mpi + " -o " + file_format + " --CR " + file_CR)
        subprocess.call("python " + NStoolsdir + "m6A_pileup_formatter.py " + " -i " + file_mpi + " -o " + file_format + " --CR " + file_CR,
            shell=True)

    print("python "+NStoolsdir+"m6A_caller.py -i "+file_format + " -o " + file_sites +" -c " + Cov +" -C "+ Counts + " -r "\
                    + minRatio +" -p "+pvalue+" -s "+multiAratio+" -R "+AGRatio +" --cutoff " + Acutoffs +" --CR "+ \
                    background+" --method " + statmethod)
    subprocess.call("python "+NStoolsdir+"m6A_caller.py -i "+file_format + " -o " + file_sites +" -c " + Cov +" -C "+ Counts + " -r "\
                    + minRatio +" -p "+pvalue+" -s "+multiAratio+" -R "+AGRatio +" --cutoff " + Acutoffs +" --CR "+ \
                    background+" --method " + statmethod,shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="run m6A sites with mRNA")

    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-i", "--NSdir", nargs="?", type=str, default=sys.stdin,help = "NSdir")
    group_required.add_argument("-q", "--fastq", nargs="?", type=str, default=sys.stdin,
                        help="fastqfiles with surfix as _1.fq;_1.fastq;_2.fq;_2.fastq")
    group_required.add_argument("-f", "--reference", nargs="?", type=str, default=sys.stdin, help="indexfile")
    group_required.add_argument("-Tf", "--transreference", nargs="?", type=str, default='None',help="transcriptom reference indexfile")

    group_required.add_argument("-pre", "--outname_prefix", nargs="?", type=str, default='default', help="--outname_prefix")
    group_required.add_argument("-o", "--outputdir", nargs="?", type=str, default=sys.stdin, help="outputdir")


    group_mappingfilter = parser.add_argument_group("mapping conditions")

    group_mappingfilter.add_argument("-a", "--anno", nargs="?", type=str, default='None', help="anno")
    group_mappingfilter.add_argument("-b", "--baseanno", nargs="?", type=str, default='None', help="baseanno")
    group_mappingfilter.add_argument("-t", "--tools", nargs="?", type=str, default='STAR', help="bowtie,bowtie2,bwa,hisat2,STAR")
    group_mappingfilter.add_argument("-T", "--Threads", nargs="?", type=str, default='1',help="number of alignment threads to launch")
    group_mappingfilter.add_argument("-mulMax", "--mulMax", nargs="?", type=str, default='1',help="suppress all alignments if > <int> exist")
    group_mappingfilter.add_argument("-m", "--mismatch", nargs="?", type=str, default='2', help="mapping mismatch")
    group_mappingfilter.add_argument("--combine", "--combine", help="whether mapping with changed reads", action="store_true")
    group_mappingfilter.add_argument("--untreated", "--untreated", help="if the input is untreated", action="store_true")
    group_mappingfilter.add_argument("--local", "--local",help="whether use local algrithm", action="store_true")

    # Filter
    group_site = parser.add_argument_group("m6A filter")
    group_site.add_argument("-c", "--coverage", dest="coverage", default='15', type=str, help="A+G coverage, default=10")
    group_site.add_argument("-C", "--count", dest="count", default='5', type=str,
                            help="A count, below which the site will not count, default=3")
    group_site.add_argument("-r", "--ratio", dest="ratio", default='0.1', type=str, help="m6A level/ratio, default=0.1")
    group_site.add_argument("-p", "--pvalue", dest="pvalue", default='0.005', type=str, help="pvalue, default=0.005")
    group_required.add_argument("-adp", "--adjustpvalue", dest="adjustpvalue", default='0.005', type=str, help="adjustpvalue, default=0.005")
    group_site.add_argument("-s", "--signal", dest="signal", default='0.8', type=str,
                            help="signal ratio, equals coverage(under A-cutoff)/coverage, default=0.9")
    group_site.add_argument("-R", "--var_ratio", dest="var_ratio", default='0.8', type=str,
                            help="the ratio cutoff of AG/Total to filter sequencing/mapping errors, default=0.8")
    group_site.add_argument("-g", "--gene_CR", dest="gene_CR", default='0.2', type=str,
                            help="conversion rate, over which a gene will be discarded, default=0.1")
    group_site.add_argument("-N", "--AG", dest="AG_number", default=0, type=int,
                            help="AG count, below which a gene will be discarded, default=0")

    # Statistics
    group_stat = parser.add_argument_group("Statistic method")
    group_stat.add_argument("--method", dest="method", default="binomial", choices=['binomial', 'fisher', 'poisson'],
                            help="statistical method: binomial, fisher exact test, or poisson, default=binomial")
    group_stat.add_argument("--CR", dest="conversion_rate", default="gene", choices=['gene', 'overall'],
                            help="conversion rate used: gene or overall, default=gene")
    group_stat.add_argument("--NA", dest="non_anno", default="ELSE",
                            choices=['ELSE', 'Median', 'Mean', 'ALL', 'discard'],
                            help="which CR to use if no aene annotation, default=ELSE")
    group_stat.add_argument("--cutoff", dest="A_cutoffs", default="3",
                            help="A-cutoffs, 1-10,15,20 or None, seperated by comma, default=3,None")



    options = parser.parse_args()
    global genome,transgenome, outputdir,tool,Threads,mulMax,mismatch,anno,baseanno,prx
    global Cov,Counts,minRatio,pvalue,adjP,multiAratio,AGRatio,geneCR,AG_number_gene,statmethod,background,NAbackground,Acutoffs
    NStoolsdir = options.NSdir+"/pipelines/"
    genome = options.reference
    transgenome = options.transreference
    outputdir = options.outputdir
    tool = options.tools
    Threads = options.Threads
    mulMax = options.mulMax
    mismatch = options.mismatch
    anno = options.anno
    baseanno = options.baseanno
    prx = options.outname_prefix
    Cov = options.coverage
    Counts = options.count
    minRatio = options.ratio
    pvalue = options.pvalue
    adjP=options.adjustpvalue
    multiAratio = options.signal
    AGRatio = options.var_ratio
    geneCR = options.gene_CR
    AG_number_gene = options.AG_number
    statmethod = options.method
    background = options.conversion_rate
    NAbackground = options.non_anno
    Acutoffs = options.A_cutoffs
    if anno == 'None':
        background = 'overall'

    run_command(options.fastq,options.combine,options.untreated,options.local)


