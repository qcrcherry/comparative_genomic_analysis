#!/usr/bin/env python3
"""
Created on Fri Feb 17 13:52:31 2023

@author: qcr
"""

from collections import defaultdict
import sys,os
import argparse
from Bio import SeqIO
from multiprocessing import Pool
import subprocess




# arguments for the script
def get_arguments():
    parser = argparse.ArgumentParser(
        description="""check the completeness of specific genes in bacterial contigs""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group('required arguments')

    # input fasta file
    required.add_argument('-q', '--query', action='store',
                          required=True,
                          help='reference sequence')
    required.add_argument('-p', '--prefix', action='store',
                          help='output files prefix')

    required.add_argument('-s', '--subject', action='store',nargs='+',
                        help="subject sequence files",
                        required=True)
    parser.add_argument('-t', '--thread', action='store',default=8,
                        help="numbers of thread used")
    # output verbosity
    parser.add_argument(
        "-v",
        "--verbose",
        const=1,
        default=0,
        type=int,
        nargs="?",
        help="Increase verbosity: 0 = only warnings, 1 = info, 2 = debug. No number means info. Default is no verbosity.")


    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    return args


def merge_fasta(seqs):
    rec=[]
    for i in seqs:
        name=os.path.basename(i).split('.')[0]
        for seq_record in SeqIO.parse(i, "fasta"):
            seq_record.id=str(name)+'|'+seq_record.name
            seq_record.name=''
            rec.append(seq_record)
    with open('merged.fasta','w') as f:
        SeqIO.write(rec,f,'fasta')

def merge_fasta_new(seq):
    records = []
    name=os.path.basename(seq).split('.')[0]
    for seq_record in SeqIO.parse(seq, "fasta"):
        seq_record.id=str(name)+'|'+seq_record.name
        seq_record.name=''
        records.append(seq_record)
    return records

def run_nucmer():
    args = get_arguments()
    command = ["nucmer","-t",args.thread, "-p",args.prefix,args.query,'merged.fasta']  
    subprocess.run(command,check=True)  

def run_show_coords():
    args = get_arguments()
    command = ["show-coords","-c","-l","-H","-T",args.prefix+'.delta'] 
    output = subprocess.run(command,check=True, capture_output=True, text=True) 
    with open (args.prefix+'.coords','w') as f:
        f.writelines(output.stdout)

def main():
    args = get_arguments()
    #merge_fasta(args.subject)
    #多线程合并
    merged_records = []
    with Pool(processes=int(args.thread)) as pool:
        for records in pool.map(merge_fasta_new, args.subject):
            merged_records += records
    with open("merged.fasta", "w") as output_handle:
        SeqIO.write(merged_records, output_handle, "fasta")
    #运行MUMER
    run_nucmer()
    run_show_coords()
    #处理结果
    dic_info=defaultdict(dict)
    with open (args.prefix+'.coords','r') as f:
        for line in f:
            info_list=line.strip().split('\t')
            unique_id=info_list[-1].split('|')[0]+info_list[-2]
            dic_info[unique_id]['query']=info_list[-2]
            dic_info[unique_id]['subject']=info_list[-1].split('|')[0]
            dic_info[unique_id]['segments']=dic_info[unique_id].get('segments',0)+1
            dic_info[unique_id]['max_coverage']=max(dic_info[unique_id].get('max_coverage',0),float(info_list[9]))
            dic_info[unique_id]['covered_positions']=dic_info[unique_id].get('covered_positions',defaultdict(dict))
            dic_info[unique_id]['query_length']=int(info_list[7])
            for i in range(int(info_list[0]), int(info_list[1])+1):
                if i not in dic_info[unique_id]['covered_positions']:
                    dic_info[unique_id]['covered_positions'][i] = True
                else:
                    dic_info[unique_id]['covered_positions'][i] = dic_info[unique_id]['covered_positions'][i] or True
    with open(args.prefix+'_output.tsv','w') as f:
        f.writelines('\t'.join(['query','subject','overall_coverage','segments','max_coverage','query_length']))
        for item in dic_info.values():
            overall_coverage=sum(item['covered_positions'].values()) / item['query_length']*100
            f.writelines('\n'+'\t'.join([item['query'],item['subject'],str(overall_coverage),str(item['segments']),str(item['max_coverage']),str(item['query_length'])]))



if __name__ == '__main__':
    main()
