#!/usr/bin/env python
#coding=utf-8
'''
Created on Dec 22, 2019
author: Dongfang Hu
Email:df2019@gmail.com
'''
import os
import sys
import glob
import argparse
import ConfigParser
parser=argparse.ArgumentParser(prog='10Xoutdirsummary',usage='%(prog)s [opthions] [value]',description='This program is used to produce outdir summary!')
parser.add_argument('-IN','--input',help='the input file',metavar='')
parser.add_argument('-OU','--outdir',help='the output dir',metavar='')
argv=vars(parser.parse_args())

#argparse
if argv['input'] == None:
    raise Exception('You should provide input file')
else:
    input=argv['input'].strip()
if argv['outdir'] == None:
    outdir=os.getcwd()
else:
    outdir=argv['outdir'].strip()
    sample = outdir.split('/')[-1]
input = open(input,"r")
output = open(outdir+"/Cell.summary.txt","w")

#Sequencing,Mapping,Cells
header = input.readline().strip().split(",")
info = input.readline().strip().replace('",',':').replace('%,','%:').replace('"','').split(":")
output.write("[Sequencing]"+"\t"+sample+"\n"+header[3]+"\t"+info[3]+'\n'+"Base of Reads(G)"+'\t'+str(round(int(info[3].replace(',',''))*300/1000000000,2))+'\n'+header[4]+"\t"+info[4]+'\n'+header[5]+"\t"+info[5]+'\n'+header[6]+"\t"+info[6]+'\n'+header[7]+"\t"+info[7]+'\n'+header[8]+"\t"+info[8]+'\n'+"[Mapping]"+"\t"+sample+"\n"+header[9]+"\t"+info[9]+'\n'+header[10]+"\t"+info[10]+'\n'+header[11]+"\t"+info[11]+'\n'+header[12]+"\t"+info[12]+'\n'+header[13]+"\t"+info[13]+'\n'+header[14]+"\t"+info[14]+'\n'+header[15]+"\t"+info[15]+"\n"+"[Cells]"+"\t"+sample+"\n"+header[0]+"\t"+info[0]+'\n'+header[16]+"\t"+info[16]+'\n'+header[1]+"\t"+info[1]+'\n'+header[17]+"\t"+info[17]+'\n'+header[2]+"\t"+info[2]+'\n'+header[18]+"\t"+info[18]+'\n')
