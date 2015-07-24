import sys
import os
import getopt
import numpy as np

def get_unique_pep(array): #this function is to find unique peptides in the psmarray
    dic={}
    for i in range(0,len(array)):
        pep=array[i][0]
        if pep not in dic:
            dic[pep]=1
    return dic

def main():
    unique_pepdic=get_unique_pep(proteinarray)
    for pep in unique_pepdic.keys():
        psm_quant=[]
        for j in range(0,len(proteinarray)):
            if proteinarray[j][0]==pep:
                psm_quant.append(proteinarray[j][2:])
    
        a=np.array(psm_quant,dtype=float)
        if method=="mean":
            mean=np.mean(a,axis=0)
            mean_round=[ '%.2f' % elem for elem in mean ]
            newline="%s\t%s\t%s\n" % (pep,protein,'\t'.join(mean_round))
            output.write(newline)
        
        elif method=="median":
            median=np.median(a,axis=0)
            median_round=[ '%.2f' % elem for elem in median ]
            newline="%s\t%s\t%s\n" % (pep,protein,'\t'.join(median_round))
            output.write(newline)

if __name__=='__main__':
    ################  Default  ################
    method="median"
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! wrong command, please read the mannual in Readme.txt."
        print "Example: python mergepsm.py --in psmtable.txt -out peptidetable.txt --method median"
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['method=','in=','out='])
        for opt, arg in options:
            if opt == '--method': method=arg
            elif opt == '--in': infilename=arg
            elif opt == '--out': outfilename=arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()

    handle=open(infilename,'r')
    header=handle.readline().split('\t')
    samplesize=len(header[2:])

    output=open(outfilename,'w')
    output.write("\t".join(header))

    protein_psmarray={}
    N=0 #count number of PSMs with missing values
    for line in handle:
        row=line.strip().replace(",",".").replace('\"',"").split("\t")
        try:
            pep=row[0]
            protein=row[1]
        except IndexError:
            print row
            print "check your file format! Program aborted"
            break;
        if protein=="":#this will discard PSMs with unknown identity
            continue;
        if len(row[2:])!=samplesize:
            N+=1
            continue;
        try:
            np.array(row[2:],dtype=float)
            if protein not in protein_psmarray:
                protein_psmarray[protein]=[row]
            else:
                protein_psmarray[protein].append(row)
        except ValueError:
            N+=1
            continue;

    print "skip %d PSMs with missing quant values" % N
    print 'there are',len(protein_psmarray),'protein groups in total'
    handle.close()
    
    i=0
    for protein in protein_psmarray.keys():
        i+=1
        proteinarray=protein_psmarray[protein]
        main()
        if i%1000==0:
            print i,'protein groups processed'
    
    print len(protein_psmarray),'proteins processed'
    print 'program finished'
    print "peptides quantities were calculated as %s of PSMs quantities" % method
