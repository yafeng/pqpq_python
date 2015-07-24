import sys
import os
import getopt
import numpy as np
import scipy.cluster.hierarchy as sch

def formatoutput(array): #format output
    for i in range(0,len(array)):
        newline='\t'.join(array[i]) +'\n'
        output.write(newline)

def convert2ratio(l,n): # l is a list of numbers, n is the reference channel
    l2=[]
    for x in l:
        l2.append(float(x)/float(l[n-1]))
    return l2

def main(): #clustering and write output
    matrix=[]
    for i in range(0,len(proteinarray)):
        if calculate_ratio=="True":
            ratio_array=convert2ratio(proteinarray[i][2:],ref)
            matrix.append(ratio_array)
        else:
            matrix.append(proteinarray[i][2:])

    dataMatrix=np.array(matrix,dtype=float)
    if log_transform=="True":
        dataMatrix=np.log2(dataMatrix)

    if len(proteinarray)>1:
        d = sch.distance.pdist(dataMatrix,metric)# vector of pairwise distances
        if metric=="correlation":
            D = np.clip(d,0,2) #when using correlation, all values in distance matrix should be in range[0,2]
        else:
            D=d
        try:
            cutoff=float(t)
        except ValueError:
            print "please provide a numeric value for --t"; sys.exit()
        L = sch.linkage(D, method,metric)
        ind = sch.fcluster(L,cutoff,'distance')#distance is dissmilarity(1-correlation)
        p=np.array(proteinarray)[:,[0,1]] #slice first and second column of original data
        p=np.concatenate((p,dataMatrix),axis=1) # replace transformed data
        p=np.column_stack([p,ind]) # add cluster result to the last column
        formatoutput(p)
    else:
        p=np.array(proteinarray)[:,[0,1]]
        p=np.concatenate((p,dataMatrix),axis=1)
        p=np.column_stack([p,[0]])
        formatoutput(p)

if __name__=='__main__':
    ################  Default  ################
    method = 'complete'
    metric = 'correlation'
    t=0.4 #default distance cutoff for clustering by pearson correlation
    log_transform='True'
    calculate_ratio='True'
    ref=1
    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! wrong command, please read the mannual in Readme.txt."
        print "Example: python clusterpeptide.py --i heavy_pepdata.txt --o heavy_pepcluster.txt --metric correlation --method complete --t 0.4"
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['metric=',
                                                             'method=',
                                                             'ratio=',
                                                             'reference=',
                                                             'log=',
                                                             't=',
                                                             'in=',
                                                             'out='])
        for opt, arg in options:
            if opt == '--metric': metric=arg
            elif opt == '--method': method=arg
            elif opt == '--t': t=arg
            elif opt == '--in':infilename=arg
            elif opt == '--out':outfilename=arg
            elif opt == '--log': log_transform=arg
            elif opt == '--reference':ref=arg
            elif opt == '--ratio':calculate_ratio=arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()

    print metric,"metric is used"
    print "linking method is",method
    print "log transformation",log_transform
    print "distance threshold of breaking into seperate clusters is",t
    
    handle=open(infilename,'r')
    output=open(outfilename,'w')
    header=handle.readline().strip().split("\t")
    samplesize=len(header[2:])
    
    if ref>samplesize:
        print "ERROR! Reference channel greater than sample size!";sys.exit()
    header.append("cluster")
    output.write("\t".join(header)+"\n")

    protein_peparray={}
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
            if protein not in protein_peparray:
                protein_peparray[protein]=[row]
            else:
                protein_peparray[protein].append(row)
        except ValueError:
            N+=1
            continue;

    print "skip %d PSMs with missing quant values" % N
    print 'there are',len(protein_peparray),'protein groups in total'
    handle.close()
    
    i=0
    for protein in protein_peparray.keys():
        i+=1
        proteinarray=protein_peparray[protein]
        main()
        if i%1000==0:
            print i,'proteins processed'
    print len(protein_peparray),'proteins processed'
    print 'program finished'    
