import sys
import os
import getopt
import numpy
import scipy.cluster.hierarchy as sch

def formatoutput(array): #format output
    for i in range(0,len(array)):
        newline='\t'.join(array[i]) +'\n'
        output.write(newline)

def main(): #clustering and write output
    if len(psm_array)>1:
        matrix=[]
        for i in range(0,len(psm_array)):
            matrix.append(psm_array[i][9:])
        
        dataMatrix=numpy.log(numpy.array(matrix,dtype=float))
        d = sch.distance.pdist(dataMatrix,metric)# vector of pairwise distances
        if metric=="correlation":
            D = numpy.clip(d,0,2) #when using correlation, all values in distance matrix should be in range[0,2]
        else:
            D=d
        L = sch.linkage(D, method,metric)
        cutoff=0.5*max(L[:,2])
        ind = sch.fcluster(L,cutoff,'distance')#distance is dissmilarity(1-correlation)
        p=numpy.array(psm_array)
        p=numpy.column_stack([p,ind])
        formatoutput(p)
    else:
        p=numpy.array(psm_array)
        p=numpy.column_stack([p,[0]])
        formatoutput(p)

if __name__=='__main__':
    ################  Default  ################
    method = 'average'
    metric = 'euclidean'
    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! wrong command, please read the mannual in Readme.txt."
        print "Example: python pqpq2.py --i psmdata.txt --o pqpqout .txt"
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['metric=',
                                                             'method=',
                                                             'i=','o='])
        for opt, arg in options:
            if opt == '--metric': metric=arg
            elif opt == '--method': method=arg
            elif opt == '--i':infilename=arg
            elif opt == '--o':outfilename=arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()

    print metric,"metric is used"
    
    handle=open(infilename,'r')
    header=handle.readline().split('\t')
    protein_psmarray={}
    for line in handle:
        row=line.strip().split("\t")
        protein=row[2]
        if protein not in protein_psmarray:
            protein_psmarray[protein]=[row]
        else:
            protein_psmarray[protein].append(row)
    
    print 'there are',len(protein_psmarray),'unique proteins in total'
    handle.close()

    header.append("cluster")
    firstline='\t'.join(header)+'\n'
    output=open(outfilename,'w')
    output.write(firstline)
    
    i=0
    for protein in protein_psmarray.keys():
        i+=1
        psm_array=protein_psmarray[protein]
        main()
        if i%1000==0:
            print i,'proteins processed'
    print len(protein_psmarray),'proteins processed'
    print 'program finished'
