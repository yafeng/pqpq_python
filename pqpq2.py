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
=======
import numpy
import scipy.cluster.hierarchy as sch

def getgene(infile): #get gene list out of infile 
    dic={}
    for line in infile:
        gene=line.split('\t')[0]
        if gene not in dic:
            dic[gene]=1

    return dic.keys();
def extractpep(gene,infile):#extract all pep for one gene
    array=[]
    for line in infile:
        if line.split('\t')[0]==gene:
            array.append(line[:-1].split('\t')); 
    return array

def main():
    infile=open(sys.argv[1],'r')
    pep_array=extractpep(gene,infile)
    infile.close()
    if len(pep_array)>1:
        matrix=[]
        for i in range(0,len(pep_array)):
            matrix.append(pep_array[i][3].split(','))

        dataMatrix=numpy.array(matrix,dtype=float)
        D = sch.distance.pdist(dataMatrix,'correlation')# vector of pairwise distances
        L = sch.linkage(D, 'single','correlation')
        ind = sch.fcluster(L,0.3,'distance')#distance is dissmilarity(1-correlation)
        p=numpy.array(pep_array)
        p=numpy.column_stack([p,ind])
        formatoutput(p)
    else:
        p=numpy.array(pep_array)
        p=numpy.column_stack([p,[0]])
        formatoutput(p)
        
def formatoutput(array):
    for i in range(0,len(array)):
        newline='\t'.join(array[i])+'\n'
        output.write(newline)
    
    
if __name__=='__main__':
    handle=open(sys.argv[1],'r')# take sample_pepdata.txt as input file
    handle.readline()
    newheader=['gene','pep','PSM count','foldchange','standard_dev','cluster']
    firstline='\t'.join(newheader)+'\n'
    genelist=getgene(handle)
    print 'there are',len(genelist),'genes in total'
    handle.close()
    
    output=open(sys.argv[2],'w')
    output.write(firstline)
    for i in range(0,len(genelist)):
        gene=genelist[i]
        main()
        if i%1000==0:
            print i,'genes processed'
    print len(genelist),'genes processed'
    print 'program finished'
