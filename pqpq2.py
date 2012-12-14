import sys
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
        ind = sch.fcluster(L,0.2,'distance')#distance is dissmilarity(1-correlation)
        print gene,len(pep_array)
        print D
        print L
        print ind
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
    for i in range(0,10):
        gene=genelist[i]
        main()
        if i%1000==0:
            print i,'genes processed'
    print len(genelist),'genes processed'
    print 'program finished'    
