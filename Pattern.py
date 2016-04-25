#********Library class to deal with making and manipulating CnaB library/cladogram/phyolgenic tree*************

#1. main()  			for general testing of Pattern.py functions

#2. isRepeat(test, List): 	looks through list of records for repeats

#3. listseqRec(infile):		returns list of seq records

#4. listseq(infile)		returns list of sequences from .fasta file

#5. discern_ebox(infile)	creates full list of eboxes for a fasta file then prompts user to discern real ones

#6  writeFromFasta()            ouputs hits for given fasta file with particular search parameters 

#7. norm(motifcount, seq)	normalizes motifcount

#8. findEnclosed(string, start, end)       given a start and end character, ret    urns what is in between

#9.  CnaBsearch(seq):		takes in sequence returns list of CnaBs 

#10. LPXTGsearch(seq):

#11. CnaBsearch(seq):     

from Bio.Seq import Seq
from Bio import SeqIO
import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast import Record
import numpy as np
import matplotlib.pyplot as plt

#1. for general testing
def main():
    string='abcde'
    print reverse(string)
#2. looks through list of records for repeats

def isRepeat(test, List):
    print test
    if test!=None:
	    for el in List:
		if el==test:
		    return 0 #this thing is a repeat
	    List.append(test)
	    return 1 #this is new protein
	    
#3. returns list of seq records
def listseqRec(infile):
    l=[]
    l=SeqIO.parse(infile, "fasta")
    return l

def listDef(infile):
    l=[]
    for seq_record in SeqIO.parse(infile, "fasta"):
        l.append(seq_record.description)
    return l


#4. returns list of sequences from .fasta file
def listseq(infile):
    l=[]
    for seq_record in SeqIO.parse(infile, "fasta"):
        l.append(seq_record.seq)
    return l

#5. creates full list of eboxes for a fasta file then prompts user to discern real ones
def discern_ebox(infile):
    eboxlist=[]#stores absolute location of catalytic site for each CnaB
    output=open("CnaBorder.txt", "w")#file containing inforamtion in eboxlist, order is same as infile
    for seq_record in SeqIO.parse(infile, "fasta"):
        print(seq_record)
        single=CnaBsearch(seq_record.seq)
        
        if len(single)==1 or len(single)==0:
        #if one or no eboxes, writes directly to list/file
            eboxlist.append(single[0])
            output.write(str(single[0]))
            output.write("\n")
            print "eboxlist is"+str(eboxlist)
        else:
        #promts user to discern correct ebox
            print("This is :"+ seq_record.name)
            print("Possible eboxes at positions: ")
            for ebox in  single:
                print ("Position: ", ebox)
                print "Ebox sequnce: "
                print seq_record.seq[(ebox-4): (ebox+1)]
            input1=raw_input("Among these, enter the correct location")
            correctspot=int(input1)
            eboxlist.append(correctspot)
            output.write(str(correctspot))
            output.write("\n")
    output.close()
    return eboxlist 
                
#6 ouputs hits for given fasta file with particular search parameters        
def writeFromFasta():
    print "REMEMBER TO MAKE SURE YOU ARE SEARCHING AGAINST RIGHT FORMAT"
    print "enter inpute file name"
    infile=raw_input()
    print "enter output file name"
    filename=raw_input()
    output=open(filename, 'w')
#    l=[]
    for seq_record in SeqIO.parse(infile, "fasta"):
	my_seq=seq_record.seq
        #***** change to correct search thing
        motifList=DGBox(my_seq)
        if len(motifList)>0:

        #### writing the file, can change output parameters at will
	    output.write('\n')
	    output.write(seq_record.description)
            for hit in motifList:
	        output.write('\n')
                motifmsg1="Motif Location: "+str(hit)
                output.write(motifmsg1)
                output.write('\n')
	        #motifmsg2="Surrounding Sequence: "+ str(seq_record.seq[(motifList[hit]):(motifList[hit]+10)])
	        #output.write(motifmsg2)
            output.write("\n \n \n")

    output.close()




#7. normalizes motifcount
def norm(motifcount, seq):
    return motifcount/(1.0*len(seq))

#9. given a start and end character, returns what is in between. LIMITATIONS: start and end can only appear once in the string
#used to find name of organism from alignment data, name is enclosed by '[' and ']'
def findEnclosed(string, start, end):
    a=string.find(start)
    #print "start is ", str(a)
    b=string.find(end)
    #print "end is ", str(b)
    enclosed=string[a+1: b]
    #print enclosed
    if (a!=-1 and b!=-1):
        return enclosed
    else: 
	return None 
####THIS IS SLOW: OPTIMIZE LATER  
def reverse(string):
    backIndex=len(string)-1
    revString=[]
    while (backIndex>=0):
        revString.append(string[backIndex])
        backIndex-=1

    return ''.join(revString)

    
#8. takes in sequence returns list of CnaBs  
def CnaBsearch(seq):
    l=[]
    length=len(seq)
    fy=0
    ilvf=2
    e=4
    motifCount=0
    while e<length:
        if seq[e]=='E':
            if seq[ilvf]=='I' or seq[ilvf]=='L' or seq[ilvf]=='V' or seq[ilvf]=='F':
                if seq[fy]=='F' or seq[fy]=='Y':
                    l.append(e)
                    #print(str('Potential pivot at '), e)
                    #print (seq[fy:e+1])
                    motifCount+=1
        e+=1
        fy+=1
        ilvf+=1
    #print "list is :"
    #print l
    return l

#9 takes in sequence and return list of LPXTG
def LPXTGsearch(seq):
    l=[]
    length=len(seq)
    print "length is ", length
    L=0
    P=1
    T=3
    G=4
    while G<length:
	if seq[L]=='L':
	    if seq[P]=='P':
		if seq[T]=='T':
		    if seq[G]=='G':
			l.append(L)
			print seq[L:G+1]
 			print L
			print '\n'
        G+=1
	T+=1
	P+=1
	L+=1
	
    return l

#Functions for CnaA Definistion
#may implement ranking system depending on sucess rate
def CnaAsearch(seq):
    List=[]
    length=len(seq)
    print "length is ", length
    k=0
    g=2
    w=15
    l=17
    n=20
    d=35
    g2=40
    y=112
    n2=128
    while n2<length:
	if seq[k]=='K':
	    if seq[g]=='G':
 		if seq[w]=='W':
		    if seq[l]=='L':
			if seq[n]=='N':
			    if seq[d]=='D':
				if seq[g2]=='G':
				    if seq[y]=='Y':
					if seq[n2]=='N':
					    List.append[k]
					    print "list is ", List
					    print seq[k : l]
					    print '\n'
	k+=1
	g+=1
	w+=1
	l+=1
	n+=1
	d+=1
	g2+=1
	y+=1
	n2+=1
    print "l is", List					     
    return List

   
def DGBox(seq):
    list1=[]
    length=len(seq)
    d=0
    g=5
    while g<length:
        if seq[d]=='D':
            if seq[g]=='G':
                    list1.append(d)
                    
        d+=1
        g+=1
    print "list is" ,list1
    return list1
##This assumes that sequence is in phase (1st codon starts at base 0)
def findStop(seq):
    stop1=seq.find("TGA")
    if stop1!=0 and stop1%3==0:
        return seq[0:stop1]
    stop2=seq.find("TAA")
    if stop2!=0 and stop2%3==0:
        return seq[0:stop2]
    stop3=seq.find("TAG")
    if stop3!=0 and stop3%3==0:
        return seq[0:stop3]
    return seq


   ########### Reads in any motif (or stdin)###############


class Motif:
    def __init__(self, motifInfo):
        self.motif=self.parse(motifInfo)## list of AAs
        self.motifInfo=motifInfo ## input string, with each proto-AA separated by spaces
        self.length=len(self.motif)
        
    def fill(self, motif,  value): ## value is agreggate data for string
        #print "value is ", value
        #print "motif so far is ", self.printf()
        posi=[] ##list to be turned into AA
        i=0
        for letter in value:
            posi.append(letter)
         
        aa=AA(posi)
        motif.append(aa)
    ####overload compare equality:
    def __eq__ (self, other):
        i=0 #counter for amino acid number
        while i<self.length:
            if self[i]!=other[i]:
                return False
            i+=1
        return True
    #####override !=
    def __ne__(self, other):
        return not self.__eq__(other)
    ####override get item
    def __getitem__(self, index):
        return self.motif[index]

    ############   Reads in info from file or stdin. whatever
    def parse(self, motifInfo):
        aaList=motifInfo.split()
        motif=[]##to be returned 
        for spot in aaList:
            self.fill(motif, spot)
        return motif
        

    def printf(self):
        for el in self.motif:
           print str(el)
    
    def __str__(self):
        return self.motifInfo

    def search(self, seq): ##moved from motifsearch      
        alignmentList=[]##returns list of seq aa that start motif
        length=len(seq)
        motifPosi=0##position of motif
        
        while motifPosi<=(length-self.length): ##check
            #print seq[motifPosi: motifPosi+self.length]
            if self.__eq__(seq[motifPosi: motifPosi+self.length]):
                alignmentList.append(motifPosi)
            motifPosi+=1 ##incremnts what motif is compared against
        return alignmentList    

###############  Amino acid position, can hold several possibilities, or x
class AA: 
    def __init__(self, alist=[]):
        self.alist=alist
        
    def __eq__(self, other):
        #print "custom eq executed"
        #print "possiblities are", self.alist
        if self.alist[0]=='x' or self.alist[0]=='X':
            return True
        for aa in self.alist:
            #print "for self: ", aa
            #print "other is ", other
            if aa.upper()==other.upper():
                return True
        return False
    def __ne__(self, other):
        return not self.__eq__(other)
        
    def __str__(self):
        list1=str(self.alist)
        return list1         
      
if __name__ == "__main__": main()
