import Bio 
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC


def main():
    ##read in name of file with protein and dna sequence, order- protein1 dna1  
    
    name=raw_input("type in name of file with protein and dna sequences in fasta format\n")
    infile=open(name, 'r')
    pNd=Bio.SeqIO.parse(infile, "fasta")
    pStartAA=raw_input("input first amino acid of desired protein sequence\n")
    pStartString=raw_input("input amino acids surrounding the first AA. DO NOT include repeats of that amino acid\n")
    pEndAA=raw_input("input last amino acid of desired protein sequence\n")
    pEndString=raw_input("input amino acids surrounding the last AA. DO NOT include repeats of that amino acid\n")
    ##get protein sequence
    EntireProtein=next(pNd).seq

    print EntireProtein


    ##find index of start and end position
    pStartIndex=str(EntireProtein).find(pStartString)+pStartString.find(pStartAA)
    pEndIndex=str(EntireProtein).find(pEndString)+pEndString.find(pEndAA)
    
    #print "index of pStartString is ",str(EntireProtein).find(pStartString)
    #print "index of pStartAA in PstartString is ", pStartString.find(pStartAA)
    print "location of pEndString is ", str(EntireProtein).find(pEndString)
    print "location of pEndAA is ", pEndString.find(pEndAA)
    print str(EntireProtein).find(pEndString)+pEndString.find(pEndAA)

    ##get dna sequence
    EntireDNA=str(next(pNd).seq).lower()
    dStartIndex=pStartIndex*3
    dEndIndex=pEndIndex*3+3 #note: this is end index, so one more than last base)
     
    print dEndIndex
    #print dStartIndex, dEndIndex

    forwardPrimerlist=[]
    reversePrimerlist=[]

    ForwardPrimerEnd=dStartIndex+3##TEST VALUE CHANGE TO 15!!
    while ForwardPrimerEnd<dStartIndex+40:##TEST VALUE CHANGE TO 30!!
        #print "while loop executed"
        #print EntireDNA[ForwardPrimerEnd].lower()
        if EntireDNA[ForwardPrimerEnd].lower()=='g':
            if EntireDNA[ForwardPrimerEnd-1].lower!='g':
                primerSeq=EntireDNA[dStartIndex:ForwardPrimerEnd+1]
                at=primerSeq.count('a')+primerSeq.count('t')
                cg=primerSeq.count('c')+primerSeq.count('g')
                meltingTemp=at*2+cg*4
                
                forwardPrimer=Primer(primerSeq, meltingTemp, at, cg)

                forwardPrimerlist.append(forwardPrimer)
        ForwardPrimerEnd+=1

    print "****************Forward Primers********************"
    for primer in forwardPrimerlist:
        print "Actual Forward primer:"
        seq=Seq(primer.seq, IUPAC.unambiguous_dna)
        complement=seq.complement()
        print complement
        primer.printf()
        print "*******"


    #make reverse primer
    reversePrimerBeg=dEndIndex-3
    while reversePrimerBeg>dEndIndex-30:
        if EntireDNA[reversePrimerBeg].lower()=='g':
            if EntireDNA[reversePrimerBeg+1].lower()!='g':
                primerSeq=EntireDNA[reversePrimerBeg:dEndIndex]

                at=primerSeq.count('a')+primerSeq.count('t')
                cg=primerSeq.count('c')+primerSeq.count('g')
                meltingTemp=at*2+cg*4
                reversePrimer=Primer(primerSeq, meltingTemp, at, cg)
                reversePrimerlist.append(reversePrimer)
        reversePrimerBeg-=1
    
    print "*****************Reverse Primers********************* "
    for primer in reversePrimerlist:
        primer.printf()
        print "Actual revrese primer:"
        seq=Seq(primer.seq, IUPAC.unambiguous_dna)
        reverseComplement=seq.reverse_complement()
        print reverseComplement
        print "******"


class Primer:
    def __init__(self, seq, melt, at, cg):
        self.seq=seq
        self.length=len(seq)
        self.meltingTemp=melt
        self.ATContent=at*100/self.length
        self.CGContent=cg*100/self.length
        
    
    def printf(self):
        print "seq:", self.seq
        print "melting temp:", self.meltingTemp
        print "atcontent: ", self.ATContent
        


 



if __name__ == "__main__": main()
