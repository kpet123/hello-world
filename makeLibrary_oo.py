###############
# makes library using the HitQuery object
# writeLib can be altered to specifications
# writeFasta() makes a fasta file that can then be manipulated: note that these are not full seq!
##############
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast import Record
import Pattern as p


def main():
        
        print "Options: 1 for new blast search, 2 for pre-made .xml file. If you want to save copy of new blast search, uncomment beginning area to save and save_file.close"
        option=raw_input()
        if option=='1':  
		print "Type name of input .fasta file"
                filename=raw_input()
                #opens .fasta file and reads as
		fasta_string=open(filename).read()
		#creates a handle object in XML 
		print "about to run blast"
		result_handle=NCBIWWW.qblast("blastp", "nr", fasta_string, None, None, None, None, '(none)', 20.0, None, None, None, 10000)
                
                print "Save this file? (type y or n)"
                answer=raw_input()
                if answer=='y':
                	#save local copy of output file
                        print "Type name of saved file (ends with .xml)"
                        savename=raw_input()
                	save_file=open(savename, "w")
         		save_file.write(result_handle.read())
         		save_file.close()
        		result_handle.seek(0)#brings 'cursor' back to beginning

        if option=='2':
        	print "Input name of .xml file to be read"
        	filename=raw_input()
        	result_handle = open(filename)
                
	#use read for 1 query search
	#blast_record=NCBIXML.read(result_handle)
	#print(type(blast_record))

#manipulates blast data and returns alignement
         
        #parse into a bunch of blast records: 
        #flavor text: Module NCBIXML contains BlastParser, which we use here: this parses XML BLAST data into a Record.Blast object.
        
	blast_records=NCBIXML.parse(result_handle)
	

	#variables referenced in entire loop
 	queryDescrList=p.listDef("./Folder-FASTA/CnaBquery.fasta")
	iterator=0
        hitQueryList=[]##holds HitQuery objects
	eboxlist=open("CnaBorder.txt", "r")#returns real eboxes in list

        
        #iterates through a series of Record.Blast objects and creates list of HitQuerys
	for blast_record in blast_records:#layer 1, different comparison proteins
	    ###progress report	    

	    absoluteEboxLocation=int(eboxlist.readline())+1#check mechanism
	    queryName=queryDescrList[iterator]
	    iterator+=1
            print "Making library for "+queryName
 	    newQuery=HitQuery(queryName, absoluteEboxLocation, blast_record)
            print "     Number of hits: "+str(len(newQuery.hitlist))
	   
            hitQueryList.append(newQuery)
            
        
       
        writeLib(hitQueryList)
  


    	    
#### Writes to file
def writeLib(hitQueryList):
	print "input name of library file to be created"
	libname=raw_input()
	f=open(libname, 'w')


	for query in hitQueryList:
	    print query.queryName
	    for hit in query.hitlist:
		print "hit written"
		entry=hit.writeFasta() 
		f.write(entry)
		#f.write("\n\n\n")
	f.close()



#this is own alignment class, holds a list of additional stuff to Biopython class
#eventually modify with inheritance?

class HitQuery:
    def __init__(self, queryName,  absEboxLocation, blast_record):
        self.queryName=queryName
        #self.querySeq=querySeq
	self.hitlist=[]
        self.blast_record=blast_record
        self.absEboxLoc=absEboxLocation
	self.fillList(self.hitlist)

    def fillList(self, hitlist):
	for alignment in self.blast_record.alignments:
	    title=alignment.title
	    #print "********Alignment**********", alignment.title
	    for hsp in alignment.hsps:#hsp is the specific alignment section
		queryList=p.CnaBsearch(hsp.query)#eboxes in hsp.query
                sbjctList=p.CnaBsearch(hsp.sbjct)#eboxes in hsp.subject
		for count in queryList:
		     #print count
		     if (count in sbjctList) and (hsp.query_start+count==self.absEboxLoc):
	                 hit=Hit(title, self.queryName, hsp, count)
                         hitlist.append(hit)

     

class Hit: 
    def __init__(self, alignTitle, queryTitle,  hsp, eboxLoc):
        self.hsp=hsp
        self.eboxLoc=eboxLoc
	self.queryTitle=queryTitle
	self.description=alignTitle
	self.gi=alignTitle[15]#depends on correct format
	self.species=self.parseSpecies(alignTitle)#depends on correct format
        

    ####shows first 150 alignments in hsp
    def showtotalAlignment(self):
	query150="Query: "+str(self.hsp.query[0:150])+"\n"
	sbjct150="Sbjct: "+str(self.hsp.sbjct[0:150])
        return query150+sbjct150
    
    ####### returns subjct sequence
    def returnHitSeq():
        return self.hsp.sbjct

    ####shows ebox plus 15 aa
    def showEboxAlignment(self):
	 queryEbox=self.hsp.query[(self.eboxLoc-4):(self.eboxLoc)+15]
         subjectEbox=self.hsp.sbjct[(self.eboxLoc-4):(self.eboxLoc+15)]
	 eboxAlign="Query Ebox: "+queryEbox+"\n"+"Sbjct Ebox: "+subjectEbox
         return eboxAlign
    
    def parseSpecies(self, title):
	return p.findEnclosed(title, "[", "]")

    def writeFasta(self):
        header= ">"+ self.description+"\n"
        seq=self.hsp.sbjct+"\n\n"
        return str(header+seq)
    ####desire output for library, for query information, use in conjuction with QueryHits write() function
    ####FORMAT: str(attribute)+"\n" 
    def write(self):
        
        ######STANDARD LIB###########
	descr="Description-  \n"+str(self.description)+"\n"
        qDescr="Query-  \n"+ self.queryTitle+'\n'
	eLoc="Ebox location- \n "+str(self.eboxLoc)+'\n'
	ealign="Ebox alignment lignment- \n"+self.showEboxAlignment()+"\n"
	talign="Total alignemnt- \n"+self.showtotalAlignment()+"\n"
	return str(descr+qDescr+eLoc+ealign+talign)
       


	

          



if __name__ == "__main__": main()


