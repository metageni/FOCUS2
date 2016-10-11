#FOCUS2: agile and sensitive classification of metagenomics data using a reduced database | version 0.1 (2015)
#---------------------------------------------------------------------------------------------------------------------------------------
#(c)            Silva, G. G. Z.,B. E. Dutilh, and R. A. Edwards: 
#		FOCUS2: agile and sensitive classification of metagenomics data using a reduced database (not submitted)
#website: 	https://edwards.sdsu.edu/FOCUS2
#---------------------------------------------------------------------------------------------------------------------------------------

# python focus2r.py -q FASTA/FASTQ -b focus2_binning -dir output/
import os,sys,random

options= "FOCUS2(R)\n"\
      "--------------------------------------------------------------------------------------------------------------------\n"\
      "Options:\n"\
      "         -h          ------: print help\n"\
      "         -q          string: folder with multiple FASTA/FASTQ files\n"\
      "         -b          file: binning file for '-q' from FOCUS2\n"\
      "         -dir        string: output directory\n"\
      "         -mi         float:  minimum identity (default 60 %)\n"\
      "         -ml         int:    minimum alignment (default 45 nucleotides)\n"\
      "         -e          float:  e-value (default 0.00001)\n"\
      "         -t          int:    number of threads (default 1)\n"\
      "         -o          string: project name (default 'my_project')\n"\
      "--------------------------------------------------------------------------------------------------------------------\n"\
      "example> python focus2r.py -q FASTA/FASTQ -b focus2_binning -dir output/"

myproject={"-mi":"60","-ml":"45","-n":"0","-o":"myproject","-dir":"out","-a":"blastn","-e":"0.00001","-t":"1"}

#gets the user parameters and add in the hash
def setParameters():
    r=1
    if "/" in sys.argv[0]:
        myproject["dir"]="/".join(sys.argv[0].split("/")[:-1])+"/"
    else:
        myproject["dir"]=""
        
    userParameters=sys.argv[1:]

    if "-h" in sys.argv:
        pass
    else:
        for i in range(0,len(userParameters),2):
            try:
                myproject[userParameters[i]]=userParameters[i+1]
            except:
                if userParameters[i] in myproject:
                    print "Please inform a value for "+userParameters[i]
                    r=0
                else:
                    print userParameters[i]+" is not a valid parameter"
                    r=0
    return r


#returns the path for a given program name
def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

#check if Biopython is installed
def checkBio():
    n=0
    try:
        from Bio import SeqIO
        n=1
    except:
        pass
    return n

def get_unmapped():
    #get sequence IDs for the mapped
    o=open(myproject["-o"]+"__unmapped.fasta","w+")
    c=0
    for fasta_name in [i for i in os.listdir(myproject["-q"]) if i.split(".")[-1].lower() in ["fna","fasta"]]:
        h={}
        #fasta_name=.split("/")[-1]
        f=open(myproject["-b"])
        for line in f:
            if fasta_name == line.split("\t")[0]:
                h[line.split("\t")[1]]=0
        f.close()
        
        #create fasta with unmapped sequences
        from Bio import SeqIO
        

        
        for seq_record in SeqIO.parse(myproject["-q"]+"/"+fasta_name, "fasta"):
            if seq_record.id not in h:
                o.write(">"+fasta_name+"_unmapped___"+seq_record.id+"\n"+str(seq_record.seq)+"\n")
                c+=1
    o.close()

    return c

def align():
    fasta_name=myproject["-o"]+"__unmapped.fasta";fasta_name=fasta_name.split("/")[-1]
    os.system("hs-blastn align -db genomes/focus2r -window_masker_db genomes/mycount.counts.obinary -query "+myproject["-o"]+"__unmapped.fasta"+" -out "+myproject["-dir"]+"/"+fasta_name+".txt -outfmt 6 -evalue "+myproject["-e"]+" -num_threads "+myproject["-t"])
    

def get_taxon():
    f=open("genomes/taxon.xls")
    taxon={}
    for line in f:
        line=line.replace("\n","").replace("'","").replace('"',"").split("\t")
        taxon[line[-1].split("___")[-1]]="\t".join(line)
    f.close()
    return taxon

taxon=get_taxon()

#parses the alignment output and puts them into a hash (result)
def besthit():
    result={}
    fasta_name=myproject["-o"]+"__unmapped.fasta";fasta_name=fasta_name.split("/")[-1]
    f=open(myproject["-dir"]+"/"+fasta_name+".txt")

    h={}
    for line in f:
        line=line.split("\t")

        #if it passes in the min. identity and min. aligment
        if float(line[2])>=float(myproject["-mi"]) and int(line[3])>=int(myproject["-ml"]):
            sequenceID=line[0]
            subject=line[1].split("___")[0]
            evalue=line[-2]
            mi=line[2]
            ml=line[3]

            if sequenceID not in h:
                h[sequenceID]=[[evalue,mi,ml],[subject]]# 1st time the query is found
            else:
                #keeps only subjects that have the same evalue of the 1st hit (top-hit)
                if (evalue == h[sequenceID][0][0]) and (subject not in  h[sequenceID][1]):
                    h[sequenceID][1]+=[subject]
    f.close()

    result[fasta_name]=h
    return result

  
#Writes the binning for the results (sample/sequence_ID/taxa/evalue/ml/mi)
def bin_data(result):
    o=open(myproject["-dir"]+"/"+myproject["-o"]+"_focus2_binning.xls","w+")
    #o.write("Sample\tSequence_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tStrain\tE-value\n")
    for i in result:
        for j in result[i]:
            for mytaxon in result[i][j][1]:
                o.write(j.split("___")[0]+"\t"+j.split("___")[1]+"\t"+taxon[mytaxon]+"\t"+result[i][j][0][0]+"\n")
    o.close()

def join_bins():
    fasta_name=myproject["-o"]+"__unmapped.fasta";fasta_name=fasta_name.split("/")[-1]
    os.system("cat "+myproject["-b"]+" "+myproject["-dir"]+"/"+myproject["-o"]+"_focus2_binning.xls "+" > "+myproject["-dir"]+"/"+myproject["-b"].split("/")[-1]+"__FOCUS2R.xls")
    os.system("rm "+myproject["-dir"]+"/"+myproject["-o"]+"_focus2_binning.xls")
    os.system("rm "+myproject["-dir"]+"/"+fasta_name+".txt")
    os.system("rm "+fasta_name)
    
def main():
    print "1) Getting unmapped sequences"
    check=get_unmapped()
    if check ==0:
        print "All the sequences were mapped by FOCUS2 already :)"
    else:
        print "2) Aligning sequences to database"
        align()
        print "3) Parsing Alignment"
        result=besthit()
        bin_data(result)
        join_bins()
        print "\nDone :)\nPlease check the folder '"+myproject["-dir"]+"' for the output"
    
print "FOCUS2(R)"
p=0;error=""
if "-h" in sys.argv[1:]:
    print options
elif setParameters() == 1:
    if "-q" not in sys.argv[1:]:
        error+="-q     file: FASTA/FASTQ file\n";p+=1
    if "-b" not in sys.argv[1:]:
        error+="-b     file: binning file for '-b' from FOCUS2\n";p+=1        
    if "-dir" not in sys.argv[1:]:
        error+="-dir   string: output directory\n";p+=1
    if which("hs-blastn") == None:
        error+="-please install hsblastn\n";p+=1
    if checkBio()==0:
        error+="-please install Biopython\n";p+=1
        
    if p!=0:
        print "You are missing the following parameters:\n"+error
    else:
        outputfocus=myproject["-q"].split("/")
        if len(outputfocus[-1])==0:
            outputfocus=outputfocus[-2]
        else:
            outputfocus=outputfocus[-1]
        main()
