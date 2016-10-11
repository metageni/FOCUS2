#FOCUS2: agile and sensitive classification of metagenomics data using a reduced database | version 0.1 (2015)
#---------------------------------------------------------------------------------------------------------------------------------------
#(c)            Silva, G. G. Z.,B. E. Dutilh, and R. A. Edwards: 
#		FOCUS2: agile and sensitive classification of metagenomics data using a reduced database (not submitted)
#website: 	https://edwards.sdsu.edu/FOCUS2
#---------------------------------------------------------------------------------------------------------------------------------------

# python focus2.py -q input2/ -dir out2/
import os,sys,random

options= "FOCUS2: agile and sensitive classification of metagenomics data using a reduced database\n"\
      "--------------------------------------------------------------------------------------------------------------------\n"\
      "Options:\n"\
      "         -h          ------: print help\n"\
      "         -q          string: folder with multiple FASTA/FASTQ files\n"\
      "         -dir        string: output directory\n"\
      "         -o          string: project name (default 'my_project')\n"\
      "         -mi         float:  minimum identity (default 60 %)\n"\
      "         -ml         int:    minimum alignment (default 45 nucleotides)\n"\
      "         -k          int:    k-mer frequency used on FOCUS (default: 7) (6/7)\n"\
      "         -n          int:    normalize counts minimum alignment (0:False/1:True)(default: 0)\n"\
      "         -t          int:    number of threads (default 1)\n"\
      "         -e          float:  e-value (default 0.00001)\n"\
      "         -a          string: aligner (blastn/hsblastn) (default: hsblastn)\n"\
      "         -s          int:    split profiling in different levels (0:False/1:True)(default: 1)\n"\
      "         -bootstrap  int:    resamples the data to have more confidence in the results (0:False/1:True)(default: 0)\n"\
      "         -ns         int:    number of resampling per sample (default: 10)\n"\
      "         -b          float:  % of sequences to resample (default: 80.0)\n"\
      "--------------------------------------------------------------------------------------------------------------------\n"\
      "example> python focus2.py -q input/ -dir output/"

myproject={"-k":"7","-mi":"60","-ml":"45","-n":"0","-o":"myproject","-dir":"out","-a":"hsblastn","-s":"1","-e":"0.00001","-t":"1","-bootstrap":"0","-b":"80","-ns":"100"}

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

taxon={}
def run_focus_parse():
    os.system("python focus/focus.py -s 1 -q "+myproject["-q"]+" -k "+myproject["-k"])

    f=open(outputfocus+"__STAMP_tabular.spf")
    #head of the focus output which contains the file names (we need it later for the alignment process)
    head=f.readline().replace("\n","").replace("\r","")
    head=head.split("\t")[8:]

    #hash where keys are the IDS for the analyzed samples and values are lists with the focus prediction
    samples_prediction={x:[] for x in head}
    
    #reads the other lines on the focus output
    for line in f:
        line=line.replace("\n","").replace("\r","").replace('"',"")
        
        genomeID=line.split("\t")[7].split("___")[-1]
        #save the taxonomy of each focus prediction for later
        taxon[genomeID]="\t".join(line.split("\t")[:8])
        
        line=line.split("\t")[7:]#strain name + its profile for all the samples
        
        c=0
        #stores the profile for genome X in sample Y
        for j in line[1:]:
            if float(j)!=0:
                samples_prediction[head[c]]+=[genomeID+".fna"]
            c+=1
    f.close()

    return samples_prediction

def run_subsample():
    def subsample(fasta,number_of_samples):
        myfile="".join(file(fasta).readlines()).split(">")[1:]
        number_sequences=len(myfile)
        
        for ID in xrange(1,number_of_samples+1):
            o=open(fasta+"___subsample_"+str(ID)+".fasta","w+")
            {o.write(">"+x):0 for x in random.sample(myfile,int(number_sequences*(float(myproject["-b"])/100)))}
            o.close()
            
    for myfile in [i for i in os.listdir(myproject["-q"]) if i.split(".")[-1].lower() in ["fna","fasta"]]:
        path=myproject["-q"]+"/"
        print "Resampling "+myfile+" "+myproject["-ns"]+" times"
        subsample(path+myfile,int(myproject["-ns"]))

def run_focus_bootstrap():
    run_subsample()
    os.system("python focus/focus.py -bootstrap 1 -s 1 -q "+myproject["-q"]+" -k "+myproject["-k"])
    os.system("rm "+myproject["-q"]+"/*___subsample_*.fasta")
    #create a hash with taxa -> samples -> abundance
    sub_profile={}
    f=open(outputfocus+"__STAMP_tabular.spf")
    #original metagenome ids
    IDS=f.readline().replace("\n","").split("\t")[8:]
    for line in f:
        line=line.replace("\n","").split("\t")
        taxa=line[7].split("___")[-1]
        taxon[taxa]="\t".join(line[:8])
        taxa=taxa+".fna"#we only care about the strain because it is the PK
        
        #we see the taxa and create a hash as value using the IDs
        sub_profile[taxa]={x.split("___subsample_")[0]:[] for x in IDS}

        #abundances of "taxa" to all metagenomes
        abundance=line[8:]
        #create a hash with taxa -> samples -> abundance
        for i in xrange(len(abundance)):
            meta=IDS[i].split("___subsample_")[0]
            sub_profile[taxa][meta]+=[abundance[i]]
    f.close()

    #creates a hash with the genomes that are present in at least 50% of the resamples
    corrected_taxa={sample.split("___subsample_")[0]:[] for sample in IDS}
    for taxa in sub_profile:
        for sample in sub_profile[taxa]:
            info=sub_profile[taxa][sample]
            if (info.count("0.0")/(len(info)*1.))<=0.2:
                corrected_taxa[sample]+=[taxa]
    os.system(outputfocus+"__STAMP_tabular.spf")
    
    return corrected_taxa
                
def produce_alignments(samples_prediction):
    path_genomes=" genomes/"
    #reads hash with the profiles to make the DB and run blastn/HS-blastn
    for fasta_name in samples_prediction:
        print "Processing:     "+fasta_name
        filename=myproject["-o"]+"__STAMP_tabular.spf__"+fasta_name

        temp=path_genomes.join(samples_prediction[fasta_name])
        #joins all the files for the sample "i"
        os.system("cat "+path_genomes+temp+" > db/"+filename)

        #formats the database
        os.system("makeblastdb -in db/"+filename+" -dbtype nucl  -out db/"+filename+".db")

        ##### Homology search
        if myproject["-a"].lower() == "blastn":
            #Running blastn
            os.system("blastn -db db/"+filename+".db -query "+myproject["-q"]+"/"+fasta_name+" -out "+myproject["-dir"]+"/"+fasta_name+".txt -outfmt 6 -evalue "+myproject["-e"]+" -num_threads "+myproject["-t"])
            
        elif myproject["-a"].lower()=="hsblastn":
            os.system("windowmasker -in db/"+filename+".db -infmt blastdb -mk_counts -out mycount.counts")
            os.system("windowmasker -in mycount.counts -sformat obinary -out mycount.counts.obinary -convert")
            os.system("hs-blastn index db/"+filename)
            os.system("hs-blastn align -db db/"+filename+" -window_masker_db mycount.counts.obinary -query "+myproject["-q"]+"/"+fasta_name+" -out "+myproject["-dir"]+"/"+fasta_name+".txt -outfmt 6 -evalue "+myproject["-e"]+" -num_threads "+myproject["-t"])
            os.system("rm mycount.counts*")
        
        #deletes all the db files related to the query
        os.system("rm db/"+filename+"* -r")

#parses the alignment output and puts them into a hash (result)
def besthit():
    result={}
    for myfile in [i for i in os.listdir(myproject["-dir"]) if i.split(".")[-1]=="txt"]:
        f=open(myproject["-dir"]+"/"+myfile)

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

        result[myfile]=h
    return result #this hash will have the result for all the FASTA/FASTQ files entered

#Writes the profiling for the results (counts/relative abundance)
def profile_data(result):
    k=result.keys()
    #writes profiling table with counts and relative abundance
    taxa={};c=0

    for i in k:
        for j in result[i]:
            L=len(result[i][j][1])
            for mytaxa in result[i][j][1]:
                if mytaxa not in taxa:
                    taxa[mytaxa]=[0]*len(k)
                if myproject["-n"] == 1:
                    taxa[mytaxa][c]+=1./L
                else:
                    taxa[mytaxa][c]+=1
        c+=1
    total_sum=[0.0]*len(k)

    for i in range(len(k)):
        for j in taxa:
            total_sum[i]+=taxa[j][i]
        
    o=open(myproject["-dir"]+"/"+myproject["-o"]+"_profiling.xls","w+")
    o.write("Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tStrain\t"+"\t".join(["Counts: "+xx for xx in k])+"\t"+"\t".join(["Relative Abundance (%): "+xx for xx in k])+"\n")

    for i in taxa:
        relative=[str((taxa[i][xx]/total_sum[xx])*100) for xx in range(len(total_sum))]
        o.write(taxon[i]+"\t"+"\t".join([str(xx) for xx in taxa[i]])+"\t"+"\t".join(relative)+"\n")
    o.close()
    
#Writes the binning for the results (sample/sequence_ID/taxa/evalue/ml/mi)
def bin_data(result):
    o=open(myproject["-dir"]+"/"+myproject["-o"]+"_binning.xls","w+")
    o.write("Sample\tSequence_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tStrain\tE-value\n")
    for i in result:
        for j in result[i]:
            for mytaxon in result[i][j][1]:
                o.write(i.replace(".txt","")+"\t"+j+"\t"+taxon[mytaxon]+"\t"+result[i][j][0][0]+"\n")
    o.close()

#gets the STAMP output and split by level
def split_stamp(input_file):
    from numpy import sum as SUM
    choice=range(7)
    for p in choice:
        h={}

        f=open(input_file)
        line=f.readline()
        head=line.split()[p]
        
        info=line.replace("\n","").split("\t")[8:]
        c=0
        for line in f:
            line=line.split("\t")
            temp_head=line[p]
            temp_info=line[8:]

            if temp_head not in h:
                h[temp_head]=[]
            h[temp_head]+=[[float(x) for x in temp_info]]
            
        f.close()

        o=open(myproject["-dir"]+"/"+head+"__"+myproject["-o"]+"_profiling.xls","w+")
        o.write(head+"\t"+"\t".join(info)+"\n")
        for i in h:
            o.write(i+"\t"+"\t".join([str(x) for x in list(SUM(h[i],axis=0))])+"\n")
        o.close()
        
def main():
    if myproject["-bootstrap"]=="0":
        print "1) Running FOCUS"
        samples_prediction=run_focus_parse()
    elif myproject["-bootstrap"]=="1":
        samples_prediction=run_focus_bootstrap()
    print "2) Aligning input data to database"
    produce_alignments(samples_prediction)
    print "3) Parsing Alignments"
    result=besthit()
    profile_data(result)
    os.system("rm "+myproject["-dir"]+"/*.txt")
    #splits the profile result creating one file per level
    if int(myproject["-s"])==1:
        split_stamp(myproject["-dir"]+"/"+myproject["-o"]+"_profiling.xls")
    bin_data(result)
    #os.system("rm *__STAMP_tabular.spf")
    print "\nDone :)\nPlease check the folder '"+myproject["-dir"]+"' for the output"

print "FOCUS2: agile and sensitive classification of metagenomics data using a reduced database"
p=0;error=""
if "-h" in sys.argv[1:]:
    print options
elif setParameters() == 1:
    if "-q" not in sys.argv[1:]:
        error+="-q     string: folder with multiple FASTA/FASTQ files\n";p+=1
    if "-dir" not in sys.argv[1:]:
        error+="-dir   string: output directory\n";p+=1
    if which("blastn") == None:
        error+="-please install blastn\n";p+=1
    if myproject["-a"]=="hsblastn" and which("hs-blastn") == None:
        error+="-please install hsblastn\n";p+=1
    if myproject["-a"]  not in ["blastn","hsblastn"]:
        error+="-a     string: aligner (blastn/hsblastn) (default: blastn)";p+=1    
    if p!=0:
        print "You are missing the following parameters:\n"+error
    else:
        outputfocus=myproject["-q"].split("/")
        if len(outputfocus[-1])==0:
            outputfocus=outputfocus[-2]
        else:
            outputfocus=outputfocus[-1]
        main()
