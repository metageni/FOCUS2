# FOCUS2 0.1 - This program is used to download and format the databases for the tool (2015)
import os

print "1. Downloading the FOCUS2 database (PATRIC genomes + k-mer count for genomes)"
c=1
try:
    os.system("wget http://edwards.sdsu.edu/FOCUS2/db.tar.gz")
except:
    c=0
    print "There is a problem in your internet connection!"
if c==1:
    print "2. Uncompressing database"
    os.system("tar -xvf db.tar.gz") #uncompress db
    os.system("mv k6.pickle focus/db/")
    os.system("mv k7.pickle focus/db/")
    os.system("rm db.tar.gz")
    print "3. Done :)"
