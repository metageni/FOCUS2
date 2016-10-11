#FOCUS2
Agile and sensitive classification of metagenomics data using a reduced database || version 0.1

(c)            Silva, G. G. Z., B. E. Dutilh, and R. A. Edwards: 
		FOCUS2: Agile and sensitive classification of metagenomics data using a reduced database (unplublished).
website: 	https://edwards.sdsu.edu/FOCUS2


Program
--------
focus2__downloadDB.py.py: Downloads the FOCUS2 database

(1) USAGE
python focus2__downloadDB.py

Program
--------
focus2.py: FOCUS2 main program

(1) USAGE
-----

Options:
	-h  none
		print help
		
	-q string
		folder with multiple FASTA/FASTQ files
		
	-dir string
		output directory
		
	-o string
		project name (default 'my_project')
		
	-mi float
		minimum identity (default 60 %)
		
	-ml int
		minimum alignment (default 45 nucleotides)
		
	-k int
		k-mer frequency used on FOCUS (default: 7) (6/7)
		
	-n int
		normalize counts minimum alignment (0:False/1:True)(default: 0)
		
	-t int
		number of threads (default 1)
		
	-e float
		e-value (default 0.00001)
		
	-a string
		aligner (blastn/hsblastn) (default: hsblastn)
		
	-s int
		split profiling in different levels (0:False/1:True)(default: 1)
		
	-bootstrap int
		resamples the data to have more confidence in the results (0:False/1:True)(default: 0)
		
	-ns int
		number of resampling per sample (default: 10)
		
	-b float
		% of sequences to resample (default: 80.0)
		
	example> python focus2.py -q input/ -dir output/
	 
(2) OUTPUT
FOCUS2 output will be add the folder selected in -dir

Program
--------
FOCUS2(R)

Options:
	-h none
		print help
		
	-q string
		folder with multiple FASTA/FASTQ files
		
	-b file
		binning file for '-q' from FOCUS2
		
	-dir string
		output directory
		
	-mi float
		minimum identity (default 60 %)
		
	-ml int
		minimum alignment (default 45 nucleotides)
		
	-e float
		e-value (default 0.00001)
		
	-t int
		number of threads (default 1)
		
	-o string
		project name (default 'my_project')
	
	example> python focus2r.py -q FASTA/FASTQ -b focus2_binning -dir output/

(2) OUTPUT
FOCUS2R output will be add the folder selected in -dir

(3) OBSERVATION
- FOCUS2(R) only uses hs-blastn as aligner

DEPENDENCIES
------------
- Python >= 2.6.X < 3.Y: http://www.python.org/download
- Jellyfish: http://www.cbcb.umd.edu/software/jellyfish
- Numpy: http://sourceforge.net/projects/numpy/files/NumPy
- SciPy: http://sourceforge.net/projects/scipy

One of the below aligners:
- BLAST: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST
- HS-BLASTN: https://github.com/chenying2016/queries

COPYRIGHT AND LICENSE
---------------------
Copyright (C) 2015-2016  Genivaldo Gueiros Z. Silva

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.