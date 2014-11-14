import sys
import string
import gzip
import os
#USAGE: python <programName> <file1Metadata> <file2Metadata>
#INPUT: the metadata for the 2 files to be compared for duplicates (%IncorrectCalls per marker)
#ASSUMPTIONS: both files must be sorted by chr:pos. file2 marker positions must be equal or greater than file1 positions.
#OUTPUT: one file for each input file specifying PASS/FAIL status for each row/marker

if sys.argv[1].endswith("gz"):
	file1 = gzip.open(sys.argv[1],"r")
else:
	file1 = open(sys.argv[1],"r")

if sys.argv[2].endswith("gz"):
	file2 = gzip.open(sys.argv[2],"r")
else:
	file2 = open(sys.argv[2],"r")

if os.path.exists(sys.argv[1]+".out"):
	outFile1 = open(sys.argv[1]+".out","r")
outFile1New = open(sys.argv[1]+".out.new","w")

if os.path.exists(sys.argv[2]+".out"):
	outFile2 = open(sys.argv[2]+".out","r")
outFile2New = open(sys.argv[2]+".out.new","w")

def readFile(fileHandle):
	line = fileHandle.readline()
	columns=[]
	chr=""
	pos=""
	if line<>"":
		columns = line.strip().split()
		if columns[0][3:].upper() == "X":
        		chr = 23
	        elif columns[0][3:].upper() == "Y":
        	        chr = 24
	        elif columns[0][3:].upper() == "M":
        	        chr = 25
	        else:
        	        chr = int(columns[0][3:])
	        pos = int(columns[1])
	
	return [chr,pos,columns]

#skip header lines on both files
file1.readline()
file2.readline()	

#read first lines from both files
[chr1,pos1,columns1] = readFile(file1)
metaChr1=""
metaPos1=""
metaColumns1=[]
if os.path.exists(sys.argv[1]+".out"):
	[metaChr1,metaPos1,metaColumns1] = readFile(outFile1)

[chr2,pos2,columns2] = readFile(file2)
metaChr2=""
metaPos2=""
metaColumns2=[]
if os.path.exists(sys.argv[2]+".out"):
	[metaChr2,metaPos2,metaColumns2] = readFile(outFile2)

while chr1<>"" and chr2<>"":
	print chr1,pos1,chr2,pos2
	if chr1<chr2  or (chr1==chr2 and pos1<pos2):
		if not metaColumns1 or (metaChr1==chr1 and metaPos1==pos1 and metaColumns1[2]<>"FAIL"):  #if this region has been processed before, FAILED marker should not be marked as PASS again
			print >>outFile1New, '\t'.join([columns1[0],columns1[1],"PASS"])
		elif metaChr1<>chr1 or metaPos1<>pos1:
			print "Metafile does not match parent file. ",chr1,pos1,metaChr1,metaPos1
			sys.exit(1)
		else:
			print >>outFile1New, '\t'.join(metaColumns1)
		[chr1,pos1,columns1] = readFile(file1)
		metaChr1=""
		metaPos1=""
		metaColumns1=[]
		if os.path.exists(sys.argv[1]+".out"):
		        [metaChr1,metaPos1,metaColumns1] = readFile(outFile1)
	elif chr1==chr2 and pos1==pos2:
		if columns1[5]<=columns2[5]:
			if not metaColumns1 or (metaChr1==chr1 and metaPos1==pos1 and metaColumns1[2]<>"FAIL"):
				print >>outFile1New, '\t'.join([columns1[0],columns1[1],"PASS"])	
			elif metaChr1<>chr1 or metaPos1<>pos1:
	                        print "Metafile does not match parent file. ",chr1,pos1,metaChr1,metaPos1
        	                sys.exit(1)
                	else:
	                        print >>outFile1New, '\t'.join(metaColumns1)


			print >>outFile2New, '\t'.join([columns2[0],columns2[1],"FAIL"])
		else:
			print >>outFile1New, '\t'.join([columns1[0],columns1[1],"FAIL"])
                        
		
			if not metaColumns2 or (metaChr2==chr2 and metaPos2==pos2 and metaColumns2[2]<>"FAIL"):
                                print >>outFile2New, '\t'.join([columns2[0],columns2[1],"PASS"])
                        elif metaChr2<>chr2 or metaPos2<>pos2:
                                print "Metafile does not match parent file. ",chr2,pos2,metaChr2,metaPos2
                                sys.exit(1)
                        else:
                                print >>outFile2New, '\t'.join(metaColumns2)
		
		[chr1,pos1,columns1] = readFile(file1)
		metaChr1=""
                metaPos1=""
                metaColumns1=[]
                if os.path.exists(sys.argv[1]+".out"):
                	[metaChr1,metaPos1,metaColumns1] = readFile(outFile1)

		[chr2,pos2,columns2] = readFile(file2)
		metaChr2=""
		metaPos2=""
		metaColumns2=[]
		if os.path.exists(sys.argv[2]+".out"):
			[metaChr2,metaPos2,metaColumns2] = readFile(outFile2)
	elif chr1>chr2 or (chr1==chr2 and pos1>pos2):
		if not metaColumns2 or (metaChr2==chr2 and metaPos2==pos2 and metaColumns2[2]<>"FAIL"):
                	print >>outFile2New, '\t'.join([columns2[0],columns2[1],"PASS"])
                elif metaChr2<>chr2 or metaPos2<>pos2:
                	print "Metafile does not match parent file. ",chr2,pos2,metaChr2,metaPos2
                        sys.exit(1)
                else:
                        print >>outFile2New, '\t'.join(metaColumns2)
		
		[chr2,pos2,columns2] = readFile(file2)	
		metaChr2=""
                metaPos2=""
                metaColumns2=[]
                if os.path.exists(sys.argv[2]+".out"):
                	[metaChr2,metaPos2,metaColumns2] = readFile(outFile2)

#print leftover rows from each file
while chr1<>"":
	if os.path.exists(sys.argv[1]+".out"):
		print >>outFile1New, '\t'.join(metaColumns1)
	else:
		print >>outFile1New, '\t'.join([columns1[0],columns1[1],"PASS"])
        [chr1,pos1,columns1] = readFile(file1)				
	metaChr1=""
        metaPos1=""
        metaColumns1=[]
        if os.path.exists(sys.argv[1]+".out"):
        	[metaChr1,metaPos1,metaColumns1] = readFile(outFile1)
while chr2<>"":
	if os.path.exists(sys.argv[2]+".out"):
                print >>outFile2New, '\t'.join(metaColumns2)
        else:
                print >>outFile2New, '\t'.join([columns2[0],columns2[1],"PASS"])
        [chr2,pos2,columns2] = readFile(file2)
        metaChr2=""
        metaPos2=""
        metaColumns2=[]
        if os.path.exists(sys.argv[2]+".out"):
        	[metaChr2,metaPos2,metaColumns2] = readFile(outFile2)
	

file1.close()
file2.close()

outFile1New.close()
outFile2New.close()

os.rename(sys.argv[1]+".out.new",sys.argv[1]+".out")
os.rename(sys.argv[2]+".out.new",sys.argv[2]+".out")
