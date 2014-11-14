import sys
import string
import gzip
import math

#filters CMS regions, extreme heterozygous positions, multi-allelic markers, <90% call rate markers,zero variance markers
#changes MIE, VQLOW and half calls to missing. changes haploid no calls on autosomal chromosomes to diploid calls
#removes phase - replace('|','/')
#NOTE: markers with missing ALT are not filtered because we don't know if the ALT is unknown or if it is a deletion

#INPUT: 
#-vcf=vcf file path, REQUIRED
#-out=output file path, REQUIRED
#-cms=cms file path; if omitted, CMS regions will not be filtered OPTIONAL
#-exthet=extreme heterozygous file path; if omitted, extreme heterozygous positions will not be filtered OPTIONAL
#-sampleStart = column no. at which sample data starts (0 based indexing)
#-maf= minimum minor allele frequency to filter on; if omitted no filtering is done based on maf OPTIONAL
#-callrate= minimum call rate to filter on OPTIONAL
#--biallelic : include bi-allelic markers only. If omitted, multi-allelic markers are included too
#--nonzeroVariance : include only those markers that have non -zero variance
#-exclude=file path containing list of markers to exclude (OPTIONAL)
#--ppc

VCF_FILENAME = ""
OUT_FILENAME = ""
CMS_FILENAME = ""
EXTHET_FILENAME = ""
MAF = 0
CALLRATE = 0
BIALLELIC = 0
NONZEROVARIANCE = 0
EXCLUDE_FILENAME = ""
PPC=0
SAMPLESTART = 0

def readInputArguments():
        global VCF_FILENAME
        global OUT_FILENAME
        global CMS_FILENAME
        global EXTHET_FILENAME
        global MAF
        global CALLRATE
 	global BIALLELIC
	global NONZEROVARIANCE
	global EXCLUDE_FILENAME
	global PPC
	global SAMPLESTART
       
	#PARSE COMMAND LINE ARGS. (2 of the above arguments are required)
        assert (len(sys.argv) >=3 ), 'Insufficient number of arguments.'

        while len(sys.argv) > 1:
                # Command line arguments containing '=' are implicitly options.
                thisArg = sys.argv.pop(1)
		print thisArg
		if not thisArg.startswith('--'):
			name,value = thisArg.split("=")  #split name value pairs
                else:
			name = thisArg
			value = ""
		if __debug__: # ...just to see what's going on.
                	print ("{},{}".format( name, value ))
		name = name.lower().strip("- ")  #strip hyphens and white spaces
                if name == "vcf":
                	VCF_FILENAME = value.strip(" ")
                elif name == "out":
	                OUT_FILENAME = value.strip(" ")
		elif name == "cms":
                        CMS_FILENAME = value.strip(" ")
		elif name == "exthet":
                        EXTHET_FILENAME = value.strip(" ")
		elif name == "samplestart":
			SAMPLESTART = int(value.lower().strip(" "))
                elif name == "maf":
                        MAF = float(value.lower().strip(" "))
                elif name == "callrate":
                        CALLRATE = float(value.lower().strip(" "))
                elif name == "biallelic":
                        BIALLELIC = 1
		elif name == "nonzerovariance":
			NONZEROVARIANCE = 1		
		elif name == "exclude":
			EXCLUDE_FILENAME = value.strip(" ")
		elif name == "ppc":
			PPC = 1
                else:
                	print "unrecognized option:",name
                        sys.exit(1)

        assert (VCF_FILENAME <> ""), 'VCF file path was not provided'
        assert (OUT_FILENAME <> ""), 'Output file path was not provided'


def isCMS(coord):
	global currentCMS
	global cmsChrNo
	if cmsChrNo == 'X':
		cmsChrNo = 23
	else:
		cmsChrNo = int(cmsChrNo)
	vcfChrNo = coord.split(":")[0][3:]
	vcfBP = int(coord.split(":")[1])	
	if vcfChrNo == 'X':
		vcfChrNo = 23
	elif vcfChrNo == 'Y':
		vcfChrNo = 24
	elif vcfChrNo == 'M':
		vcfChrNo = 25
	else:
		vcfChrNo = int(vcfChrNo)
		
	while len(currentCMS) > 1 and ((vcfChrNo == cmsChrNo and vcfBP > int(currentCMS[2])) or (vcfChrNo > cmsChrNo)):
		currentCMS = cmsFile.readline().strip().split('\t')
		if len(currentCMS) > 1:   #not end of file
			cmsChrNo = currentCMS[0][3:]
			if cmsChrNo == 'X':
				cmsChrNo = 23
			else:
				cmsChrNo = int(cmsChrNo)

		#vcfPos is <= currentCMS end pos; check if its also greater than currentCMSstartPos
	if len(currentCMS) > 1 and vcfBP >= int(currentCMS[1]):
		return 1
	else:
		return 0


def isAutosomal(chrStr):
	if chrStr in ['chrX','chrY','chrM','chr23','chr24','chr25','23','24','25']:
		return 0
	else:
		return 1

readInputArguments()

#change file paths if required
vcfFile = gzip.open(VCF_FILENAME,"r")
outputFile = gzip.open(OUT_FILENAME,"w")
logFile = gzip.open(OUT_FILENAME+".log.gz","w")
if CMS_FILENAME <> "":
	cmsFile = open(CMS_FILENAME,"r")
	currentCMS = cmsFile.readline().strip().split('\t')
	cmsChrNo = currentCMS[0][3:]
if EXTHET_FILENAME <> "":
	extHetFile = open(EXTHET_FILENAME,"r")
	extHetPos = extHetFile.read().strip().split()
if EXCLUDE_FILENAME <> "":
	excludeFile = open(EXCLUDE_FILENAME,"r")
	excludeList = excludeFile.read().strip().split()

altIndex = -1
for line in vcfFile:
	boolExtHet = 0
	boolCMS = 0
	boolExclude = 0
	boolNonPPC = 0
	multiallelic = 0
	missingCount = 0	
	altSet = set()
	altCount = 0
	refCount = 0	
	thisMaf = 0.0
	columns = line.strip().split('\t')
	if columns[0].find('#CHROM') == -1 and len(columns) > 1:   #data row

		coord = columns[0]+":"+columns[1]
		#filter flags for reuse	
		if extHetFile and coord in extHetPos:
			boolExtHet = 1
	
		if cmsFile and isCMS(coord):
			boolCMS = 1
	
		if excludeFile and coord in excludeList:
                        boolExclude = 1

                if BIALLELIC and columns[altIndex].find(',') <> -1:
                        multiallelic = 1
	
		if PPC and columns[6].find('PPC') == -1:
			boolNonPPC = 1

		assert(altIndex <> -1),"ALT Index = -1!"
		print boolNonPPC, columns[6]
		if boolExtHet == 0 and boolCMS == 0 and multiallelic == 0 and boolExclude == 0 and boolNonPPC == 0: 	#universally heterozygous coordinates and CMS and multi-allelic markers 
			for index in range(SAMPLESTART,len(columns)):	#for each sample genotype
				if columns[index].find("VQLOW") <> -1 or columns[index].find("MIE") <> -1: #VQLOW and MIE = missing
					columns[index] = columns[index].split(":")[0].replace('|','/')
					if isAutosomal(columns[0]): #force these to be diploid calls. mkvcf reports autosomal no calls as haploid sometimes.
						columns[index] = './.'
					elif columns[index].find('/')<>-1:#sex chromosomes, keep ploidy as is
						columns[index] = './.'
					else:
						columns[index] = '.'	
					missingCount += 1
					altSet.add('./.')
				else:
					columns[index] = columns[index].split(":")[0].replace('|','/')	#retain only the genotype field
					if columns[index].find(".") <> -1:   #even one missing allele(half call) changes to missing genotype
						if isAutosomal(columns[0]):
							columns[index] = './.'
						elif columns[index].find('/')<>-1:
							columns[index] = './.' 
						else:
							columns[index] = '.'
						missingCount += 1
						altSet.add('./.')
					#ADDED: Oct 15, 2013
					else:
						alleles = columns[index].split('/') #full valid call; see if it has non-ref allele
						altCount += sum(1 for x in alleles if x <> '0')
						refCount += sum(1 for x in alleles if x == '0')
						altSet.add(columns[index])
						
			missingRate = missingCount/float(len(columns[SAMPLESTART:len(columns)]))
			if missingRate < 1:	
				if altCount < refCount:
					thisMaf = altCount/float(altCount+refCount)
				else:
					thisMaf = refCount/float(altCount+refCount)
			if (CALLRATE == 0 or (CALLRATE <> 0 and missingRate <= (1-CALLRATE))) and (MAF==0 or (MAF <> 0 and thisMaf >= MAF)) and (NONZEROVARIANCE==0 or (NONZEROVARIANCE <> 0 and len(altSet)>1)):   
				outputStr = '\t'.join(columns)+'\n'
				outputFile.write(outputStr)
				
			else: 
				logText = '\t'.join(['CallRate/MAF/Variance',coord,'MissingRate',str(missingRate),'MAF',str(thisMaf),'ValueSet',str(altSet)])+'\n'
				logFile.write(logText)
		else:
			logText=''
			if boolExtHet == 1:
				logText = 'ExtremeHeterozygous'+'\t'
			if boolCMS == 1:	
				logText = logText+'CMS'+'\t'
			if multiallelic == 1:
				logText = logText+'Multi-allelic'+'\t'
			if boolExclude == 1:
				logText = logText+'ExcludedDuplicate'+'\t'
			if boolNonPPC == 1:
				logText = logText+'NonPPC'+'\t'
			logText=logText+coord+'\n'
			logFile.write(logText)	 
						
	else: #output header as is
		if columns[0].find('#CHROM') <> -1:
			columns = line.strip().split('\t')
       		        outputFile.write(line)
			altIndex = columns.index('ALT')
		else:
			outputFile.write(line)

vcfFile.close()
extHetFile.close()
cmsFile.close()	
logFile.close()
outputFile.close()
