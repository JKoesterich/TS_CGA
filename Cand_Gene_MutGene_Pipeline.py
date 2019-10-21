#Candidate Gene analysis Pipeline
import subprocess 
import argparse
import sys
import os
#import statistics
import numpy
#Add arguments to the script to allow it to be modifiable
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--help", "-h", help="print", action="help")
parser.add_argument("--files", "-f", help="This takes in the files that are going to be iterated over", nargs='*', required= True)
parser.add_argument("--anno", "-a", help="This takes in the annotation file and will pull Sift, Polyphen-2, Exac, gnomAD, and clinvar data", nargs=1)
parser.add_argument("--gtex", "-g", help="This takes in the median tpm gtex file and will gather max and median values for the different tissues for a gene", nargs = 1)
parser.add_argument("--disease", "-d", nargs ="*", help="This takes in a disease association file and will output the disease associations via keyword search, will print all disease if second file not provided. Second file is 1 column with headings of general for lowercase matching, specific for exact matching, and exclude for lowercase matches to be excluded. Will treat all cases below and between headings as part of the heading")
#parser.add_argument("--segpat", "-s", help= "This takes in the segregation pattern for the families and will output the segregation per mutation", nargs=1)
parser.add_argument("--segpat","-s", help= "This takes in the .reflist files for the families to make the segregation patterns. the files names of the reflist must have family name in it so that the family files could be connected to the segregation pattern", nargs="*")
parser.add_argument("--genes", "-e", help= "This takes in sets of 2 information, the first being a file of newline separated genes and the second item is where they are from, Ex: <file> GWAS will put GWAS_found in file name and 'found in GWAS' in mutation line", nargs="*")
parser.add_argument("--output", "-o", help= "This is the output path where the files will print to", required= True)
parser.add_argument("--genetol", "-t", help= "This takes in a file that has the pLI, mis_z, and protein cluster information", nargs=1)
argus = parser.parse_args()
#Parser assigns the long name as the object name, or the 1 character name if a longer one is not provided 

def segpattern(seg, famm, unge, linete): #Make definition for segregation pattern, taking in the mutation segregation, family segregation, unknown genotypes and the full mutation line
	unpe = []
	incl = []
	excl = []
	#fam = famm.strip('"').split(",")
	fam = famm
	for te in famm: # iterate over the family segregation pattern
		if (te.find("?") > -1): # "?" is in the segregation pattern file to represent unknown phenotype
			tem = te[:-1] # Remove the ? and add the individual to the unknown phenotype list
			unpe.append(tem)
			fam.remove(te)
			fam.append(tem)
	for found in seg: # for every individual that has the mutation
		if(found not in fam and found not in unge): # If the mutation is not in the segregation pattern and not an unknown genotype
			incl.append(found) # Add the individual to the included list
		elif(found in fam and found in unpe): # If the mutation is in the segregation pattern and is one of the unknown phenotypes
			temm = found + "?P"
			incl.append(temm) # Add ?P to the individual and add it to the included list to signify that phenotype is unknown
		elif(found in fam and found in unge): # If mutation is found in the segregation pattern and is of unknown genotype
			temm = found + "?G"
			excl.append(temm) # Add ?G and add to the excluded because the genotype is unknown and the individual might be reference allele
		elif(found not in fam and found in unge): # If the mutation is not in the segregation and has unknown genotype
			temm = found + "?G"
			incl.append(temm) # Add ?G and add to the included because the genotype might be alternate for unaffected individual
	for pat in fam: # For every individual in the family segregation pattern
		if(pat not in seg and pat not in unpe): # If the individual was not in the mutation segregation and not an unknown phenotype
			excl.append(pat) # Add the individual to excluded
		elif(pat not in seg and pat in unpe): # If the individual was not in the mutation segregation and was unknown phenotype
			temm = pat + "?P"
			excl.append(temm) # Add ?P and add it to the excluded

	exclstr = ""
	inclstr = ""
	sgpt = ""
	if(len(excl) == 0 and len(incl) == 0): # If there was no included or excluded individuals the segregation pattern is perfect
		sgpt = "Perfect"
	if(len(excl) > 0): # If there were any exlcuded indivudals 
		exclstr = "Excludes: " + ",".join(excl)
		sgpt = sgpt + " " + exclstr
	if(len(incl) > 0): # If there were any included individuals
		inclstr = "Includes: " + ",".join(incl)
		sgpt = sgpt + " " + inclstr

	#print(seg, fam)
	#print(sgpt)
	linete.append(sgpt) #Take the mutation line, add the segregation pattern to the end and return the new line
	return(linete)

def findannot(varty, chromo, varposit, linite): # Definition function for linding annotation for mutations.
# Function takes in variant type, chromosome as it would appear in annotation file, variant position, and the whole mutation line
	varpk = chromo + varposit
	varpn1 = int(varposit) -1
	varn1k = chromo + str(varpn1)
	varpch = -1
	varn1ch = -1
	if(varpk in darction.keys()):
		varpch = 1
	if(varn1k in darction.keys()):
		varn1ch = 1

	if(varty.lower() == "del" and varn1ch == 1):
		lifo = darction[varn1k]
	elif(varty.lower() == "mnp" and (varpch == 1 or varn1ch == 1)):
		if(varn1ch == 1):
			lifo = darction[varn1k]
		if(varpch == 1):
			lifo = darction[varpk]
	elif(varpch == 1):
		lifo = darction[varpk]

	else:
		for z in range(0,6):
			linite.append("Mutation not found")
			return(linite)
	
# Assuming the correct line has been found the function continues here
	
	lifter = lifo.split(";")
	siftli = []
	pol1li = []
	pol2li = []
	poltogli = []
	exacli = []
	gnomadli = []
	clinvarli = []
	clinvartot = []
	interproli = []
	for linef in lifter: # Search for the sift information
		if (linef.find("SIFT_pred=") > -1):
			pot = linef.find("=")
			sifti = linef[pot+1:]
			if(sifti == "D"):
				sift = "Deleterious"
			elif(sifti == "T"):
				sift = "Tolerated"
			elif(sifti == "."):
				sift = "."
			else:
				print("ERROR IN SIFT")
			siftli.append(sift)

		elif(linef.find("Polyphen2_HDIV_pred=") > -1): # Search for polyphen information
			pot = linef.find("=")
			poly1i = linef[pot+1:]
			if(poly1i == "D"):
				poly1 = "Damaging"
				pol1 = 3
			elif(poly1i == "P"):
				poly1 = "Possibly Damaging"
				pol1 = 2
			elif(poly1i == "B"):
				poly1 = "Benign"
				pol1 = 1		
			elif(poly1i == "."):
				poly1 = "."
				pol1 = 0
			else: 
				print("ERROR IN POLYPHEN_1")
			pol1li.append([poly1,pol1])

		elif(linef.find("Polyphen2_HVAR_pred=") > -1): # Search for polyphen information
			pot = linef.find("=")
			poly2i = linef[pot+1:]
			if(poly2i == "D"):
				poly2 = "Damaging"
				pol2 = 3
			elif(poly2i == "P"):
				poly2 = "Possibly Damaging"
				pol2 = 2
			elif(poly2i == "B"):
				poly2 = "Benign"
				pol2 = 1
			elif(poly2i == "."):
				poly2 = "."
				pol2 = 0
			else:
				print("ERROR IN POLYPHEN_2")
			pol2li.append([poly2,pol2])

		elif(linef.find("ExAC_nonpsych_ALL=") > -1):
			pot = linef.find("=")
			subse = linef[pot+1:]
			exacli.append(subse)

		elif(linef.find("gnomAD_genome_ALL=") > -1):
			pot = linef.find("=")
			subse = linef[pot+1:]
			gnomadli.append(subse)

		elif(linef.find("CLN") > -1 or linef.find("CLIN") > -1): # Add all clinvar information
			clinvarli.append(linef)

		elif(linef.find("Interpro_domain") > -1):
			le = linef.find("=")
			subse = linef[le+1:]
			interproli.append(subse)

		elif(linef.find("ALLELE_END") > -1): # Allele end is at the end of each allele, if it is a mnp there is annotation for each allele
			clinvarstr = ",".join(clinvarli) # join all the clinvar information for this allele and add it to a holder list
			clinvartot.append(clinvarstr)
			clinvarli = []
	#print(pol1li)
	#print(pol2li)
	#print(siftli)
# After all the information for all alleles are gathered
	for xen in range(0, len(pol1li)): # Go through the list of polyphen2 information
		pol1n = pol1li[xen][1]
		pol1s = pol1li[xen][0]
		pol2n = pol2li[xen][1]
		pol2s = pol2li[xen][0]
		if(pol1n > pol2n): # Take the number score assigned earlier and compare to decide how the string is written
			poltog = pol1s + "/" + pol2s
		elif(pol2n > pol1n):
			poltog = pol2s + "/" + pol1s
		elif(pol1n == pol2n):
			poltog = pol1s
		poltogli.append(poltog) # Add the new string to a new list
# Join all the lists into strings separating allele infromation by ;
	siftfin = ";".join(siftli)
	poltogfin = ";".join(poltogli)
	exacfin = ";".join(exacli)
	gnomadfin = ";".join(gnomadli)
	clinvarfin = ";".join(clinvartot)
	interprofin = ";".join(interproli)
	linite.append(siftfin) # Add all the new information to the whole mutation line and return the new line
	linite.append(poltogfin)
	linite.append(exacfin)
	linite.append(gnomadfin)
	linite.append(clinvarfin)
	linite.append(interprofin)
	
	return(linite)

if(argus.anno != None):
	darction = {}
	darta = open(argus.anno[0], 'r')
	for larn in darta:
		if(larn.find("#") == 0):
			continue
		larnn = larn.strip()
		larnt = larnn.split("\t")
		ke = larnt[0] + larnt[1]
		darction[ke] = larnt[7]
	darta.close()

if(argus.genes != None): # If the --genes/-g parameter is being used
	gensetli = []
	if(len(argus.genes) % 2 != 0): # check to see if the order is: file heading_name
		sys.exit("Error: --genes requires sets of 2, first file then heading name. Both must be supplied")
	else:
		genewq = argus.genes
		for po in range(0,len(genewq),2):
			holdgenli = []

			if(os.path.isfile(genewq[po]) == False):
				sys.exit("Error: file not found. First item of each set for --genes must be a file")
			else: # if there are no errors then open the file and assuming the genes are separated by newlines add them to a list
				dq = open(genewq[po], 'r')
				for dpq in dq:
					holdgenli.append(dpq.strip())
				dq.close()
			oulis = [genewq[po+1], holdgenli] # Add to the list the heading name and the gene list to allow for multiple sets of --genes to be done at once
			gensetli.append(oulis)
	#print(gensetli)


if(argus.segpat != None): # If --segpat is being used and is supplied a file the is family name \t comma separated family segregation
	diction = {}
	#data = open(argus.segpat[0], 'r') # Open the supplied file and make a dictionary holding the family name as key and comma separated segregation pattern as the value
	#for lin in data:
	#	linn = lin.strip()
	#	lint = linn.split("\t")
	#	diction[lint[0]] = lint[1]
	#data.close()
	for pat in argus.segpat:
		see = os.path.basename(pat)
		sq = see.find(".")
		seo = see[:sq]
		patter = open(pat, 'r')
		holind = []
		for lat in patter:
			latn = lat.strip()
			latt = latn.split("\t")
			if(latt[5] == "Aff"):
				holind.append(latt[0])
		diction[seo] = holind
		patter.close()

if(argus.disease != None): # If --disease is being used and the tab separated value file is supplied. File must have fields: Gene, disease, and mental health
	#keywrd = ["autism", "epilepsy", "intellectual", "language", "adhd", "ocd","obsessive","attention","tourette","huntington","parkinson","alzheimer", "brain"]
	keywrd = []
	speckey = []
	exkey = []
	# The list above is keywords that I chose that would be part of diseases that are relavent for Tourettes. For other diseases this list and mental health must be edited
	if(len(argus.disease) == 2):
		dapper = open(argus.disease[1], 'r')
		typer = 0
		for lio in dapper:
			if(lio.lower().find("general") > -1):
				typer = 1
				continue
			elif(lio.lower().find("specific") > -1):
				typer = 2
				continue
			elif(lio.lower().find("exclude") > -1):
				typer = 3
				continue
			if(typer == 1):
				keywrd.append(lio.lower().strip())
				continue
			elif(typer == 2):
				speckey.append(lio.strip())
				continue
			elif(typer == 3):
				exkey.append(lio.lower().strip())
				continue
		dapper.close()
	daction = {}
	dats = open(argus.disease[0], 'r')

	for lan in dats:
		lann = lan.strip()
		lani = lann.split("|")
		if(lani[0].find("Gene") > -1 and lani[1].find("Disease") > -1):
			continue
		else:
			gens = lani[0].strip()
			if (gens not in daction.keys()): # Add the gene to the dictionary if not already there
				daction[gens] = []
			lanturn = lani[1].split(";")
			for lite in lanturn:
				dises = lite.strip()
				if(len(keywrd) != 0): # Look to see if there are any general key words 	
					for chk in keywrd: # Iterate over the keywords
						if(dises.lower().find(chk) > -1): # If the key word is found in disease name
							exer = 0
							for bis in exkey: # Check for excluded words from the list
								if(dises.lower().find(bis) > -1):
									exer = 1 # If found then the disease name has a word that doesn't want to be included so it is not added to the list of diseases
									break
							if(exer == 0): # If excluded word is not found then add the disease name to the list
								daction[gens].append(dises)
				if(len(speckey) != 0): # Looking for specific key words
					for chel in speckey:
						if(dises.find(chel) > -1):
							exers = 0
							for bise in exkey:
								if(dises.lower().find(bise) > -1):
									exers = 1
									break
							if(exers == 0):
								daction[gens].append(dises)
				if(len(keywrd) == 0 and len(speckey) == 0 and len(exkey) != 0): # If only excluded keywords are provided
					exerses = 0
					for bises in exkey:
						if(dises.lower().find(bise) > -1):
							exerses = 1
							break
					if(exerses == 0):
						daction[gens].append(dises) # If the disease name does not have an exlcuded word in it then it will print the disease
				if(len(keywrd) == 0 and len(speckey) == 0 and len(exkey) == 0): # If a keyword file was not given then these 3 lists would be blank and all disease will be added to the disease line
					daction[gens].append(dises)
				
	dats.close()

if(argus.genetol != None):
	dtoltion = {}
	natol = []
	genpos = -1
	plipos = -1
	mis_zpos = -1
	propos = -1
	plich = False
	mis_zch = False
	procluch = False
	datol = open(argus.genetol[0], 'r')
	for lintol in datol:
		lintos = lintol.strip().split("\t")
		#print(lintos)
		if(lintol.lower().find("symbol") > -1 and lintol.find("pLI")):
			for na in range(0, len(lintos)):
				if(lintos[na].lower() == "symbol"):
					genpos = na
				elif(lintos[na] == "pLI"):
					plipos = na
				elif(lintos[na].lower() == "mis_z"):
					mis_zpos = na
				elif(lintos[na].lower().find("protein") > -1):
					propos = na
			if(plipos == -1):
				pass
			else:
				plich = True
				natol.append("pLI")
			if(mis_zpos == -1):
				pass
			else:
				mis_zch = True
				natol.append("mis_z")
			if(propos == -1):
				pass
			else:
				procluch = True
				natol.append("Protein cluster")
			#tollihold = tolli
		else:
			tolli = []
			gen = lintos[genpos]
			if(plich == True): 
				tolli.append(lintos[plipos])
			if(mis_zch == True):
				tolli.append(lintos[mis_zpos])
			if(procluch == True):
				tolli[x].append(lintos[propos])
			dtoltion[gen] = tolli		
	datol.close()

if(argus.gtex != None): # If --gtex and a median file is supplied 
	duction = {}
	dut = open(argus.gtex[0], 'r')
	brnpos = []
	nonbrnpos = []
	genposi = -1
	strter = 0
	for lun in dut:
		#print(lun)
		alltpm = []
		brntpm = []
		nonbrntpm = []
		#maxall = []
		#maxbrn = []
		#nonbrnmax = []
		lunn = lun.strip()
		lunt = lunn.split("\t")
		if(lun.lower().find("description") > -1): # Find the line that says description as that is the header line
			strter = 1
			for h in range(0,len(lunt)):
				if(lunt[h].lower().find("description") > -1): # Assuming that description is the column that holds the gene names
					genposi = h
				elif(lunt[h].lower().find("brain") > -1): # Find the columns that are brain tissues and add it to the list
					brnpos.append(h)
			headip = lunt[genposi+1:] # Store the tissue names for later use 
		elif(strter == 1): # For all lines after the heading line 
			for hp in range(genposi+1,len(lunt)):
				alltpm.append(float(lunt[hp])) # Add the tpm for each tissue to the all tpm list
				if(hp in brnpos):
					brntpm.append(float(lunt[hp])) # Add the tpm for the brain tissues to its list
				else:
					nonbrntpm.append(float(lunt[hp])) # Add the tpm for the non-brain tissues to its list
			#print(brntpm, nonbrntpm, alltpm)
			gol = numpy.argsort(nonbrntpm)[len(nonbrntpm)//2] # Use numpy.argsort to sort the list based on number and the len()//2 takes the middle index position which would be the median tpm
			mol = numpy.argsort(brntpm)[len(brntpm)//2]
			dol = numpy.argsort(alltpm)[len(alltpm)//2]
			maxbrn = max(brntpm) # Find the max tpm 
			maxbind = alltpm.index(maxbrn) # Find the index of the max tpm
			maxbrnstr = headip[maxbind] # Use the index to hold the tissue of the max tpm
			maxnonbrn = max(nonbrntpm) # Repeat the above 3 steps for the nonbrain, all tissues, and for the median tpms
			maxnbind = alltpm.index(maxnonbrn)
			maxnonbrnstr = headip[maxnbind]
			maxall = max(alltpm)
			maxaind = alltpm.index(maxall)
			maxallstr = headip[maxaind]
			medbrn = brntpm[mol]
			medbrnin = alltpm.index(medbrn)
			medbrnstr = headip[medbrnin]
			mednonbrn = nonbrntpm[gol]
			mednonbrnin = alltpm.index(mednonbrn)
			mednonbrnstr = headip[mednonbrnin]
			medall = alltpm[dol]
			medallin = alltpm.index(medall)
			medallstr = headip[medallin]
			if(maxbrn == 0.0): # If the max brain tissue is 0.0 then assumed that all brain tissues are 0.0. Checked for nonbrain and all and changed the output string
				maxbrnstr = "All Brain 0.0 tpm"
			#if(medbrn == 0.0):
				medbrnstr = "All Brain 0.0 tpm"
			if(maxnonbrn == 0.0):
				maxnonbrnstr = "All Non Brain 0.0 tpm"
			#if(mednonbrn == 0.0):
				mednonbrnstr = "All Non Brain 0.0 tpm"
			if(maxall == 0.0):
				maxallstr = "All tissues 0.0 tpm"
			#if(medall == 0.0):
				medallstr = "All tissues 0.0 tpm"
# Create the output strings
			maxbout = "Max Brain: " + maxbrnstr + "|" + str(maxbrn)
			maxnbout = "Max Non Brain: " + maxnonbrnstr + "|" + str(maxnonbrn)
			maxaout = "Max All Tissues: " + maxallstr + "|" + str(maxall)
			medbout = "Median Brain: " + medbrnstr + "|" + str(medbrn)
			mednbout = "Median Non Brain: " + mednonbrnstr + "|" + str(mednonbrn)
			medallout = "Median All Tissues: " + medallstr + "|" +str(medall)
			fullout = "{};{};{};{};{};{}".format(maxbout,maxnbout,maxaout,medbout,mednbout,medallout)
			duction[lunt[genposi]] = fullout


	dut.close()

for x in argus.files: # For each file that is being used 
	fil = open(x, 'r')
	xer = os.path.basename(x)
	xam = xer.find(".")
	outnam = xer[:xam]
	if(argus.segpat != None): # For each of the optional arguments add the appropriate file additions
		outnam = outnam + "_segpat"
	if(argus.anno != None):
		outnam = outnam + "_annotated"
	if(argus.genetol != None):
		outnam = outnam + "_gene_tolerance"
	if(argus.disease != None):
		outnam = outnam + "_disease_assoc"
	if(argus.gtex != None):
		outnam = outnam + "_tissue_expression"
	if(argus.genes != None):
		for chep in gensetli:
			outnam = outnam + "_" + chep[0]
	if(os.path.exists(argus.output) == False): # Look for problems with the provided output and give errors if they appear
		sys.exit("Error: output must be a valid folder")
	elif(os.path.isfile(argus.output) == True):
		sys.exit("Error: output must be a folder to print to not a file")
	if(argus.output[-1] == "/"):
		outnammid = argus.output + outnam
	else:
		outnammid = argus.output + "/" + outnam
	outnamfin = outnammid + ".tsv"
	#print(outnamfin)
	outfile = open(outnamfin, 'w')
	#####
	gn = -1
	cr = -1
	vrty = -1
	vrps = -1
	gnty = -1
	holgen = ""
	holchrom = ""
	holvarpos = ""
	for line in fil: # For each line in the file
		#linen = line.strip()
		linen = line.rstrip()
		#linet = linen.split("\t")
		linet = linen.split("\t")

		if(line.find("CMD") > -1):
			pass
		elif(line.lower().find("gene") > -1 and line.lower().find("varposition") > -1): # Look for the heading line
			for q in range(0,len(linet)): # Save the index of each of the wanted headings
				if(linet[q].lower() == "gene"): 
					gn = q
				elif(linet[q].lower() == "chr"):
					cr = q
				elif(linet[q].lower() == "vartype"):
					vrty = q
				elif(linet[q].lower() == "varposition"):
					vrps = q
				elif(linet[q].lower() == "genotype"):
					gnty = q

			if(argus.segpat != None): # If --segpat was used, Add the segregation pattern to the heading line
				for fa in diction.keys():
					if (x.lower().find(fa.lower()) > -1):
						famseg = diction[fa]
						outfams = "Segregation information: " + ",".join(famseg)
				linet.append(outfams)

			if(argus.anno != None): # If doing annotation information add the following to the heading line
				linet.append("Sift score")
				linet.append("Polyphen-2 score")
				linet.append("Exac AF")
				linet.append("gnomAD AF")
				linet.append("Clinvar information")
				linet.append("Interpro_domain information")

			if(argus.genetol != None):
				for naf in natol:
					linet.append(naf)

			if(argus.disease != None): # If doing disease association add the following columns to the heading line
				linet.append("Disease association:Keyword")

			if(argus.gtex != None):
				linet.append("GTEx tpm information")

			if(argus.genes != None):
				for tepper in gensetli:
					linet.append(tepper[0])

			midlinet = "\t".join(linet) # Join the list of header line after adding all the columns and write out the heading line
			#print(midlinet)
			outfile.write(midlinet)
			outfile.write("\n")
			#####
		else: # For all the mutation lines
			if(len(linet) < 5):
				continue
			gen = linet[gn]
			chrom = linet[cr] 
			vartype = linet[vrty]
			varpos = linet[vrps]
			geno = linet[gnty] # If there is no segregation pattern then gnty would still be -1 and not report an index error
			if(gen != "" and chrom != ""):
				holgen = gen
				holchrom = chrom
			else:
				gen = holgen
				chrom = holchrom
			if(argus.segpat != None): # For finding segregation pattern
				gentype = []
				holds = []
				unknogen = [] 
				genop = geno.strip('"').split(";") # Each mutation is separated by semi colons so separate on the semi colon
				for genops in genop:
					genoper = genops.split("|") # The individuals, bp change, and aa change are separated by |
					check = any(che.isdigit() for che in genoper[0]) # For each character in the first element check to see if there is at least 1 number in the string
					if(check): # If the first element is the indviduals, the index will be held
						cc = 0
					else: # If it isn't then the first element is N or B so the second element will hold the individuals and that index will be held
						cc = 1
					lev = genoper[cc]
					lever = lev.split(",") # Separate the individuals by comma 
					for comma in lever:
						if(comma.find("-") > -1): # If the indviduals are in a range separate the first and second number and use range to pull the individual numbers
							dash = comma.find("-")
							rnst = int(comma[:dash])
							rnen = int(comma[dash+1:])
							rang = list(range(rnst,rnen+1))
							for lippp in rang:
								gentype.append(str(lippp))
								holds.append(str(lippp))
						else: # If there is no range then add the individual to the list
							gentype.append(comma)
							holds.append(comma)
					if(genoper[cc+1].find("^") > -1): # If the bp change is unknown "^" then the numbers preceding it will be held in the unknown genotype list as well as the regular list
						for lippers in holds:
							unknogen.append(lippers)
						holds = []
					else:
						holds = []
				templin = segpattern(gentype,famseg, unknogen, linet) # Send the information to the segregation pattern function
				linet = templin
				#print(linet)
			if(argus.anno != None):
				tikki = findannot(vartype, chrom, varpos, linet) # Send the information to the annotation function
				linet = tikki
				#print(tikki)
				
			if(argus.genetol != None):
				if(gen not in dtoltion.keys()):
					tekey = dtoltion.keys()[0]
					num = len(dtoltion[tekey])
					for fn in range(0,num):
						linet.append("Gene not found")
				else:
					teval = dtoltion[gen]
					for fnfon in teval:
						if(fnfon == ""):
							linet.append(".")
						else:
							linet.append(fnfon)
				#print(linet)

			if(argus.disease != None): # If the disease association is being done
				if (gen not in daction.keys()):
					keystr = "Gene not found"
				else:
					keyli = daction[gen] # Pull the information from the keyword and mental health list and combine them into 2 strings
					if(len(keyli) == 0):
						keystr = "."
					else:
						keystr = ";".join(keyli)

				linet.append(keystr)
				

			if(argus.gtex != None):
				if (gen not in duction.keys()): # If the gene is not in the gtex file print the error
					harp = "Gene not found"
				else: # Otherwise get the gene information and add it to the mutation list
					harp = duction[gen] 
				linet.append(harp)


			if(argus.genes != None):
				for tems in gensetli: # For each of the gene list files if the gene is part of the list add that it was found in the heading listed, otherwise add a "."
					if gen in tems[1]:
						stry = "Found"
					else:
						stry = "."
					linet.append(stry)

			midlinet = "\t".join(linet) # After all of the different arguments have been done, combine the list into a single list and write out that line
			#print(midlinet)
			outfile.write(midlinet)
			outfile.write("\n")
			#####
	fil.close() # After all the lines in a file have been done, close the file and close the output written file
	outfile.close()
	#####
	finstring = x + " file has been completed {} of {}".format(argus.files.index(x)+1, len(argus.files))
	print(finstring)
