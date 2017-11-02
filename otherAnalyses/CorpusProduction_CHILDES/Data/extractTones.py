import os,codecs
# extract tone type counts from CHILDES files

toneList = [1,2,3,4,5,6]


def readChildes(text):
	date = ""
	try:
		date = [x for x in text.split("\n") if x.startswith("@Date")][0][6:].strip()
	except:
		pass
	parts = [x for x in text.split("\n") if x.startswith("@Participants")][0][14:].strip().replace("\t",' ')
	
	partsX = [x.strip() for x in parts.split(",")]
	childID = "CHI"
	if parts.count("Target_Child")>0:
		childID = [x[:x.index(" ")] for x in partsX if x.count("Target_Child")>0][0]
	
	id = [x[4:].strip() for x in text.split("\n") if x.startswith("@ID:") and x.count("Target_Child")>0][0]
	#nan , eng|HKU|CHI|2;05.01|male|||Target_Child|||
#	ids = id.split("|")

	language,corpus,code,age,sex,group,SES,role,education,custom,x = id.split("|")
	
	speaker = "CHI"
	counts = dict(zip(toneList,[0 for x in toneList]))
	wordCount = 0
	for l in text.split("\n"):
		if l.startswith("*"):
			speaker = l[1:l.index(":")].strip()
		if l.startswith("%mor") and speaker != childID:
			for i in toneList:
				try:
					counts[i] += l.count(str(i))
				except:
					counts[i] = l.count(str(i))
			wordCount += l.count("|")
	countsString = "\t".join([str(counts[x]) for x in toneList])
	
	out = [date,parts,childID,language,corpus,code,age,sex,group,SES,role,education,wordCount,countsString]
	
	print "\t".join([str(x) for x in out])
	#return counts
		

header = "date	parts	childID	language	corpus	code	age	sex	group	SES	role	education	wordCount"
header += '\t'+'\t'.join(["T"+str(i) for i in toneList])

print header


o = open("/Users/pplsuser/Documents/MPI/ClimateAndLanguage/CHILDES/Data/Cantonese/HKU/20500.cha")
d = o.read()
o.close()
readChildes(d)


for root, dirs, files in os.walk("Cantonese", topdown=False):
	for name in files:
		# only larry king transcripts
		if name.count("meta")==0 and name.count("DS_Store")==0 and name.endswith(".cha"):
			#print name,root
			o = codecs.open(os.path.join(root, name),'r',encoding='utf-8')
			d = o.read()
			o.close()
			readChildes(d)