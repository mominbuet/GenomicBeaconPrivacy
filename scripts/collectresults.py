def byteify(input):
	"""
		http://stackoverflow.com/questions/956867/how-to-get-string-objects-instead-of-unicode-ones-from-json-in-python
		    converts unicode to byte strings. for use with json dictionaries
	"""
	if isinstance(input, dict):
		return {byteify(key):byteify(value) for key,value in input.iteritems()}
	elif isinstance(input, list):
		return [byteify(element) for element in input]
	elif isinstance(input, unicode):
		return input.encode('utf-8')
	else:
		return input

import json
import sys

outfile=sys.argv[1]

hittable={}
for infile in sys.argv[2:]:
	with open(infile,'r') as fp:
		rawd=json.load(fp)
		cleand=byteify(rawd)
		responsemap=[(x["beacon"]["name"],x["response"]) for x in cleand]
		chr=cleand[0]["query"]["chromosome"]
		pos=cleand[0]["query"]["position"]
		allele=cleand[0]["query"]["allele"]

		for (key,value) in responsemap:
			if value:
				#if key=="PGP":
				#	print key,value,chr,pos,allele
				if key not in hittable:
					hittable[key]=1
				else:
					hittable[key]=hittable[key]+1
			else:
				if key not in hittable:
					hittable[key]=0
				#if key=="1000 Genomes Project":
				#if key=="PGP":
				#if key=="Known VARiants":
				#if key=="AMPLab":
				#if key=="Wellcome Trust Sanger Institute":
				#	print key,value,chr,pos,allele

ofp=open(outfile,'w')

for (key,value) in hittable.iteritems():
	print key,value
	ofp.write(",".join([key,str(value)]))
	ofp.write("\n")

