#!/usr/bin/env python3
# $1 = ensemb archive name eg. feb2021
# $2 = ourdir of csv results file
# $3..$n = chromosomes to retrieve

import sys
import re

import requests
from bs4 import BeautifulSoup as bs

import pandas as pd


ENSEMBL_ARCHIVE = sys.argv[1]
OUTDIR = sys.argv[2]
CHROMS = sys.argv[3:]


print(f"CHROMS {CHROMS}\n ENSEMBL_ARCHIVE {ENSEMBL_ARCHIVE}")

BASE_URL = f"http://{ENSEMBL_ARCHIVE}.archive.ensembl.org/Homo_sapiens/Location/Chromosome?chr="

def query_chrLen(CHR, BASE_URL):
	r = requests.get(BASE_URL + str(CHR))

	# connection status
	if r.status_code != 200:
		return (r.status_code for x in range(4))

	soup = bs(r.text, "html.parser")

	info = soup.head.title.contents[0]
	pat = r"Chromosome (?P<Chr>[0-9MTXY]{1,2}): (?P<Start>[0-9]+)-(?P<End>[0-9,]+) [^0-9]+ (?P<Version>[0-9]+)"
	re.compile(pat)
	m = re.search(pat, info)
	
	# Not get all the fields
	if not m:
		return (CHR, 0, 0, "notFound" )

	r_chr = "chr" + m.group("Chr")
	r_start = int(re.sub(",", "", m.group("Start")))
	r_end = int(re.sub(",", "", m.group("End")))
	r_ver =  m.group("Version")
		
	return (r_chr, r_start, r_end, r_ver)

data = []

for ch in CHROMS:
	print(f"querying chr{ch}\n")

	CHR, START, END, VERSION  = query_chrLen(ch, BASE_URL)

	# Something weird that could happen
	if "chr" + str(ch) != str(CHR):
		VERSION = f"ERROR! request chr{ch}, response chr{CHR}"
	data.append([CHR, START, END, VERSION])

to_csv = pd.DataFrame(data, columns=["Chr", "Start", "End", "Version"])
print(to_csv.to_string(index=False))
to_csv.to_csv(OUTDIR + "chrLens.csv", index= False)

