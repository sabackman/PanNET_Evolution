#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

SNP concordance heat maps found in supplementary figures.

"""

import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd

def load_data(folder):
	files = os.listdir(folder)
	samples = {}
	sites = set()
	for f in files:
		samples[f.split("_")[0]] = {}
	for f in files:
		sample = f.split("_")[0]
		with open(os.path.join(folder, f), "r") as infile:
			for line in infile:
				if line.strip() == "":
					continue
				words = line.strip().split("\t")
				site = words[0]+"-"+words[1]
				if not site in sites:
					for s in samples.keys():
						samples[s][site] = ""
					sites.add(site)
				
				reads = words[4].upper()
				counts = {}
				for allele in ["A", "T", "C", "G"]:
					counts[allele] = reads.count(allele)
				samples[sample][site] = max(counts, key=counts.get)

	return samples, sites
				
def process_data(samples, sites, chromosomes=["chr1"]):
	rownames = samples.keys()
	colnames = samples.keys()
	
	columns = {}
	
	for sample1 in colnames:
		column = []
		for sample2 in rownames:
			concordant=0
			discordant=0
			missing=0
			for site in sites:
				if site.split("-")[0] not in chromosomes:
					continue
				if samples[sample1][site]=="" or samples[sample2][site]=="":
					missing+=1
					continue
				if samples[sample1][site]==samples[sample2][site]:
					concordant += 1
				else:
					discordant += 1
			column.append(float(concordant)/(concordant+discordant))
			print(concordant, discordant, missing)
		columns[sample1] = column
	return pd.DataFrame(columns, index=rownames)

samples, sites = load_data("PanNET3_pileups_small")
results=process_data(samples, sites, chromosomes=["chr1", "chr2", "chr3", "chr6", 
												  "chr8", "chr9", "chr10", "chr11"
												  "chr15", "chr16", "chr21", "chr22"])

results = results.rename(mapper={'7382':"M3b", 
								 '7385':"M4a", 
								 '5621':"P1a", 
								 '7370':"M1", 
								 '7372':"M2", 
								 '7386':"M4b", 
								 '7379':"M3a", 
								 '765-3':"P1b",
       '765-5':"P1c"}, axis=0)

results = results.rename(mapper={'7382':"M3b", 
								 '7385':"M4a", 
								 '5621':"P1a", 
								 '7370':"M1", 
								 '7372':"M2", 
								 '7386':"M4b", 
								 '7379':"M3a", 
								 '765-3':"P1b",
       '765-5':"P1c"}, axis=1)

samples, sites = load_data("PanNET4_pileups_small")
results2=process_data(samples, sites, chromosomes=["chr1", "chr2", "chr3", "chr6", 
												  "chr8", "chr10", "chr11"
												  "chr15", "chr16", "chr18", "chr21", "chr22"])

results2=results2.rename(mapper={'320':"P1a",
								 '2586':"M1b",
								 '2580':"P1b",
								 '320-lever':"M1c",
								 '320-10':"P1c", 
								 '2600':"M2b",
								 '2597':"M2a",
								 '6680':"M3",
       '2584':"M1a"}, axis=0)
results2=results2.rename(mapper={'320':"P1a",
								 '2586':"M1b",
								 '2580':"P1b",
								 '320-lever':"M1c",
								 '320-10':"P1c", 
								 '2600':"M2b",
								 '2597':"M2a",
								 '6680':"M3",
       '2584':"M1a"}, axis=1)

sns.heatmap(results, vmin=0.9, annot=True)
plt.savefig("PanNET3_concordance.png", dpi=300)
plt.figure()
sns.heatmap(results2, vmin=0.9, annot=True)
plt.savefig("PanNET4_concordance.png", dpi=300)