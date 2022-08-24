import pandas as pd
import json
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
import sys

"""
arg 1: any-quality mappings tsv
arg 2: consider only genes (pox viruses and herpes viruses) OR only whole genomes (all other) [GENES/GENOMES]
arg 3: output path
"""

top_signals = 20 #pick top signals from each sample

if sys.argv[2] in ["Genes", "GENES", "genes"]:
	output_path = sys.argv[3]+"/genes_"
elif sys.argv[2] in ["Genomes", "genomes", "GENOMES"]:
	output_path = sys.argv[3]+"/genomes_"
else:
	raise Exception("ERROR: Genes/genomes option not specified.")

"""
optional argument - plots genome cov only for one virus name (all contigs)
"""
virus = ""

with open("virdict_detail.json") as f:
	virdict = json.load(f)

def get_virname(code):
	return virdict[code.replace(";","")]["usual_name"]#+" clinical: "+virdict[code.replace(";","")]["clinical_typing"]

def get_genomelength(code):
	return virdict[code.replace(";","")]["genome_length"]

def read_cov(cov_tsv, sample_name):

	# read allquality
	df_allqual = pd.read_csv(sys.argv[1]+cov_tsv, sep="\t")
	df_allqual.columns = ["contig", "pos", "depth"]

	df_genomecov = pd.DataFrame()
	df_allqual_master = pd.DataFrame()
	for index, data in df_allqual.groupby("contig"):
		if sys.argv[2] in ["Genes", "GENES", "genes"]:
			if "GENE" not in index:
				continue
		elif sys.argv[2] in ["Genomes", "genomes", "GENOMES"]:
			if "GENE" in index:
				continue
		else:
			raise Exception("ERROR: Genes/genomes option not specified.")
			
		genomelength = int(get_genomelength(index))
		virname = get_virname(data.iloc[0]["contig"])
		covered_bases = data["depth"].sum()
		goodcov = len(data[data["depth"]>=10])
		goodcov_fraction = float(goodcov)/float(genomelength)
		normalized_coverage = float(covered_bases)/float(genomelength)
		if normalized_coverage>=1:
			df_genomecov = df_genomecov.append({"sample": sample_name, "contig":index, "normalized_coverage":normalized_coverage, "10x_fraction":goodcov_fraction ,"virus_name":virname, "contig_name":index.split(":")[0]+"_"+virname}, ignore_index=True)

	return df_genomecov

df = pd.DataFrame()
for i in [[str(i)+"R_allmap.tsv", "R"+str(i)] for i in range(1,17)]:
	try:
		df = pd.concat([df, read_cov(i[0], i[1])])
	except:
		pass

df = df[["sample", "contig_name", "10x_fraction"]]
df = df.sort_values(by = ["10x_fraction"], ascending = False)
df = df.pivot(index="contig_name", columns="sample", values="10x_fraction")
df = df.fillna(0)

# select top n signals from each sample
selected_contigs = []
for column in df.columns:
	df1 = df.sort_values(by = [column], ascending = False)
	df1 = df1[:top_signals]
	selected_contigs.append(list(df1.index))
selected_contigs = list(set(list(selected_contigs[0])))
df = df[df.index.isin(selected_contigs)]

# sort values by first column
df = df.sort_values(by = [df.columns[0]], ascending = False)
df_vals = df.values.tolist()
fig = px.imshow(df_vals, x = list(df.columns), y = list(df.index), labels=dict(x="Sample", y="Virus", color = "Fraction of genome covered >=10x"))
fig.update_layout(height=len(df)*25, width=len(df.columns)*60)
fig.update_xaxes(side="top")
fig.write_html("respiratory_panel_coverage.html")


