import pandas as pd
import json
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
import sys

"""
arg 1: quality mappings tsv
arg 2: any-quality mappings tsv
arg 3: consider only genes (pox viruses and herpes viruses) OR only whole genomes (all other) [GENES/GENOMES]
arg 4: output path
"""


if sys.argv[3] in ["Genes", "GENES", "genes"]:
	output_path = sys.argv[4]+"/genes_"
elif sys.argv[3] in ["Genomes", "genomes", "GENOMES"]:
	output_path = sys.argv[4]+"/genomes_"
else:
	raise Exception("ERROR: Genes/genomes option not specified.")


virus = ""

with open("virdict_detail.json") as f:
	virdict = json.load(f)

def get_virname(code):
	return virdict[code.replace(";","")]["usual_name"]#+" clinical: "+virdict[code.replace(";","")]["clinical_typing"]

def get_genomelength(code):
	return virdict[code.replace(";","")]["genome_length"]

# real quality mappings
df_qualmap = pd.read_csv(sys.argv[1], sep="\t")
df_qualmap.columns = ["contig", "pos", "depth"]

# read allquality
df_allqual = pd.read_csv(sys.argv[2], sep="\t")
df_allqual.columns = ["contig", "pos", "depth"]

df_genomecov = pd.DataFrame()
if virus == "":
	df_allqual_master = pd.DataFrame()
	for index, data in df_allqual.groupby("contig"):
		if sys.argv[3] in ["Genes", "GENES", "genes"]:
			if "GENE" not in index:
				continue
		elif sys.argv[3] in ["Genomes", "genomes", "GENOMES"]:
			if "GENE" in index:
				continue
		else:
			raise Exception("ERROR: Genes/genomes option not specified.")
			
		genomelength = int(get_genomelength(index))
		virname = get_virname(data.iloc[0]["contig"])
		covered_bases = data["depth"].sum()
		normalized_coverage = float(covered_bases)/float(genomelength)
		df_genomecov = df_genomecov.append({"contig":index, "normalized_coverage":normalized_coverage, "virus_name":virname}, ignore_index=True)

		data = data[data["depth"]>=5]
		if float(len(data))>=float(0.1*genomelength):
			df_allqual_master = pd.concat([df_allqual_master, data])

df_allqual = df_allqual_master

df_allqual["virname"] = df_allqual.contig.apply(get_virname)

if virus != "":
	df_allqual = df_allqual[df_allqual["virname"].str.contains(virus, case = False)]

subplot_titles = []
for index, data in df_allqual.groupby("contig"):
	virname = get_virname(data.iloc[0]["contig"])+" "+index.split(":")[0]
	subplot_titles.append(virname)
	
fig = make_subplots(rows=len(list(df_allqual["contig"].unique())), cols=1, subplot_titles = subplot_titles)
cntr = 1
for index, data in df_allqual.groupby("contig"):
	#process qualmap
	virname = get_virname(data.iloc[0]["contig"])
	genomelength = int(get_genomelength(data.iloc[0]["contig"]))
	data = data.sort_values(by = ["pos"])
	data.index = data.pos
	data = data.reindex(list(range(0, genomelength)),fill_value=0)
	
	
	#process allqual
	qualmap = df_qualmap[df_qualmap["contig"]==index]
	qualmap = qualmap.sort_values(by = ["pos"])
	qualmap.index = qualmap.pos
	qualmap = qualmap.reindex(list(range(0, genomelength)),fill_value=0)
	
	
	fig.add_trace(go.Scatter(x=qualmap.index, y=qualmap["depth"], fill = "tozeroy", name = virname), row= cntr, col = 1)
	fig.add_trace(go.Scatter(x=data.index, y=data["depth"], mode = "lines", line=dict(color='black', width=1), name = virname), row= cntr, col = 1)
	fig.update_xaxes(title_text="Genomic position", row=cntr, col=1)
	fig.update_yaxes(title_text="Depth of coverge (log)", row=cntr, col=1)
	
	cntr += 1
	
fig.update_layout(height=200*len(subplot_titles), width=1200, title_text = "Coverage plot")
fig.update_layout(showlegend=False)
fig.update_yaxes(type="log")
fig.write_image(output_path+"coverage_plot.png")
fig.write_html(output_path+"coverage_plot.html")

df_genomecov = df_genomecov.sort_values(by = "normalized_coverage", ascending = False)
fig2 = px.bar(df_genomecov, x='virus_name', y='normalized_coverage', color = "contig")
fig2.update_yaxes(type="log")
fig2.write_image(output_path+"average_coverage.png")
fig2.write_html(output_path+"average_coverage.html")

df_allqual.to_csv(output_path+"allqual_coverage.tsv", sep = "\t", index = False)
df_genomecov.to_csv(output_path+"average_coverage.tsv", sep = "\t", index = False)
