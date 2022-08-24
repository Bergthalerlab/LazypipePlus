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
		df_uniq = df_qualmap[df_qualmap["contig"]==index]
		uniqmap_covered_bases = df_uniq["depth"].sum()
		uniqmap_goodcov = len(df_uniq[df_uniq["depth"]>=10])
		
		goodcov = len(data[data["depth"]>=10])
		goodcov_fraction = float(goodcov)/float(genomelength)
		uniq_goodcov_fraction = float(uniqmap_goodcov)/float(genomelength)
		normalized_coverage = float(covered_bases)/float(genomelength)
		uniq_normalized_coverage = float(uniqmap_covered_bases)/float(genomelength)
		df_genomecov = df_genomecov.append({"contig":index, "normalized_coverage":normalized_coverage, "uniq_normalized_coverage":uniq_normalized_coverage, "uniq_10x_fraction": uniq_goodcov_fraction, "10x_fraction":goodcov_fraction ,"virus_name":virname, "contig_name":index.split(":")[0]+"_"+virname}, ignore_index=True)

		t_data = data[data["depth"]>=5]
		if float(len(t_data))>=float(0.1*genomelength):
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
	
fig.update_layout(height=300*len(subplot_titles), width=1200, title_text = "Coverage plot")
fig.update_layout(showlegend=False)
#fig.update_yaxes(type="log")
fig.write_image(output_path+"coverage_plot"+virus+".png")
fig.write_html(output_path+"coverage_plot"+virus+".html")

if virus !="":
	quit()

fig2 = make_subplots(rows=2, cols=1, subplot_titles = ["Average depth", "Fraction of genome covered >=10x"])
df_genomecov = df_genomecov[df_genomecov["normalized_coverage"]>=1]
df_genomecov = df_genomecov.sort_values(by = "normalized_coverage", ascending = False)
allnames = list(df_genomecov["virus_name"].unique())
coldict = {}
for name, col in zip(allnames, px.colors.qualitative.Dark24*5):
	coldict.update({name:col})
df_genomecov["color"] = df_genomecov["virus_name"].apply(lambda x: coldict[x])
marker_colors = list(df_genomecov["color"])

for index, data in df_genomecov.groupby("virus_name"):
	fig2.append_trace(go.Bar(x=data['contig_name'], y=data['normalized_coverage'], opacity=0.5,showlegend = False, marker_color = coldict[index]), 1,1)
	fig2.append_trace(go.Bar(x=data['contig_name'], y=data['uniq_normalized_coverage'], name = index, marker_color = coldict[index]), 1,1)


#fig2.append_trace(go.Bar(x=df_genomecov['contig_name'], y=df_genomecov['normalized_coverage'], marker_color = marker_colors), 1,1)
df_genomecov = df_genomecov.sort_values(by = "10x_fraction", ascending = False)
for index, data in df_genomecov.groupby("virus_name"):
	fig2.append_trace(go.Bar(x=data['contig_name'], y=data['10x_fraction'], showlegend = False, opacity=0.5, name = index, marker_color = coldict[index]), 2,1)
	fig2.append_trace(go.Bar(x=data['contig_name'], y=data['uniq_10x_fraction'], showlegend = False, marker_color = coldict[index]), 2,1)


fig2.update_layout(height=2300, width=1800, barmode='overlay')
#fig2.update_yaxes(type="log")
fig2.update_xaxes(title_text="Virus")
fig2.update_yaxes(title_text="Depth of coverge", row = 1, col = 1)
fig2.update_yaxes(title_text="Genome fraction", row = 2, col = 1)


fig2.write_image(output_path+"average_coverage.png")
fig2.write_html(output_path+"average_coverage.html")

df_allqual.to_csv(output_path+"allqual_coverage.tsv", sep = "\t", index = False)
df_genomecov.to_csv(output_path+"average_coverage.tsv", sep = "\t", index = False)
