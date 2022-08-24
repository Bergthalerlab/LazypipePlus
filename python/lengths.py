import sys
import pandas
import json
import plotly.express as px


with open("virdict_detail.json") as f:
	virdict = json.load(f)

lengths = []
for key, val in virdict.items():
	if "GENE" in key:
		continue
	else:
		if val["genome_length"] < 50000:
			lengths.append(val["genome_length"])

fig = px.violin(y=lengths, points='all')
#fig.update_yaxes(type="log")
fig.show()
