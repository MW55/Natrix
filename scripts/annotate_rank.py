import urllib.request
from bs4 import BeautifulSoup
import pandas as pd

file = snakemake.input[0]

data = pd.read_csv(file)

ids = data.tail(1)

result = {}

LINK_TEMPLATE = "https://www.datapunk.net/opus23/utopia/display.pl?"

for tax in ids:
	if tax == "Unnamed: 0":
	    continue
	_id = int(ids[tax])
	link = f"{LINK_TEMPLATE}{_id}"
	result[tax] = []

	with urllib.request.urlopen(link) as response:
		html_doc = response.read().strip()
		soup = BeautifulSoup(html_doc, 'html.parser')
		TAGS = soup.body.find(text="TAGS")
		tags_anchor = next(soup.body.find(text="TAGS").parents)
		tags_table = tags_anchor.find_next("table")
		tag_entries = tags_table.find_all("td")
		for entry in tag_entries:
			styles = entry.get('style').split(";")
			if("font-weight: bold") in styles:
				result[tax].append(entry.get_text().strip())

for res in result.keys():
    result[res] = [";".join(result[res])]

result["Unnamed: 0"] = "tags"
tags = pd.DataFrame.from_dict(result)
data = pd.concat([data,tags])
data = data.set_index("Unnamed: 0")
#data = data.rename(columns={"Unnamed: 0": ""})

data.to_csv(snakemake.output[0])
