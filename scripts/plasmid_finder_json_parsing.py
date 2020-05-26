import sys
import json
import csv 
filename = sys.argv[1]

with open(filename) as data_json: 
    data=json.load(data_json)

# Accessing the strain and assembler details. 
filename = str(data['plasmidfinder']['user_input']['filename(s)'])
fn = filename.split('/')
strain = fn[-2]
assembler_raw = fn[-1].strip(']')
assembler= assembler_raw.strip("\'")

# Acessing the run details.
run_stats = data['plasmidfinder']['run_info']

# Accessing the results data. 
results = data['plasmidfinder']['results']["Enterobacteriaceae"]["enterobacteriaceae"]

print(results)

# Making header
header_list = []
header_list.append('Strain')
header_list.append('Assembley')
for run_stat in run_stats:
    header_list.append(run_stat.title())
header_list.append('Plasmid')
header_list.append('Identity')
header_list.append('HSP_length')
header_list.append('Template_length')
header_list.append('Position_in_ref')
header_list.append('Contig_name')
header_list.append('Position_in_contig')
header_list.append('Note')
header_list.append('Accession')
header_list.append('Coverage')
header_list.append('Hit_ID')
    
write_results = open('plasmid_finder.csv', "a")
csv_writer = csv.writer(write_results)

with open('plasmid_finder.csv') as f: 
    line = f.readline()
    if line == '':
        csv_writer.writerow(header_list)

# Making content list.
if results != "No hit found":
    for result in results:
        print(result)
        content_list = []
        content_list.append(strain)
        content_list.append(assembler)
        for value in run_stats.values():
            content_list.append(value) 
        hit = data['plasmidfinder']['results']["Enterobacteriaceae"]["enterobacteriaceae"][result]
        print(hit)
        for value in hit.values():
            content_list.append(value)
        csv_writer.writerow(content_list)
else:
    content_list = []
    content_list.append(strain)
    content_list.append(assembler)
    for value in run_stats.values():
            content_list.append(value)
    content_list.append("No plasmid found")
    csv_writer.writerow(content_list)

