import csv
import sys

filename = sys.argv[1]
i=filename.split("/")
assembler= i[-3]
output_filename = assembler +"_mlplasmids.tab"
strain= i[-4]
d =open(filename,"r").readlines()

write_results = open(output_filename, "a")
csv_writer = csv.writer(write_results)

Header_list=["Prob_chromosome","Prob_plasmid","Prediction","Contig_name","Contig_length"]

with open(output_filename) as f: 
    line = f.readline()
    if line == '':
        csv_writer.writerow(Header_list)
csv_writer.writerow([strain])
csv_writer.writerow(d[1:])
