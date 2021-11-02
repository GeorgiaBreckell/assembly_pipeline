# Script to take all my metrics and output a single CSV file which will emiliate the need for the google sheet. 
# First metric is number of short ORFs. 
import csv
import sys
import numpy as np 
import pandas as pd 

ORF_filename = sys.argv[1]
genome_stats = sys.argv[2]
#mlplasmids= sys.argv[3]
plasmidfinder= sys.argv[3]
concordant_reads= sys.argv[4]
non_mapping_reads= sys.argv[5]
socru= sys.argv[6]

##Setting up output file 
output_filename= "Assembly_metrics.tab"
#write_results = open(output_filename, "a")
#csv_writer = csv.writer(write_results)


# Adding ORF data. 
ORF_split = ORF_filename.split("/")
strain=ORF_split[-4]
assembler=ORF_split[-3]
print("Strain is " + strain)
print("Assembler is " + assembler)


ORF=pd.read_table(ORF_filename, header=None)
col1=ORF.iloc[:,0]
col2=ORF.iloc[:,1]
num_rows=len(ORF.index)
short_ORF= ORF.iloc[:,0]/ORF.iloc[:,1]
short_count=short_ORF[short_ORF<0.90].count()
short_perc=(short_count/num_rows)*100

print("There are " + str(short_count) + " short ORFs")
print(str(short_perc)+ "% of all ORFs")
#adding contig info #####NEED TO CHANGE THE SEQSTATS TO -T AND THEN CHANGE THIS. 

with open(genome_stats) as g:
    line=g.readlines()
    for i in range(len(line)):
        if (i==1): 
         j=line[1].split("\t")
         contigs=j[3]
print(j)
print("There are " + str(contigs) + " contigs")

# adding mlplasmids 

#with open(mlplasmids) as m:
#    mlplas_count=0
#    l=m.readlines()
#    for g in range(len(l)):
#        if (g!=0):
#            k=l[g].split("\t")
#            if (k[2]=='"Plasmid"'):
#                mlplas_count += 1 
#print(str(mlplas_count) + " plasmids by mlplasmids")

mlplas_count=5

#adding plasmid_finder from json parsed text file
with open(plasmidfinder) as p:
    pfplas_count = 0 
    n=p.readlines()
    plasmid_list1=[]
    plasmid_list2=[]
    for h in range(len(n)):
        if (h!=0):
            o=n[h].split(",")
            if (o[4] == 'No plasmid found\n'):
                break;
            else:
                b=o[9].split(" ")
                plasmid_list1.append(b[0])
                for i in plasmid_list1:
                    if i not in plasmid_list2:
                        pfplas_count += 1 
                        plasmid_list2.append(i)
print(str(pfplas_count)+ " plasmids by plasmidfinder")

#Unidentified_contigs_by_method
mlplas_unID= (int(contigs)-1)-int(mlplas_count)
pfplas_unID= (int(contigs)-1)-int(pfplas_count)

print(str(mlplas_unID) + " UNid'd contigs by mlplasmids")
print(str(pfplas_unID)+ " UNid'd contigs by plasmidfinder")
# adding concordant reads 
with open(concordant_reads) as c:
    d=c.readlines()
    concord=d[0].strip("\n")
    frac_concord=int(concord)/10000 
print(str(frac_concord) + " Concord mapping reads") 

with open(non_mapping_reads) as c:
    d=c.readlines()
    non_mapping=d[0].strip("\n")
    frac_non_mapping=int(non_mapping)/10000 
print(str(frac_non_mapping) + " Unmapped reads")


with open(socru) as c:
    s=c.readlines()
    d=s[0].split("\t")
    j=d[3:]
    #j[-1]=j[-1].strip("\n")
    socru_state=" "
    if (j == ['1','2','3','4','5','6','7\n']):
        socru_state="Correct"
    elif (j== ["1'",'2','3','4','5','6','7\n']):
        socru_state="Known"
    else: 
        socru_state="Incorrect"
        

print(socru_state)

metrics=["strain","assembler","contigs","mlplasmids","plasmidfinder","short_ORF","perc_short_ORF","illumina_concords","illumina_non_mapping","socru_arrangement","socru_correct","socru_known"]

#creating dict for output to CSV 
results={}
results.update({"Strain": strain , "Assembler": assembler,"Contigs": contigs, "mlplasmids": mlplas_count,"Plasmidfinder": pfplas_count ,"UnID_mlplas": mlplas_unID,"UnID_plasfinder": pfplas_unID ,"Short_ORF": short_count ,"Perc_short_ORF": short_perc , "Illumina_concord": frac_concord,"Illumina_non_mapping": frac_non_mapping, "Socru_state": socru_state })

open("Assembly_metrics.tab",'a')

with open("Assembly_metrics.tab", "r+") as tabfile: 
    fieldnames= ["Strain","Assembler","Contigs","mlplasmids","Plasmidfinder","UnID_mlplas","UnID_plasfinder","Short_ORF","Perc_short_ORF","Illumina_concord","Illumina_non_mapping","Socru_state"]
    writer=csv.DictWriter(tabfile,fieldnames=fieldnames,delimiter='\t')
    i= tabfile.readlines()
    #print(i[0])
    if (i == []):
        writer.writeheader()
    writer.writerow(results)
    
#with open("Assembly_metrics.tab","r+") as header_check: 
#    fieldnames= ["Strain","Assembler","Contigs","mlplasmids","Plasmidfinder","UnID_mlplas","UnID_plasfinder","Short_ORF","Perc_short_ORF","Illumina_concord","Illumina_non_mapping","Socru_state"]
#    writer=csv.DictWriter(header_check,fieldnames=fieldnames,delimiter='\t')
#    i=header_check.readlines()
#    #print(i[0])
#    if (i[0] != " "):
#        writer.writerow(results) 



#print(results)
