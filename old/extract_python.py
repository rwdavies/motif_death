## do not generalize for now...
import gzip
import math
workdir="/well/myers/rwdavies/primates/mapping/reddeer/"
file1 = workdir + "/SRR4013902_1.fastq.gz"
file2 = workdir + "/SRR4013902_2.fastq.gz"
nmax=math.inf

f1 = gzip.open(file1, 'rb')
f2 = gzip.open(file2, 'rb')

out_files_1=[]
out_files_2=[]
names_out=[]

to_break=0;
count=0
while to_break == 0:
    count += 1
    if count > nmax:
        to_break=1
    ## get quartests of lines
    a1=f1.readline()
    if a1 == b"":
        to_break=1
        break
    a2=f1.readline()
    a3=f1.readline()
    a4=f1.readline()
    b1=f2.readline()
    b2=f2.readline()
    b3=f2.readline()
    b4=f2.readline()
    ## which output file?
    c=str.split(str(a1), " ")[1]
    flowcell=str.split(c, ":")[0]
    lane=str.split(c, ":")[1]
    out_name=flowcell + "_" + lane
    ##
    if (out_name in names_out) == False:
        out_f1=workdir + flowcell + "-" + lane + "_1.fastq.gz"
        out_f2=workdir + flowcell + "-" + lane + "_2.fastq.gz"
        out_files_1.append(gzip.open(out_f1, 'wb'))
        out_files_2.append(gzip.open(out_f2, 'wb'))    
        names_out.append(out_name)
    i_file=names_out.index(out_name)
    out_files_1[i_file].write(a1)
    out_files_1[i_file].write(a2)
    out_files_1[i_file].write(a3)
    out_files_1[i_file].write(a4)
    out_files_2[i_file].write(b1)
    out_files_2[i_file].write(b2)
    out_files_2[i_file].write(b3)
    out_files_2[i_file].write(b4)


        
for out_name in names_out:
    i_file=names_out.index(out_name)
    out_files_1[i_file].close()
    out_files_2[i_file].close()


