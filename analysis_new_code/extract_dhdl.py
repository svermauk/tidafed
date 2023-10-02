# Python code to get the data
# Shivani Verma
# 7th June, 2022

fhand = open('ti001.out')
count=0
for line in fhand:
    if line.startswith(' DV/DL  ='):
        count=count+1
        if(count%2!=0):
            dvdl=line.strip()
            value=float(dvdl[8:])
            print(value)
