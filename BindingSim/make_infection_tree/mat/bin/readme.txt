gentree.exe 
4 arguments: infile, smpno, starttime, endtime.
infile: input csv file name wothout '.csv'
usage: gentree [infile] [smpno] [starttime] [endtime]
example:
if you have a csv file 'voutput_small.csv' and the first line is the list of the column names
then use the command, gentree voutput_small 30 10 200, you will get the output

output:
Two types of the output files: 1. the tree figure, 2. the tree text files 
Tree figure file name: [infile]_tree_[smpno].jpg
Tree text file name: [infile]_tree_[smpno].branch.nex