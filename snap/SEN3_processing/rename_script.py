"""
This routine shortens the filenames of Sentinel-3 files and stores the actual filename in a textfile

"""

import sys
import os

direc = 'E:/Lake_Erie_HAB/Data/remote_sensing_data/Sentinel/2016'
fnames = os.listdir(direc)


sr = 0; # serial number
new_name_list = []
for fname in fnames:
    direc1 = direc+'/'+fname
    fname1 = os.listdir(direc1)
    x = fname1[0].split('_')

    # assign new name
    sr = sr+1
    new_name_list.append(x[0]+'_'+x[1]+'_'+x[2]+'_'+x[3]+'_'+x[7]+'_'+str(sr)+'.SEN3')

    # rename the file
    old_name = direc1+'/'+fname1[0]
    new_name = direc1+'/'+new_name_list[sr-1]
    os.rename(old_name,new_name)

    # move the file one directory up
    old_path = new_name
    new_path = direc+'/'+new_name_list[sr-1]
    os.replace(old_path,new_path)
    

# write original name and corresponding new names in a textfile
wfname = 'correspondence_original_new_name.txt'
wfilename = direc+'/'+wfname

fid = open(wfilename,'w+')
fid.write('New_file_name\t\t\tOriginal_name\n')
for ind in range(len(fnames)):
    fid.write(new_name_list[ind]+'\t\t\t'+fnames[ind]+'\n')
fid.close()

