#
import os
import glob
#
# command line is line
#
# -dispose clear images for next image (so titles don't overwrite)
# -delay creates a 2 second delay between images
# *.gif is create a give movie
# gm convert -dispose previous -delay 200 *.png heat_wave_recent.gif
# 
# how do we put in a file list for the *.png ?
#
filepath='/home/flbahr/heat_content/'
filepre='heat_wave_sst_'
fileend='_refine_wcofs.png'
files=glob.glob(filepath+filepre+'*'+fileend)
files.sort(key=os.path.getmtime)
lf=len(files)
thelist=[]
for k in range(lf-20,lf):
    thelist.append(files[k])
   # thelist.append('{}.png'.format(k))
command='gm convert -dispose previous -delay 200 {} /home/flbahr/heat_content/heat_wave_recent.gif'.format(' '.join(thelist))
os.system(command)
