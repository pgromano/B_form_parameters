from glob import glob

f = open('settime.inp', 'w')
files = sorted(glob('*_*0.xtc'))
for file in files:
	f.write('%s\n' %file[3:-4])
f.close()
