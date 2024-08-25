import os, fileinput, sys

fnum_dict = {2012: [12360, 12630, 12631, 12362],
             2015, [20174, 20178, 20179, 20180],
             2018, [22570, 22580, 22583, 22586]}

particles = ["p", "He", "O", "Fe"]

def runSim(particle, fnum):
    for line in fileinput.input('steeringcard_dstsim',inplace=1):
	arg_string = particle + ' ' + str(fnum)
	if 'arguments' in line:
	    line = line.replace(line, 'arguments = '+arg_string+'\n')
	sys.stdout.write(line)
    os.system('condor_submit steeringcard_dstsim')
    print('running for', particle, fnum)


for year in fnum_dict.keys():
    fnum = fnum_dict[year]
    for j in len(fnum):
        runSim(particles[j], fnum[j])
