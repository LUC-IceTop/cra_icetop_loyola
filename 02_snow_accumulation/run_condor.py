import os, fileinput, sys, argparse

if __name__ == "__main__":
    # Set up command line options
    parser = argparse.ArgumentParser(description=
                                     "Process burnsample hdf5 files into numpy arrays.")
    parser.add_argument("-i", "--indir", dest="indir", type=str, default=None,
                        help="Input directory with burnsample hdf5 files \
                        e.g. /data/user/gagrawal/dst_data/10yburnsample")
    parser.add_argument("-o", "--output", dest="output", type=str, default=None,
                        help="Output directory for numpy arrays")
    parser.add_argument("-l", "--log", dest="log", type=str,
                        default=None, help="Directory for condor output/error/log files")
    args = parser.parse_args()

    # Create directories where necessary
    #if not os.path.isdir(args.indir):
        #Throw exception?
    if not os.path.isdir(args.output):
        print('Creating output directory {}'.format(args.output))
        os.makedirs(args.output)
    if not os.path.isdir(args.log):
        print('Creating condor output directory {}'.format(args.log))
        os.makedirs(args.log)

    for line in fileinput.input('steeringcard_dst',inplace=1):
        if 'output =' in line:
            line = line.replace(line,'output = '+ args.log +'/condor.out.$(ClusterId)\n')
        elif 'error =' in line:
            line = line.replace(line,'error = '+ args.log +'/condor.err.$(ClusterId)\n')
        elif 'log =' in line:
            line = line.replace(line,'log = '+ args.log +'/condor.log.$(ClusterId)\n')
        sys.stdout.write(line) 
        
    for run_year in range(2011, 2022):         
        for line in fileinput.input('steeringcard_dst',inplace=1):
            arg_string = '{year} {i} {o}'.format(year = run_year, i = args.indir, o = args.output)
            if 'arguments' in line:
                line = line.replace(line,'arguments = '+arg_string+'\n')
            sys.stdout.write(line)
        os.system('condor_submit steeringcard_dst') 
        print('running for', run_year)