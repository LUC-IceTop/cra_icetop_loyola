#!/usr/bin/env python3
import numpy as np
import h5py, glob
import os
import sys
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
parser.add_argument('run_year', type=int)
parser.add_argument('indir', type=str)
parser.add_argument('outdir', type=str)
opts = parser.parse_args()

run_year = opts.run_year 
path_to_hdf5_files = opts.indir
path_to_save_files = opts.outdir 

if not os.path.isdir(path_to_save_files):
        print('Creating output directory {}'.format(path_to_save_files))
        os.makedirs(path_to_save_files)

globals()['s125_'+str(run_year)] = []
globals()['stations_'+str(run_year)] = []
globals()['zenith_'+str(run_year)] = []
globals()['azimuth_'+str(run_year)] = []
globals()['ShowerPlaneZen_'+str(run_year)] = []
globals()['ShowerPlaneAz_'+str(run_year)] = []
for file in glob.glob('{}/l3_data_run_config_{}_*.hdf5'.format(path_to_hdf5_files, run_year)):
    f = h5py.File(file, 'r')
    globals()['s125_'+str(run_year)].extend(f['LaputopParams']['s125'])
    globals()['stations_'+str(run_year)].extend(f['NStations']['value'])
    globals()['zenith_'+str(run_year)].extend(f['Laputop']['zenith'])
    globals()['azimuth_'+str(run_year)].extend(f['Laputop']['azimuth'])
    globals()['ShowerPlaneZen_'+str(run_year)].extend(f['ShowerPlane']['zenith'])
    globals()['ShowerPlaneAz_'+str(run_year)].extend(f['ShowerPlane']['azimuth'])    
np.save('{}/s125_{}'.format(path_to_save_files, str(run_year)), globals()['s125_'+str(run_year)])
np.save('{}/stations_{}'.format(path_to_save_files, str(run_year)), globals()['stations_'+str(run_year)])
np.save('{}/zenith_{}'.format(path_to_save_files, str(run_year)), globals()['zenith_'+str(run_year)])
np.save('{}/azimuth_{}'.format(path_to_save_files, str(run_year)), globals()['azimuth_'+str(run_year)])
np.save('{}/ShowerPlaneZen_{}'.format(path_to_save_files, str(run_year)), globals()['ShowerPlaneZen_'+str(run_year)])
np.save('{}/ShowerPlaneAz_{}'.format(path_to_save_files, str(run_year)), globals()['ShowerPlaneAz_'+str(run_year)])
print(run_year, len(globals()['s125_'+str(run_year)]), len(globals()['stations_'+str(run_year)])) 


globals()['IceTop_reco_succeeded_'+str(run_year)] = []
globals()['Laputop_FractionContainment_'+str(run_year)] = []
globals()['StationDensity_passed_'+str(run_year)] = []
globals()['exists_'+str(run_year)] = []
globals()['BetaCutPassed_'+str(run_year)] = []
globals()['IceTopMaxSignalAbove6_'+str(run_year)] = []
globals()['IceTopMaxSignalInside_'+str(run_year)] = []
globals()['IceTopNeighbourMaxSignalAbove4_'+str(run_year)] = []
globals()['IceTop_StandardFilter_'+str(run_year)] = []
for file in glob.glob('{}/l3_data_run_config_{}_*.hdf5'.format(path_to_hdf5_files, run_year)):
    f = h5py.File(file, 'r')
    globals()['IceTop_reco_succeeded_'+str(run_year)].extend(f['IT73AnalysisIceTopQualityCuts']['IceTop_reco_succeeded'])
    globals()['Laputop_FractionContainment_'+str(run_year)].extend(f['IT73AnalysisIceTopQualityCuts']['Laputop_FractionContainment'])
    globals()['StationDensity_passed_'+str(run_year)].extend(f['IT73AnalysisIceTopQualityCuts']['StationDensity_passed'])
    globals()['exists_'+str(run_year)].extend(f['IT73AnalysisIceTopQualityCuts']['exists'])
    globals()['BetaCutPassed_'+str(run_year)].extend(f['IT73AnalysisIceTopQualityCuts']['BetaCutPassed'])
    globals()['IceTopMaxSignalAbove6_'+str(run_year)].extend(f['IT73AnalysisIceTopQualityCuts']['IceTopMaxSignalAbove6'])
    globals()['IceTopMaxSignalInside_'+str(run_year)].extend(f['IT73AnalysisIceTopQualityCuts']['IceTopMaxSignalInside'])
    globals()['IceTopNeighbourMaxSignalAbove4_'+str(run_year)].extend(f['IT73AnalysisIceTopQualityCuts']['IceTopNeighbourMaxSignalAbove4'])
    globals()['IceTop_StandardFilter_'+str(run_year)].extend(f['IT73AnalysisIceTopQualityCuts']['IceTop_StandardFilter'])
np.save('{}/IceTop_reco_succeeded_{}'.format(path_to_save_files, str(run_year)), globals()['IceTop_reco_succeeded_'+str(run_year)])
np.save('{}/Laputop_FractionContainment_{}'.format(path_to_save_files, str(run_year)), globals()['Laputop_FractionContainment_'+str(run_year)])
np.save('{}/StationDensity_passed_{}'.format(path_to_save_files, str(run_year)), globals()['StationDensity_passed_'+str(run_year)])
np.save('{}/exists_{}'.format(path_to_save_files, str(run_year)), globals()['exists_'+str(run_year)])
np.save('{}/BetaCutPassed_{}'.format(path_to_save_files, str(run_year)), globals()['BetaCutPassed_'+str(run_year)])
np.save('{}/IceTopMaxSignalAbove6_{}'.format(path_to_save_files, str(run_year)), globals()['IceTopMaxSignalAbove6_'+str(run_year)])
np.save('{}/IceTopMaxSignalInside_{}'.format(path_to_save_files, str(run_year)), globals()['IceTopMaxSignalInside_'+str(run_year)])
np.save('{}/IceTopNeighbourMaxSignalAbove4_{}'.format(path_to_save_files, str(run_year)), globals()['IceTopNeighbourMaxSignalAbove4_'+str(run_year)])
np.save('{}/IceTop_StandardFilter_{}'.format(path_to_save_files, str(run_year)), globals()['IceTop_StandardFilter_'+str(run_year)])
print(run_year, len(globals()['IceTop_reco_succeeded_'+str(run_year)]), len(globals()['exists_'+str(run_year)]))