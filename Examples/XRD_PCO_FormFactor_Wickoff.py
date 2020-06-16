import sys
sys.path.append('../')

from SingleXRD.CristalStructure import No139
from SingleXRD.Diffractometer import Diffractometer

import numpy as np

# Build unit cell. Using atomic number as atomic factor
PCO = No139(0.399573,1.248272)

PCO.add_atom_wickoff('e',59,label = 'Pr',z = 0.35038)
PCO.add_atom_wickoff('a',29,label = 'Cu')
PCO.add_atom_wickoff('c',8,label = 'O1')
PCO.add_atom_wickoff('d',8,label = 'O2')

# Orientation of normal and in plane of the cristal
hkl_z = (1,0,0)
pqr_x = (0,1,0)

# Put the sample in the diffractometer with the right orientation
dm = Diffractometer(cristal_structure = PCO, surface_normal_hkl = hkl_z,azimuth_pqr = pqr_x)

# Show some peaks values
HKL = (1,0,0)
print('hkl : {0:d}{1:d}{2:d}, 2T : {3:.4f}, Source : {4:.4f}, Detector : {5:.4f}, Phi : {6:.4f}, Factor : {7:.3f}'.format(*HKL,dm.hkl_2_two_theta(HKL),*dm.hkl_2_theta_source_theta_detector_phi(HKL),PCO.structure_factor(HKL)))

HKL = (2,0,0)
print('hkl : {0:d}{1:d}{2:d}, 2T : {3:.4f}, Source : {4:.4f}, Detector : {5:.4f}, Phi : {6:.4f}, Factor : {7:.3f}'.format(*HKL,dm.hkl_2_two_theta(HKL),*dm.hkl_2_theta_source_theta_detector_phi(HKL),PCO.structure_factor(HKL)))

HKL = (3,0,1)
print('hkl : {0:d}{1:d}{2:d}, 2T : {3:.4f}, Source : {4:.4f}, Detector : {5:.4f}, Phi : {6:.4f}, Factor : {7:.3f}'.format(*HKL,dm.hkl_2_two_theta(HKL),*dm.hkl_2_theta_source_theta_detector_phi(HKL),PCO.structure_factor(HKL)))
