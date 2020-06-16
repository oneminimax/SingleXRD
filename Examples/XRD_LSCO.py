import sys
sys.path.append('../')

from SingleXRD.CristalStructure import Tetragonal
from SingleXRD.Diffractometer import Diffractometer

import numpy as np

# Build unit cell. Using atomic number as atomic factor
LSCO = Tetragonal(0.3775,1.323)
# LSCO.add_atom((0,0,0),1,label = 'LSCO')

# Orientation of normal and in plane of the cristal
hkl_z = (1,0,0)
pqr_x = (0,1,0)

# Put the sample in the diffractometer with the right orientation
dm = Diffractometer(cristal_structure = LSCO, surface_normal_hkl = hkl_z,azimuth_pqr = pqr_x)

# Show some peaks values
HKL = (2,0,0)
print('hkl : {0:d}{1:d}{2:d}, 2T : {3:.4f}, Source : {4:.4f}, Detector : {5:.4f}, Phi : {6:.4f}, Factor : {7:.3f}'.format(*HKL,dm.hkl_2_two_theta(HKL),*dm.hkl_2_theta_source_theta_detector_phi(HKL),LSCO.structure_factor(HKL)))

HKL = (3,1,0)
print('hkl : {0:d}{1:d}{2:d}, 2T : {3:.4f}, Source : {4:.4f}, Detector : {5:.4f}, Phi : {6:.4f}, Factor : {7:.3f}'.format(*HKL,dm.hkl_2_two_theta(HKL),*dm.hkl_2_theta_source_theta_detector_phi(HKL),LSCO.structure_factor(HKL)))

HKL = (3,0,1)
print('hkl : {0:d}{1:d}{2:d}, 2T : {3:.4f}, Source : {4:.4f}, Detector : {5:.4f}, Phi : {6:.4f}, Factor : {7:.3f}'.format(*HKL,dm.hkl_2_two_theta(HKL),*dm.hkl_2_theta_source_theta_detector_phi(HKL),LSCO.structure_factor(HKL)))
