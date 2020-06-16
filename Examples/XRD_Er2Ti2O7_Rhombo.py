import sys
sys.path.append('../')

from SingleXRD.CristalStructure import Cubic, Rhombohedral
from SingleXRD.Diffractometer import Diffractometer
from SingleXRD.Geometry import two_theta_2_wave_vector_length
from SingleXRD import wave_length

import numpy as np

ETO_prim = Rhombohedral(0.715491,60)

ETO_prim.add_atom((0.625000, 0.125000, 0.125000),68,'Er')
ETO_prim.add_atom((0.125000, 0.125000, 0.625000),68,'Er')
ETO_prim.add_atom((0.125000, 0.125000, 0.125000),68,'Er')
ETO_prim.add_atom((0.125000, 0.625000, 0.125000),68,'Er')
ETO_prim.add_atom((0.125000, 0.625000, 0.625000),22,'Ti')
ETO_prim.add_atom((0.625000, 0.625000, 0.125000),22,'Ti')
ETO_prim.add_atom((0.625000, 0.625000, 0.625000),22,'Ti')
ETO_prim.add_atom((0.625000, 0.125000, 0.625000),22,'Ti')
ETO_prim.add_atom((0.000000, 0.000000, 1.000000),8,'O')
ETO_prim.add_atom((0.250000, 0.250000, 0.250000),8,'O')
ETO_prim.add_atom((0.706533, 0.293467, 0.293467),8,'O')
ETO_prim.add_atom((0.956533, 0.543467, 0.543467),8,'O')
ETO_prim.add_atom((0.293467, 0.293467, 0.706533),8,'O')
ETO_prim.add_atom((0.706533, 0.293467, 0.706533),8,'O')
ETO_prim.add_atom((0.293467, 0.706533, 0.293467),8,'O')
ETO_prim.add_atom((0.706533, 0.706533, 0.293467),8,'O')
ETO_prim.add_atom((0.293467, 0.706533, 0.706533),8,'O')
ETO_prim.add_atom((0.956533, 0.543467, 0.956533),8,'O')
ETO_prim.add_atom((0.543467, 0.956533, 0.543467),8,'O')
ETO_prim.add_atom((0.543467, 0.956533, 0.956533),8,'O')
ETO_prim.add_atom((0.543467, 0.543467, 0.956533),8,'O')
ETO_prim.add_atom((0.956533, 0.956533, 0.543467),8,'O')

hkl_z = (1,1,1)

dm = Diffractometer(cristal_structure = ETO_prim, surface_normal_hkl = hkl_z)
# pqr_x = (1,0,0)

# # Show some peaks values
HKL = (2,2,2)
print('hkl : {0:d}{1:d}{2:d}, 2T : {3:.4f}, Source : {4:.4f}, Detector : {5:.4f}, Phi : {6:.4f}, Factor : {7:.3f}'.format(*HKL,dm.hkl_2_two_theta(HKL),*dm.hkl_2_theta_source_theta_detector(HKL),0,np.abs(ETO_prim.structure_factor(HKL))**2))

HKL = (6,4,4)
print('hkl : {0:d}{1:d}{2:d}, 2T : {3:.4f}, Source : {4:.4f}, Detector : {5:.4f}, Phi : {6:.4f}, Factor : {7:.3f}'.format(*HKL,dm.hkl_2_two_theta(HKL),*dm.hkl_2_theta_source_theta_detector(HKL),0,np.abs(ETO_prim.structure_factor(HKL))**2))

HKL = (4,4,2)
print('hkl : {0:d}{1:d}{2:d}, 2T : {3:.4f}, Source : {4:.4f}, Detector : {5:.4f}, Phi : {6:.4f}, Factor : {7:.3f}'.format(*HKL,dm.hkl_2_two_theta(HKL),*dm.hkl_2_theta_source_theta_detector(HKL),0,np.abs(ETO_prim.structure_factor(HKL))**2))



# Build unit cell. Using atomic number as atomic factor
# STO = Cubic(0.3905)
# STO.add_atom((0,0,0),38,label = 'Sr')
# STO.add_atom((0.5,0.5,0.5),22,label = 'Ti')
# STO.add_atom((0.5,0.5,0),8,label = 'O1')
# STO.add_atom((0.5,0,0.5),8,label = 'O2')
# STO.add_atom((0,0.5,0.5),8,label = 'O3')

# # Orientation of normal and in plane of the cristal
# hkl_z = (0,0,1)
# pqr_x = (1,0,0)

# # Put the sample in the diffractometer with the right orientation
# dm = Diffractometer(cristal_structure = STO, surface_normal_hkl = hkl_z,azimuth_pqr = pqr_x)

# # Show some peaks values
# HKL = (0,0,2)
# print('hkl : {0:d}{1:d}{2:d}, 2T : {3:.4f}, Source : {4:.4f}, Detector : {5:.4f}, Phi : {6:.4f}, Factor : {7:.3f}'.format(*HKL,dm.hkl_2_two_theta(HKL),*dm.hkl_2_theta_source_theta_detector_phi(HKL),STO.structure_factor(HKL)))

# HKL = (1,0,3)
# print('hkl : {0:d}{1:d}{2:d}, 2T : {3:.4f}, Source : {4:.4f}, Detector : {5:.4f}, Phi : {6:.4f}, Factor : {7:.3f}'.format(*HKL,dm.hkl_2_two_theta(HKL),*dm.hkl_2_theta_source_theta_detector_phi(HKL),STO.structure_factor(HKL)))

# HKL = (0,1,3)
# print('hkl : {0:d}{1:d}{2:d}, 2T : {3:.4f}, Source : {4:.4f}, Detector : {5:.4f}, Phi : {6:.4f}, Factor : {7:.3f}'.format(*HKL,dm.hkl_2_two_theta(HKL),*dm.hkl_2_theta_source_theta_detector_phi(HKL),STO.structure_factor(HKL)))
