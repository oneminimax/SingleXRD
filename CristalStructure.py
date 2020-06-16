import numpy as np

from SingleXRD.Wyckoff import Position

class GeneralStructure(object):

    """ General notation :

    a     : lattice parameters (ex : a1, a2, etc)
    angle : angles between lattice base vector (ex : angle3 is angle between a1 and a2)
    b     : reciprocal lattice parameters (ex : b1, b2, etc)
    r     : direct space vector
    pqr   : components of a direct space vector in the lattice base
    g     : reciprocal space vector
    hkl   : components of a reciprocal space vector in the reciprocal lattice base
    """

    base_wyckoff_positions = '(0,0,0)'
    wyckoff_positions = {
        'a' : '(x,y,z)'
    }
    def __init__(self,a_lengths,a_angles):

        '''
        a_lengths : lattice parameter in nm
        a_angles : lattice angles in deg (converted in rad) (23),(13),(12)

        '''

        self.a_lengths = np.array(a_lengths)
        self.a_angles = np.deg2rad(np.array(a_angles))
        self.angles_matrix = np.array([[0,self.a_angles[2],self.a_angles[1]],[self.a_angles[2],0,self.a_angles[0]],[self.a_angles[1],self.a_angles[0],0]])

        self._precalc_aa_dot_mat()
        self._precalc_bb_dot_mat()      

        self.atoms = list()

    def __str__(self):

        return 'General structure : \n a = {0:.4f}, b = {1:.4f}, c = {2:.4f}\n alpha = {3:.4f}, beta = {4:.4f}, gamma = {5:.4f}'.format(*self.a_lengths,*np.rad2deg(self.a_angles))

    def list_atoms(self):

        for atom in self.atoms:
            print(atom)

    def get_atoms(self):

        return self.atoms

    def add_atom_wickoff(self,letter,form_factor,label = '',**kwargs):

        wyckoff_position = Position.generate_from_string(self.wyckoff_positions[letter],self.base_wyckoff_positions)
        coordinates = wyckoff_position(**kwargs)
        for coord in coordinates:
            self._add_atom(coord,form_factor,label)

    def add_atom(self,coordinates,form_factor,label = ''):

        self._add_atom(coordinates,form_factor,label)

    def _add_atom(self,coordinates,form_factor,label = ''):

        self.atoms.append(Atom(coordinates,label = label,form_factor = form_factor))

    def structure_factor(self,hkl):

        hkl = np.array(hkl)
        G = self.g_length(hkl)

        atoms = self.get_atoms()

        form_factors = np.zeros((len(atoms),))
        positions = np.zeros((3,len(atoms)))

        for i, atom in enumerate(atoms):
            form_factors[i] = atom.get_form_factor(G)
            positions[:,i] = atom.coordinates

        factor = np.sum(form_factors * np.exp(1j*2*np.pi*np.dot(hkl,positions)))

        return factor

    def debye_waller_factor(self,hkl):

        hkl = np.array(hkl)
        G = self.g_length(hkl)

    def a_length(self,i):

        '''
        Get lattice parameter, cyclic function
        '''

        return self.a_lengths[i%3]

    def a_a_angle(self,i,j):

        '''
        Get angle between lattice vectors
        '''

        return self.angles_matrix[i%3,j%3]

    def a_cell_volume(self):

        '''
        Unit cell volume
        '''

        return np.prod(self.a_lengths)*np.sqrt(1 + 2*np.prod(np.cos(self.a_angles)) - np.sum(np.cos(self.a_angles)**2))

    def r_length(self,pqr):

        '''
        Length of a direct vector in PQR coordinates
        '''

        return np.sqrt(self.r_r_dot(pqr,pqr))

    def r_r_angle(self,pqr,pqrp):

        '''
        Get angle between two direct vectors
        '''

        return np.rad2deg(np.arccos(self.r_r_dot(pqr,pqrp) / np.sqrt( self.r_r_dot(pqr,pqr) * self.r_r_dot(pqrp,pqrp) )))

    def b_length(self,i):


        '''
        Get reciprocal lattice parameter, cyclic function
        '''

        i = i%3
        return np.sqrt(self.bb_dot_mat[i,i])

    def b_b_angle(self,i,j):

        '''
        Get angle between reciprocal lattice vectors
        '''

        return np.rad2deg(np.arccos(self.bb_dot_mat[i,j] * np.sqrt(1/(self.bb_dot_mat[i,i]*self.bb_dot_mat[j,j]))))

    def g_length(self,hkl):

        '''
        Get length of reciprocal vector
        '''

        return np.sqrt(self.g_g_dot(hkl,hkl))

    def g_g_angle(self,hkl1,hkl2):

        '''
        Get angle between two reciprocal vectors
        '''

        return np.rad2deg(np.arccos(self.g_g_dot(hkl1,hkl2) / np.sqrt( self.g_g_dot(hkl1,hkl1) * self.g_g_dot(hkl2,hkl2) )))

    def r_r_angle(self,pqr1,pqr2):

        '''
        Get angle between two direct vectors
        '''

        return np.rad2deg(np.arccos(self.r_r_dot(pqr1,pqr2) / np.sqrt( self.r_r_dot(pqr1,pqr1) * self.r_r_dot(pqr2,pqr2) )))

    def g_r_angle(self,hkl,pqr):

        '''
        Get angle between two reciprocal/direct vectors
        '''

        return np.rad2deg(np.arccos(self.g_r_dot(hkl,pqr) / np.sqrt( self.g_g_dot(hkl,hkl) * self.r_r_dot(pqr,pqr) )))

    def hkl_planes_spacing(self,hkl):

        '''
        Spacing between plane of familly hkl
        '''

        return 2*np.pi/self.g_length(hkl)

    def _precalc_aa_dot_mat(self):

        self.aa_dot_mat = self.a_lengths[:,None] * self.a_lengths[None,:] * np.cos(self.angles_matrix)

    def _precalc_bb_dot_mat(self):

        aa_dot_mat = self.aa_dot_mat

        aaaa_cdc = aa_dot_mat[:,:,None,None] * aa_dot_mat[None,None,:,:] - aa_dot_mat[:,None,None,:] * aa_dot_mat[None,:,:,None]

        eijk = np.zeros((3, 3, 3))
        eijk[0, 1, 2] = eijk[1, 2, 0] = eijk[2, 0, 1] = 1
        eijk[0, 2, 1] = eijk[2, 1, 0] = eijk[1, 0, 2] = -1

        self.bb_dot_mat = (np.pi/self.a_cell_volume())**2 * np.einsum('imn,jmn',eijk,np.einsum('jpq,mpnq',eijk,aaaa_cdc))

    def g_g_dot(self,hkl1,hkl2):

        hkl1 = np.array(hkl1)
        hkl2 = np.array(hkl2)

        prod = np.sum(hkl1[:,None]*hkl2[None,:]*self.bb_dot_mat)

        return prod

    def r_r_dot(self,pqr1,pqr2):

        pqr1 = np.array(pqr1)
        pqr2 = np.array(pqr2)

        prod = np.sum(pqr1[:,None]*pqr2[None,:]*self.aa_dot_mat)

        return prod

    def g_r_dot(self,hkl,pqr):

        hkl = np.array(hkl)
        pqr = np.array(pqr)

        prod = 2 * np.pi * np.dot(hkl,pqr)

        return prod

    def g_g_cross(self,hkl1,hkl2):

        hkl1 = np.array(hkl1)
        hkl2 = np.array(hkl2)

        eijk = np.zeros((3, 3, 3))
        eijk[0, 1, 2] = eijk[1, 2, 0] = eijk[2, 0, 1] = 1
        eijk[0, 2, 1] = eijk[2, 1, 0] = eijk[1, 0, 2] = -1

        hh = hkl1[:,None] * hkl2[None,:]

        pqr = (2 * np.pi)**2 / self.a_cell_volume() * np.einsum('ijk,ij',eijk,hh)

        return pqr

    def r_r_cross(self,pqr1,pqr2):

        pqr1 = np.array(pqr1)
        pqr2 = np.array(pqr2)

        eijk = np.zeros((3, 3, 3))
        eijk[0, 1, 2] = eijk[1, 2, 0] = eijk[2, 0, 1] = 1
        eijk[0, 2, 1] = eijk[2, 1, 0] = eijk[1, 0, 2] = -1

        pp = pqr1[:,None] * pqr2[None,:]

        hkl = self.a_cell_volume()/(2*np.pi) * np.einsum('ijk,ij',eijk,pp)

        return hkl

    def proj_g_on_g(self,hkl1,hkl2):

        hkl1 = np.array(hkl1)
        hkl2 = np.array(hkl2)

        return self.g_g_dot(hkl1,hkl2)/self.g_g_dot(hkl2,hkl2) * hkl2

    def proj_r_on_r(self,pqr1,pqr2):

        pqr1 = np.array(pqr1)
        pqr2 = np.array(pqr2)

        return self.r_r_dot(pqr1,pqr2)/self.r_r_dot(pqr2,pqr2) * pqr2

    def proj_g_on_r(self,hkl,pqr):

        hkl = np.array(hkl)
        pqr = np.array(pqr)

        return self.g_r_dot(hkl,pqr)/self.r_r_dot(pqr,pqr) * pqr

    def proj_r_on_g(self,pqr,hkl):

        pqr = np.array(pqr)
        hkl = np.array(hkl)

        return self.g_r_dot(hkl,pqr)/self.g_g_dot(hkl,hkl) * hkl

    def r_to_g(self,pqr):

        pqr = np.array(pqr)

        return 1/(2*np.pi)*np.dot(self.aa_dot_mat,pqr)

    def g_to_r(self,hkl):

        hkl = np.array(hkl)

        return 1/(2*np.pi)*np.dot(self.bb_dot_mat,hkl)

    def g_g_angle_on_plane_hkl(self,hkl1,hkl2,hkl_plane):

        '''
        Get angle between two reciprocal vectors projected on a plane perpendicular to vector hkl_plane
        '''

        hkl_plane = np.array(hkl_plane)
        hkl1 = np.array(hkl1)
        hkl2 = np.array(hkl2)

        hkl1_onplane = hkl1 - self.proj_g_on_g(hkl1,hkl_plane)
        hkl2_onplane = hkl2 - self.proj_g_on_g(hkl2,hkl_plane)

        angle = self.g_g_angle(hkl1_onplane,hkl2_onplane)
        direction = np.sign(self.g_g_dot(self.g_g_cross(hkl1,hkl2),hkl_plane))

        return angle * direction

    def g_r_angle_on_plane_hkl(self,hkl,pqr,hkl_plane):

        '''
        Get angle between two vectors (reciprocal and direct) projected on a plane perpendicular to vector hkl_plane
        '''

        hkl_plane = np.array(hkl_plane)
        hkl = np.array(hkl)
        pqr = np.array(pqr)

        hkl_onplane = hkl - self.proj_g_on_g(hkl,hkl_plane)
        pqr_onplane = pqr - self.g_to_r(self.proj_r_on_g(pqr,hkl_plane))

        angle = self.g_r_angle(hkl_onplane,pqr_onplane)
        direction = np.sign(self.g_g_dot(self.g_g_cross(hkl,self.r_to_g(pqr)),hkl_plane))

        return angle * direction

class Rhombohedral(GeneralStructure):
    def __init__(self,a_length,angle):

        super().__init__((a_length,a_length,a_length),(angle,angle,angle))

class Hexagonal(GeneralStructure):
    def __init__(self,a12_length,a3_length):

        super().__init__((a12_length,a12_length,a3_length),(90,90,120))

class Monoclinic(GeneralStructure):
    def __init__(self,a1_length,a2_length,a3_length,angle2):

        super().__init__((a1_length,a2_length,a3_length),(90,angle2,90))

class Tetragonal(GeneralStructure):
    def __init__(self,a12_length,a3_length):

        super().__init__((a12_length,a12_length,a3_length),(90,90,90))

class Orthorhombic(GeneralStructure):
    def __init__(self,a1_length,a2_length,a3_length):

        super().__init__((a1_length,a2_length,a3_length),(90,90,90))

class Cubic(GeneralStructure):
    def __init__(self,a_length):

        super().__init__((a_length,a_length,a_length),(90,90,90))

class CubicFaceCentered(Cubic):
    
    def add_atom(self,coordinates,form_factor,label = ''):

        for v1 in [np.array((0,0,0)),np.array((0.5,0.5,0)),np.array((0,0.5,0.5)),np.array((0.5,0,0.5))]:
            self._add_atom(coordinates + v1,form_factor,label)

class CubicBodyCentered(Cubic):
    
    def add_atom(self,coordinates,form_factor,label = ''):

        for v1 in [np.array((0,0,0)),np.array((0.5,0.5,0.5))]:
            self._add_atom(coordinates + v1,form_factor,label)

class Diamond(Cubic):

    def add_atom(self,coordinates,form_factor,label = ''):

        for v1 in [np.array((0,0,0)),np.array((0.5,0.5,0)),np.array((0,0.5,0.5)),np.array((0.5,0,0.5))]:
                for v2 in [np.array((0,0,0)),np.array((0.25,0.25,0.25))]:
                    self._add_atom(coordinates + v1 + v2,form_factor,label)

class No139(Tetragonal): # I 4 / m m m
    base_wyckoff_positions = '(0,0,0) (1/2,1/2,1/2)'
    wyckoff_positions = {
        'a' : '(0,0,0)',
        'b' : '(0,0,1/2)',
        'c' : '(0,1/2,0) (1/2,0,0)',
        'd' : '(0,1/2,1/4) (1/2,0,1/4)',
        'e' : '(0,0,z) (0,0,-z)',
        'f' : '(1/4,1/4,1/4) (3/4,3/4,1/4) (3/4,1/4,1/4) (1/4,3/4,1/4)',
        'g' : '(0,1/2,z) (1/2,0,z) (0,1/2,-z) (1/2,0,-z)',
        'h' : '(x,x,0) (-x,-x,0) (-x,x,0) (x,-x,0)',
        'i' : '(x,0,0) (-x,0,0) (0,x,0) (0,-x,0)',
        'j' : '(x,1/2,0) (-x,1/2,0) (1/2,x,0) (1/2,-x,0)',
        'k' : '(x,x+1/2,1/4) (-x,-x+1/2,1/4) (-x+1/2,x,1/4) (x+1/2,-x,1/4) (-x,-x+1/2,3/4) (x,x+1/2,3/4) (x+1/2,-x,3/4) (-x+1/2,x,3/4)',
        'l' : '(x,y,0) (-x,-y,0) (-y,x,0) (y,-x,0) (-x,y,0) (x,-y,0) (y,x,0) (-y,-x,0)',
        'm' : '(x,x,z) (-x,-x,z) (-x,x,z) (x,-x,z) (-x,x,-z) (x,-x,-z) (x,x,-z) (-x,-x,-z)',
        'n' : '(0,y,z) (0,-y,z) (-y,0,z) (y,0,z) (0,y,-z) (0,-y,-z) (y,0,-z) (-y,0,-z)',
        'o' : '(x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (-x,y,-z) (x,-y,-z) (y,x,-z) (-y,-x,-z) (-x,-y,-z) (x,y,-z) (y,-x,-z) (-y,x,-z) (x,-y,z) (-x,y,z) (-y,-x,z) (y,x,z)'
    }

class No206(Cubic): # I a -3
    base_wyckoff_positions = '(0,0,0) (1/2,1/2,1/2)'
    wyckoff_positions = {
        'a' : '(0,0,0) (1/2,0,1/2) (0,1/2,1/2) (1/2,1/2,0)',
        'b' : '(1/4,1/4,1/4) (1/4,3/4,3/4) (3/4,3/4,1/4) (3/4,1/4,3/4)',
        'c' : '(x,x,x) (-x+1/2,-x,x+1/2) (-x,x+1/2,-x+1/2) (x+1/2,-x+1/2,-x) (-x,-x,-x) (x+1/2,x,-x+1/2) (x,-x+1/2,x+1/2) (-x+1/2,x+1/2,x)',
        'd' : '(x,0,1/4) (-x+1/2,0,3/4) (1/4,x,0) (3/4,-x+1/2,0) (0,1/4,x) (0,3/4,-x+1/2) (-x,0,3/4) (x+1/2,0,1/4) (3/4,-x,0) (1/4,x+1/2,0) (0,3/4,-x) (0,1/4,x+1/2)',
        'e' : '(x,y,z) (-x+1/2,-y,z+1/2) (-x,y+1/2,-z+1/2) (x+1/2,-y+1/2,-z) (z,x,y) (z+1/2,-x+1/2,-y) (-z+1/2,-x,y+1/2) (-z,x+1/2,-y+1/2) (y,z,x) (-y,z+1/2,-x+1/2) (y+1/2,-z+1/2,-x) (-y+1/2,-z,x+1/2) (-x,-y,-z) (x+1/2,y,-z+1/2) (x,-y+1/2,z+1/2) (-x+1/2,y+1/2,z) (-z,-x,-y) (-z+1/2,x+1/2,y) (z+1/2,x,-y+1/2) (z,-x+1/2,y+1/2) (-y,-z,-x) (y,-z+1/2,x+1/2) (-y+1/2,z+1/2,x) (y+1/2,z,-x+1/2)',
    }

class No221(Cubic): # P m -3 m
    base_wyckoff_positions = '(0,0,0)'
    wyckoff_positions = {
        'a' : '(0,0,0)',
        'b' : '(1/2,1/2,1/2)',
        'c' : '(0,1/2,1/2) (1/2,0,1/2) (1/2,1/2,0)',
        'd' : '(1/2,0,0) (0,1/2,0) (0,0,1/2)',
        'e' : '(x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x)',
        'f' : '(x,1/2,1/2) (-x,1/2,1/2) (1/2,x,1/2) (1/2,-x,1/2) (1/2,1/2,x) (1/2,1/2,-x)',
        'g' : '(x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (x,x,-x) (-x,-x,-x) (x,-x,x) (-x,x,x)',
        'h' : '(x,1/2,0) (-x,1/2,0) (0,x,1/2) (0,-x,1/2) (1/2,0,x) (1/2,0,-x) (1/2,x,0) (1/2,-x,0) (x,0,1/2) (-x,0,1/2) (0,1/2,-x) (0,1/2,x)',
        'i' : '(0,y,y) (0,-y,y) (0,y,-y) (0,-y,-y) (y,0,y) (y,0,-y) (-y,0,y) (-y,0,-y) (y,y,0) (-y,y,0) (y,-y,0) (-y,-y,0)',
        'j' : '(1/2,y,y) (1/2,-y,y) (1/2,y,-y) (1/2,-y,-y) (y,1/2,y) (y,1/2,-y) (-y,1/2,y) (-y,1/2,-y) (y,y,1/2) (-y,y,1/2) (y,-y,1/2) (-y,-y,1/2)',
        'k' : '(0,y,z) (0,-y,z) (0,y,-z) (0,-y,-z) (z,0,y) (z,0,-y) (-z,0,y) (-z,0,-y) (y,z,0) (-y,z,0) (y,-z,0) (-y,-z,0) (y,0,-z) (-y,0,-z) (y,0,z) (-y,0,z) (0,z,-y) (0,z,y) (0,-z,-y) (0,-z,y) (z,y,0) (z,-y,0) (-z,y,0) (-z,-y,0)',
        'l' : '(1/2,y,z) (1/2,-y,z) (1/2,y,-z) (1/2,-y,-z) (z,1/2,y) (z,1/2,-y) (-z,1/2,y) (-z,1/2,-y) (y,z,1/2) (-y,z,1/2) (y,-z,1/2) (-y,-z,1/2) (y,1/2,-z) (-y,1/2,-z) (y,1/2,z) (-y,1/2,z) (1/2,z,-y) (1/2,z,y) (1/2,-z,-y) (1/2,-z,y) (z,y,1/2) (z,-y,1/2) (-z,y,1/2) (-z,-y,1/2)',
        'm' : '(x,x,z) (-x,-x,z) (-x,x,-z) (x,-x,-z) (z,x,x) (z,-x,-x) (-z,-x,x) (-z,x,-x) (x,z,x) (-x,z,-x) (x,-z,-x) (-x,-z,x) (x,x,-z) (-x,-x,-z) (x,-x,z) (-x,x,z) (x,z,-x) (-x,z,x) (-x,-z,-x) (x,-z,x) (z,x,-x) (z,-x,x) (-z,x,x) (-z,-x,-x)',
        'n' : '(x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y,x,-z) (-y,-x,-z) (y,-x,z) (-y,x,z) (x,z,-y) (-x,z,y) (-x,-z,-y) (x,-z,y) (z,y,-x) (z,-y,x) (-z,y,x) (-z,-y,-x) (-x,-y,-z) (x,y,-z) (x,-y,z) (-x,y,z) (-z,-x,-y) (-z,x,y) (z,x,-y) (z,-x,y) (-y,-z,-x) (y,-z,x) (-y,z,x) (y,z,-x) (-y,-x,z) (y,x,z) (-y,x,-z) (y,-x,-z) (-x,-z,y) (x,-z,-y) (x,z,y) (-x,z,-y) (-z,-y,x) (-z,y,-x) (z,-y,-x) (z,y,x)',
    }

class No225(Cubic): # F m -3 m
    base_wyckoff_positions = '(0,0,0) (0,1/2,1/2) (1/2,0,1/2) (1/2,1/2,0)'
    wyckoff_positions = {
        'a' : '(0,0,0)',
        'b' : '(1/2,1/2,1/2)',
        'c' : '(1/4,1/4,1/4) (1/4,1/4,3/4)',
        'd' : '(0,1/4,1/4) (0,3/4,1/4) (1/4,0,1/4) (1/4,0,3/4) (1/4,1/4,0) (3/4,1/4,0)',
        'e' : '(0,1/4,1/4) (0,3/4,1/4) (1/4,0,1/4) (1/4,0,3/4) (1/4,1/4,0) (3/4,1/4,0)',
        'f' : '(x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (x,x,-x) (-x,-x,-x) (x,-x,x) (-x,x,x)',
        'g' : '(x,1/4,1/4) (-x,3/4,1/4) (1/4,x,1/4) (1/4,-x,3/4) (1/4,1/4,x) (3/4,1/4,-x) (1/4,x,3/4) (3/4,-x,3/4) (x,1/4,3/4) (-x,1/4,1/4) (1/4,1/4,-x) (1/4,3/4,x) ',
        'h' : '(0,y,y) (0,-y,y) (0,y,-y) (0,-y,-y) (y,0,y) (y,0,-y) (-y,0,y) (-y,0,-y) (y,y,0) (-y,y,0) (y,-y,0) (-y,-y,0)',
        'i' : '(1/2,y,y) (1/2,-y,y) (1/2,y,-y) (1/2,-y,-y) (y,1/2,y) (y,1/2,-y) (-y,1/2,y) (-y,1/2,-y) (y,y,1/2) (-y,y,1/2) (y,-y,1/2) (-y,-y,1/2)',
        'j' : '(0,y,z) (0,-y,z) (0,y,-z) (0,-y,-z) (z,0,y) (z,0,-y) (-z,0,y) (-z,0,-y) (y,z,0) (-y,z,0) (y,-z,0) (-y,-z,0) (y,0,-z) (-y,0,-z) (y,0,z) (-y,0,z) (0,z,-y) (0,z,y) (0,-z,-y) (0,-z,y) (z,y,0) (z,-y,0) (-z,y,0) (-z,-y,0)',
        'k' : '(x,x,z) (-x,-x,z) (-x,x,-z) (x,-x,-z) (z,x,x) (z,-x,-x) (-z,-x,x) (-z,x,-x) (x,z,x) (-x,z,-x) (x,-z,-x) (-x,-z,x) (x,x,-z) (-x,-x,-z) (x,-x,z) (-x,x,z) (x,z,-x) (-x,z,x) (-x,-z,-x) (x,-z,x) (z,x,-x) (z,-x,x) (-z,x,x) (-z,-x,-x)',
        'l' : '(x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y,x,-z) (-y,-x,-z) (y,-x,z) (-y,x,z) (x,z,-y) (-x,z,y) (-x,-z,-y) (x,-z,y) (z,y,-x) (z,-y,x) (-z,y,x) (-z,-y,-x) (-x,-y,-z) (x,y,-z) (x,-y,z) (-x,y,z) (-z,-x,-y) (-z,x,y) (z,x,-y) (z,-x,y) (-y,-z,-x) (y,-z,x) (-y,z,x) (y,z,-x) (-y,-x,z) (y,x,z) (-y,x,-z) (y,-x,-z) (-x,-z,y) (x,-z,-y) (x,z,y) (-x,z,-y) (-z,-y,x) (-z,y,-x) (z,-y,-x) (z,y,x)'
    }



class Atom(object):

    def __init__(self,coordinates,label = '',form_factor = 'ITC'):

        self.coordinates = np.mod(np.array(coordinates),1)
        self.form_factor = form_factor
        if form_factor == 'ITC':
            self.formFactorFct = ITCFct(label)
        self.label = label

    def __str__(self):

        return '{0:s} : ({1:.3f},{2:.3f},{3:.3f}) : f = {4:.4f}'.format(self.label,*self.coordinates,self.form_factor)

    def get_form_factor(self,q):

        if self.form_factor == 'ITC':
            return self.formFactorFct(q)
        else:
            return self.form_factor

def GA_GB_GC_to_cristalStructure(GA,GB,GC):

    VG = np.abs(np.dot(GA,np.cross(GB,GC)))

    Ra = 2*np.pi/VG * np.cross(GB,GC)
    Rb = 2*np.pi/VG * np.cross(GC,GA)
    Rc = 2*np.pi/VG * np.cross(GA,GB)

    a_length,b_length,a3_length = np.sqrt(np.sum(Ra**2)),np.sqrt(np.sum(Rb**2)),np.sqrt(np.sum(Rc**2))
    alpha = np.rad2deg(np.arccos(np.dot(Rb,Rc)/(b_length*a3_length)))
    beta  = np.rad2deg(np.arccos(np.dot(Ra,Rc)/(a_length*a3_length)))
    gamma = np.rad2deg(np.arccos(np.dot(Ra,Rb)/(a_length*b_length)))

    return GeneralStructure([a_length,b_length,a3_length],[alpha, beta, gamma])




