import numpy as np

def two_theta_2_wave_vector_length(wave_length,two_theta):

    """ Convert two theta to wave vector length

    Inputs
    wave_length : x-ray source wave length (float)
    two_theta : angle between the source and the detector in degree (float, list, np.array)

    Output
    wave_vector_length : length of the wave vector of the reciprocal space in rad/units of wave_length input
    """

    return 4*np.pi*np.sin(np.deg2rad(two_theta)/2)/wave_length

def wave_vector_length_2_two_theta(wave_length,wave_vector_length):

    """ Convert wave_vector_length to two_theta

    Inputs
    wave_length : x-ray source wave length (float)
    wave_vector_length : length of the wave vector of the reciprocal space in rad/units of wave_length input  (float, list, np.array)

    Output
    two_theta : angle between the source and the detector in degree (float, list, np.array)
    """

    return 2*np.rad2deg(np.arcsin(wave_length*wave_vector_length/(4*np.pi)))

def theta_source_theta_detector_2_wave_vector_x_z(wave_length,theta_source,theta_detector):

    """ Convert pairs of theta source and theta detector to x and z components a wave vector

    Inputs
    wave_length : x-ray source wave length (float)
    theta_source : angle of the source relative to the x axis
    theta_detector : angle of the detector relative to the x axis

    Outputs
    wvx : x component of the wave vector
    wvz : z component of the wave vector

    notes : 
    theta_source + theta_detector = two_theta
    x and z are the axis in the same plane defined by the source, the detector and the sample center
    """

    alpha, beta = theta_source_theta_detector_2_alpha_beta(theta_source,theta_detector)

    return alpha_beta_2_wave_vector_x_z(wave_length,alpha, beta)

def theta_source_two_theta_2_wave_vector_x_z(wave_length,theta_source,two_theta):

    """ Convert pairs of theta source and theta detector to x and z components a wave vector

    Inputs
    wave_length : x-ray source wave length (float)
    theta_source : angle of the source relative to the x axis
    two_theta : angle between the source and the detector in degree

    Outputs
    wvx : x component of the wave vector
    wvz : z component of the wave vector

    notes : 
    theta_source + theta_detector = two_theta
    x and z are the axis in the same plane defined by the source, the detector and the sample center
    """

    alpha = (2*theta_source - two_theta)/2
    beta = two_theta/2

    return alpha_beta_2_wave_vector_x_z(wave_length,alpha, beta)

def alpha_beta_2_wave_vector_x_z(wave_length,alpha,beta):

    """ Convert pairs of theta source and theta detector to x and z components a wave vector

    Inputs
    wave_length : x-ray source wave length (float)
    alpha : angle between the bisectrix between source-sample and detector-sample and z axes
    beta : half the angle between source-sample and detector-sample (half of two theta)

    Outputs
    wvx : x component of the wave vector
    wvz : z component of the wave vector

    notes : 
    theta_source + theta_detector = two_theta
    x and z are the axis in the same plane defined by the source, the detector and the sample center
    """

    wave_vector_length = two_theta_2_wave_vector_length(wave_length,2*beta)

    wvx = wave_vector_length*np.sin(np.deg2rad(alpha))
    wvz = wave_vector_length*np.cos(np.deg2rad(alpha))

    return wvx, wvz

def theta_source_theta_detector_phi_2_wave_vector_x_y_z(wave_length,theta_source,theta_detector,phi):

    """ Convert trios of theta source, theta detector and phi to x, y and z components a wave vector

    Inputs
    wave_length : x-ray source wave length (float)
    theta_source : angle of the source relative to the xy plane
    theta_detector : angle of the detector relative to the xy plane
    phi : angle between the rotated x axis and the plane defined by the source, the detector and the sample center

    Outputs
    wvx : x component of the wave vector
    wvy : x component of the wave vector
    wvz : z component of the wave vector

    notes : 
    theta_source + theta_detector = two_theta
    theta_source_theta_detector_phi_2_wave_vector_x_y_z(...,0) = theta_source_theta_detector_2_wave_vector_x_z(...)
    """

    wvxy,wvz = theta_source_theta_detector_2_wave_vector_x_z(wave_length, theta_source, theta_detector)
    
    wvx = wvxy * np.cos(np.deg2rad(phi))
    wxy = wvxy * np.sin(np.deg2rad(phi))

    return wvx, wxy, wvz

def theta_source_theta_detector_2_alpha_beta(theta_source,theta_detector):

    """ Convert pairs of theta source and theta detector to alpha and beta angle

    Inputs
    theta_source : angle of the source relative to the x axis
    theta_detector : angle of the detector relative to the x axis

    Outputs
    alpha : angle between the bisectrix between source-sample and detector-sample and z axes
    beta : half the angle between source-sample and detector-sample (half of two theta)
    """

    beta  = (theta_source + theta_detector)/2
    alpha = (theta_source - theta_detector)/2

    return alpha, beta

def xi_zeta_2_alpha_phi(xi,zeta):

    # Tilt stage adjustmen to normal direction.

    xi_rad = np.deg2rad(xi)
    zeta_rad = np.deg2rad(zeta)

    alpha = np.rad2deg(np.arctan(np.sqrt(np.tan(xi_rad)**2 + np.tan(zeta_rad)**2)))
    phi = np.rad2deg(np.arctan2(np.tan(xi_rad),np.tan(zeta_rad)))

    return alpha, phi



