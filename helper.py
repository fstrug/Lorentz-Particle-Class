# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 10:05:00 2021

@author: Frank
"""
import numpy as np

############################Define Helper functions
#Calculates pmag of a four vector
def calc_pmag_4_vector(input_four_vector):
    pmag = np.sqrt(input_four_vector[1]**2 + input_four_vector[2]**2 + input_four_vector[3]**2)
    return(pmag)

#Calculates pt of a four vector, transverse to the beampipe (z_hat)
def calc_pt_4_vector(input_four_vector):
    pt = np.sqrt(input_four_vector[1]**2 + input_four_vector[2]**2)
    return(pt)

#When decaying in the rest frame of some parent particle, the magnitude of the
#momentom for either of it's daughter's can be derived from momentum conservation
def calc_mag_p1(M, m1, m2):
    num = np.sqrt((M**2 - (m1 + m2)**2)*(M**2 - (m1 - m2)**2))
    denom = 2*M
    return(num/denom)

#Produces the boost matrix along the z-axis for some beta
def boost_matrix_z(beta):
    gamma = 1/np.sqrt(1-beta**2)
    boost = [[gamma, 0, 0, gamma*beta],
             [0,1,0,0],
             [0,0,1,0],
             [gamma*beta, 0, 0, gamma]]
    return(boost)


#Produces the rotation matrix according to the Rodrigues Rotation Algorithm.
#This matrix rotates a vector counter-clockwise by an angle theta about some
#rotation axis.
def rotation_matrix(rotation_axis, theta):
    norm_axis = rotation_axis / np.linalg.norm(rotation_axis)
    #if rotation_axis = [0,0,0], return identity matrix
    if rotation_axis[0]==0 and rotation_axis[1]==0 and rotation_axis[2]==0:
        return(np.identity)
    K = [[0, -norm_axis[2], norm_axis[1]],
         [norm_axis[2], 0, -norm_axis[0]],
         [-norm_axis[1], norm_axis[0], 0]]
    #Convert to numpy array for speed and numpy function compatibility
    K = np.array(K)
    K_squared = np.matmul(K, K)
    K = np.sin(theta) * K    
    K_squared = np.array(K_squared)
    K_squared = (1-np.cos(theta))*K_squared
    R = np.identity(3) + K + K_squared
    return(R)

#Rotates a 3-vector about some axis by an angle theta
def rotate_vector(input_vector, rotation_axis, theta):
    R = rotation_matrix(rotation_axis, theta)
    rotated_vector = np.matmul(R, input_vector)
    return(rotated_vector)

#Finds the matrix for some vector which brings us into it's rest frame.
#Takes particle four vector in lab frame and returns the rest-frame four vector 
#of the particle as well as the matrix that rotates and boosts into its
#rest frame. This matrix can be inversed to bring particles decayed in it's
#rest frame to the lab frame
def boost_4_vector_rest(vector):
    #Find the rotation matrix to change coordinate systems
    z_hat = [0,0,1]
    norm_4_vector = vector[1:] / np.linalg.norm(vector[1:])
    rotation_axis = np.cross(z_hat, norm_4_vector)

    #When boosting into rest frame, if the 4_vector points along z
    #the cross product will give [0,0,0]: we must control for this
    if rotation_axis[0]==0 and rotation_axis[1]==0 and rotation_axis[2]==0:
        rotation_axis = [0,0,1]
        theta = 0
    else:
        theta = -np.arccos(np.dot(z_hat, norm_4_vector))    
    R = rotation_matrix(rotation_axis, theta)
    #Bring R matrix to 4x4 space
    R = np.pad(R, ((1, 0),(1,0)), 'constant')
    R[0,0] = 1
    
    #Create lorentz matrix
    beta = - calc_pmag_4_vector(vector) / vector[0]
    boost_matrix = boost_matrix_z(beta)
    #Rotate lorentz matrix into vector's coordinate system
    boost_matrix_rotated_prime = np.matmul(boost_matrix, R)
    vector_rest_frame = np.matmul(boost_matrix_rotated_prime, vector)
    return(vector_rest_frame, boost_matrix_rotated_prime)   


#Calculates the angles of the four vector
#Polar angle is defined with respect to the z-axis, and the azimuthal angle
#is defined in the xy plane.
def calc_theta_phi(vector):
    theta = np.arccos(vector[3]/calc_pmag_4_vector(vector))
    cos_theta = vector[3]/calc_pmag_4_vector(vector)
    phi = np.arctan2(vector[1],vector[2])
    return(theta, phi, cos_theta)
    
###Check entire state of particles for available decay chains
def can_state_decay(state, mass_states):
    all_available_states = []
    for particle in state:
        particle_decay_states = particle.check_decay(mass_states)
        if len(particle_decay_states) > 0:
            all_available_states.append(particle_decay_states)
    
    if len(all_available_states) == 0:
        return False
    else:
        return True
    
#Calculates the invariant mass of a jet
def calc_inv_mass_jet(jet):
    #Some of the neutrino jets may return a negative under the square root
    #due to rounding errors.
    try: 
        inv_mass = np.sqrt(jet[0]**2 - (jet[1]**2 + jet[2]**2 + jet[3]**2))
    #In this event return that the invariant mass is 0
    except:
        return(0)
    return(inv_mass)

#calculates invariant mass of some dijet system
#takes 4-vectors of jets as a list for inputs
def calc_inv_mass_2jets(jet1, jet2):
    vector_dijet = np.add(jet1, jet2)
    inv_mass = np.sqrt(vector_dijet[0]**2 - (vector_dijet[1]**2 + vector_dijet[2]**2 + vector_dijet[3]**2))
    return(inv_mass)