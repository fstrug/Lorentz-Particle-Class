# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 10:09:45 2021

@author: Frank
"""
############################Import Libraries
import numpy as np
import math
import time 
from helper import *


############################Create particle class object
class particle_lorentz:
    def __init__(self, four_vector, tag = None):
        self.four_vector = four_vector
        #Momentum Attributes
        self.Px = four_vector[1]
        self.Py = four_vector[2]
        self.Pz = four_vector[3]
        self.Pmag = np.sqrt(four_vector[1]**2 + four_vector[2]**2 + four_vector[3]**2)
        self.Pt = np.sqrt(four_vector[1]**2 + four_vector[2]**2)
        #Particle Energy
        self.E = four_vector[0]
        self.M = np.sqrt(self.E**2 - self.Pmag**2)
        if math.isnan(self.M):
            self.M = 0 
        #Particle Velocity (not defined for null particles with four vector [0,0,0,0])
        try:
            self.V = self.Pmag / self.E
        except:
            pass
        self.tag = tag
            
    #Combines two particle object into a singular particle
    #returns di-jet system particle object
    def combine_jets(self, particle):
        combined_jets = self.four_vector + particle.four_vector
        combined_jets_object = particle_lorentz(combined_jets)
        return(combined_jets_object)
   
############################Create particle decayer object
class particle_Decay:
    def __init__(self, mass_dictionary, decay_dictionary):
        self.mass_dictionary = mass_dictionary
        self.decay_dictionary = decay_dictionary
    
    def select_decay_branch(self, particle):
        try:
            decay_modes = self.decay_dictionary[particle.tag]
        except:
            decay_modes = None
            return(decay_modes)
        #choose a decay mode according to probability
        r = np.random.rand()
        x_prob_cdf = []
        for i in range(len(decay_modes)):
            x_prob_cdf.append(decay_modes[i][0]+np.sum(x_prob_cdf)) 
        index = next(x for x, val in enumerate(x_prob_cdf)
                                  if val > r)
        decay_mode = decay_modes[index]
        
        #Parse through and return final state particles
        particles_decay_branch = []
        for i in range(len(decay_mode)-1):
            particles_decay_branch.append(decay_mode[i+1])
        return(particles_decay_branch)
    
    def two_body_decay(self, particle, mode = None):
        #Select decay mode
        if mode == None:
            particles = self.select_decay_branch(particle)
        else: 
            particles = mode
        ##Calculate four vectors
        daughter_1_M = np.random.normal(mass_dictionary[particles[0]][0], mass_dictionary[particles[0]][1])
        daughter_2_M = np.random.normal(mass_dictionary[particles[1]][0], mass_dictionary[particles[1]][1])
        ###Create daughter1
        ###Construct daughter_1 four vector
        daughter_1_pmag = calc_mag_p1(particle.M, daughter_1_M, daughter_2_M)
        #Decay isotropically in parent frame
        phi = np.random.uniform(0, 2*np.pi)
        theta = np.arccos(np.random.uniform(-1, 1))
        #Find px, py, pz, and E
        daughter_1_px = daughter_1_pmag * np.sin(theta) * np.cos(phi)
        daughter_1_py = daughter_1_pmag * np.sin(theta) * np.sin(phi)
        daughter_1_pz = daughter_1_pmag * np.cos(theta)
        daughter_1_E = np.sqrt(daughter_1_M**2 + daughter_1_pmag**2)
            
        #create four vector
        daughter_1_four_vector = np.array([daughter_1_E, daughter_1_px, daughter_1_py, daughter_1_pz])
        
        ###Construct daughter_2 four vector
        daughter_2_pmag = daughter_1_pmag
        daughter_2_E = np.sqrt(daughter_2_M**2 + daughter_2_pmag**2)
        daughter_2_four_vector = np.array([daughter_2_E, -daughter_1_px, -daughter_1_py, -daughter_1_pz])
        ###Check for boost
        if particle.V != 0:
            parent_rest_four_vector, boost_particle_parent_rest = boost_4_vector_rest(particle.four_vector)
            boost_lab = np.linalg.inv(boost_particle_parent_rest)
            #Boost daughters into lab frame
            daughter_1_lab_four_vector = np.matmul(boost_lab, daughter_1_four_vector)
            daughter_1_particle = particle_lorentz(daughter_1_lab_four_vector, tag = particles[0])
                
            daughter_2_lab_four_vector = np.matmul(boost_lab, daughter_2_four_vector)
            daughter_2_particle = particle_lorentz(daughter_2_lab_four_vector, tag = particles[1])
           
        #if parent is at rest, no boosting is neccessary
        else:
            daughter_1_particle = particle_lorentz(daughter_1_four_vector, tag = particles[0])
            daughter_2_particle = particle_lorentz(daughter_2_four_vector, tag = particles[1])
            
        return(daughter_1_particle, daughter_2_particle) 
    
   
############################Create event object
class event_lorentz:
    def __init__(self, particles):
        self.particles = particles
        
        self.tags = []
        for particle in particles:
            self.tags.append(particle.tag)
        self.n_W = self.tags.count("W")
        self.n_top = self.tags.count("top")
        self.n_bottom = self.tags.count("bottom")
        self.n_electron = self.tags.count("electron")
        self.n_muon = self.tags.count("muon")
        self.n_tau = self.tags.count("tau")
        self.n_neutrino = self.tags.count("neutrino")
    
class events_lorentz:
    def __init__(self, events):
        self.events = events
        
    def cut_nparticles(self, particle, n):
        events_pass_cut = []
        for event in self.events:
            if event.tags.count(particle) >= n:
                events_pass_cut.append(event)
        return(events_lorentz(events_pass_cut))
      
############################Standard Model Parameters

decay_dictionary = {
    "top": [[1.0, "W", "bottom"],
            [0.0, "fake", "fake", "fake"]],
    "W"  : [[0.33, "electron", "neutrino"],
            [0.33, "muon", "neutrino"],
            [0.34, "tau", "neutrino"]]
    }
mass_dictionary = {
    "top"   : [172.76, 3],
    "W"     : [80.369, 2],
    "bottom": [4.18, 0.1],
    "electron": [0.000511,0.00005],
    "muon" : [0.105, 0.01],
    "tau" : [1.776,0.1],
    "neutrino" : [0,0]
    }
############################Simulation Example/Debug
if True:
    number_of_simulations = 5
    events = []
    top = particle_lorentz([172.76,0,0,0], "top")
    decayer = particle_Decay(mass_dictionary, decay_dictionary)
    for i in range(number_of_simulations):
        particles = [top]
        decay_possible = True
        while decay_possible == True:
            #Check that any particle can decay
            decay_possible = False
            for particle in particles:
                if decayer.select_decay_branch(particle) != None:
                    decay_possible = True
            for particle in particles:
                if decayer.select_decay_branch(particle) == None:
                    continue
                particle_1, particle_2 = decayer.two_body_decay(particle)
                particles.append(particle_1)
                particles.append(particle_2)
                particles.remove(particle)
        event = event_lorentz(particles)
        events.append(event)
    events = events_lorentz(events)
    
# for i, event in enumerate(events_lorentz.events):
#     print("Event #", i)
#     mother_particle = particle_lorentz([0,0,0,0])
#     for particle in event.particles:
#         print(particle.tag, particle.M)
#         mother_particle = mother_particle.combine_jets(particle)
#     print(mother_particle.M)
    
events_muons = events.cut_nparticles("muon", 1)
for i, event in enumerate(events_muons.events):
    print("Event #", i)
    mother_particle = particle_lorentz([0,0,0,0])
    for particle in event.particles:
        print(particle.tag, particle.M)
        mother_particle = mother_particle.combine_jets(particle)
    print(mother_particle.M)
    

        
                
                
            














