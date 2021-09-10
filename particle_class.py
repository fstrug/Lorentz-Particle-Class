# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 10:09:45 2021

@author: Frank
"""
############################Import Libraries
import numpy as np
import time 
from helper import *

############################Create particle class object

class particle_lorentz:
    def __init__(self, four_vector, tag, decay_dictionary):
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
        #Particle Velocity
        self.V = self.Pmag / self.E
        self.tag = tag
        self.decay_dictionary = decay_dictionary
        
    #Takes list of avaiable mass states and see which pairs are accessible
    def check_decay(self, mass_states):
        available_states = []
        #we want to check every combination of candidates without repeating
        for i in range(len(mass_states)):
            for j in range(len(mass_states[i:])):
                #check if final state is accessible
                if self.M >= mass_states[i+j] + mass_states[i+j]:
                    #if mass state is accesible, add it to list of available
                    #states
                    available_states.append([mass_states[i], mass_states[i+j]])
       
        return(available_states)
    
    def select_decay_branch(self):
        decay_modes = self.decay_dictionary[self.tag]
        #choose a decay mode according to probability
        r = np.random.rand()
        r = .31
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
    
    def two_body_decay(self, mode = None):
        #Select decay mode
        if mode == None:
            particles = select_decay_branch()
        else: 
            particles = mode
        
        
        
    
    #Combines two particle object into a singular particle
    #returns di-jet system particle object
    def combine_jets(self, particle):
        combined_jets = self.four_vector + particle.four_vector
        combined_jets_object = particle_lorentz(combined_jets)
        return(combined_jets_object)
    
############################Standard Model Parameters

decay_dictionary = {
    "top": [(1.0, "W", "bottom"),
            (0.0, "fake", "fake", "fake")]
    }
mass_dictionary = {
    "top"   : 173.29,
    "W"     : 80.369,
    "bottom": 1.14
    }
    
top = particle_lorentz([10,0,0,0], "top")














