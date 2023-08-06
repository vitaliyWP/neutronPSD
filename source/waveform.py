#! /usr/bin python

import numpy as np

class waveformsFile:
    

class Waveform:
    def __init__(self, data):
        self.data = data
    
    def max(self):
        return np.max(self.data)
    
    def locmax(self):
        return np.argmax(self.data)
    
    def get_integral(self):
        return np.sum(self.data[ np.where(self.data > 0) ])
    
    def get_fps(self):
        gradself.data
