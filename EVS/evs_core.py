# -*- coding: utf-8 -*-
"""
Created on Tue May 25 10:46:58 2021

@author: Thomas Maynadié
"""

from EVS import evs_cosmology
from EVS import evs_wrapper
from EVS import evs_threading
from EVS import evs_computation

import threading

import numpy as np
import os
from datetime import datetime

from enum import Enum

EVS_ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
EVS_RESULTS_PATH = EVS_ROOT_DIR + "\\Results"


    
class EVSErrCode(Enum):
    EVS_ERROR = 1
    EVS_STATUS = 0

def debug(msg, errorCode=EVSErrCode.EVS_STATUS):
    HEADER = datetime.now().strftime("%H:%M:%S.%f")[:-3] + " >> "
    
    if errorCode == EVSErrCode.EVS_ERROR: 
        ERROR_MSG = "[EVS ERROR] "        
    if errorCode == EVSErrCode.EVS_STATUS: 
        ERROR_MSG = "[EVS STATUS] "
    
    print(HEADER + ERROR_MSG + msg)
    

class EVS(object):
    def __init__(self):
        self.cosmo = evs_cosmology.EVSCosmology(self)
        self.wrapper = evs_wrapper.EVSWrapper()
        self.computationHelper = evs_computation.EVSComputation()
        
    """ Performs EVS calculations""" 
    def evs_calculation(self, NM_MAX):
        debug("Stating computation : ", EVSErrCode.EVS_STATUS)
        
        self.cosmo.updateParameters()
        
        debug("Parameters : " + str(self.cosmo.getParameters()), EVSErrCode.EVS_STATUS)
        
        self.wrapper.setup_evs_computation()

        results = self.computationHelper.compute(self, NM_MAX)
        
        debug("Computation done : ", EVSErrCode.EVS_STATUS)
        
        return results
    
    