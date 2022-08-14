# -*- coding: utf-8 -*-
"""
Created on Tue May 17 09:07:20 2022

@author: nkp68
"""

#keep stack question toy data here

#create random numbers array and retest
import pandas as pd
import numpy as np
trainset = pd.DataFrame(np.random.randint(0,100,size=(50, 50)))
testset = pd.DataFrame(np.random.randint(0,100,size=(50, 50)))
Samples=['Normal', 'Disease']
trainset['Sample'] = [i for i in Samples for _ in range(25)]
testset['Sample'] = [i for i in Samples for _ in range(25)]
trainclasses, testclasses = multiclasstest.mulitanalysis(trainset, testset, 5)