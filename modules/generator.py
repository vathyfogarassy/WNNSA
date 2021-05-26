# -*- coding: utf-8 -*-
"""
Created on Thu Feb 01 00:39:44 2018

@author: punck
"""

import numpy as np
import scipy.stats as ss

import pandas as pd

class Generator:
    """A random dataset generator class"""
    
    def Binomial(self, n, p, size):
        """
        Dataset of random binomial variables with probability of p and size of size
        n = number of trials
        p = probability of success
        size = output shape
        """
        
        if not isinstance(size, tuple):
            size = (size, 1)
        
        y = np.random.binomial(n, p, size)
        columns = []
        length = size[1]
        for i in range(length):
            columns.append('x{}'.format(i+1))            
        df = pd.DataFrame(y, columns=columns)
        return df
        
    def Bernoulli(self, p, size):
        """
        Dataset of Bernoulli random variables with probability of p and size of size
        n = number of trials
        p = probability of success
        size = output shape
        """
    
        if not isinstance(size, tuple):
            size = (size, 1)        
        
        y = np.random.binomial(1, p, size)
        columns = []
        length = size[1]
        for i in range(length):
            columns.append('x{}'.format(i+1))            
        df = pd.DataFrame(y, columns=columns)
        return df
    
    def Normal(self, mu, sigma, size):
        """
        Dataset of variables of normal distribution with probability of p and size of size
        n = number of trials
        p = probability of success
        size = output shape
        """
        
        if not isinstance(size, tuple):
            size = (size, 1)        
        
        y = np.random.normal(mu, sigma, size)     
        columns = []
        length = size[1]
        for i in range(length):
            columns.append('x{}'.format(i+1))        
        df = pd.DataFrame(y, columns=columns)
        return df
    
    def Uniform(self,low, high, size):
        """
        [low, high)
        """
        if not isinstance(size, tuple):
            size = (size, 1)        
        
        y = np.random.rand(*size) * (high - low) + low
        columns = []
        length = size[1]
        for i in range(length):
            columns.append('x{}'.format(i+1))        
        df = pd.DataFrame(y, columns=columns)
        return df
    
    def DiscreteUniform(self, low, high, size):
        """
        [low, high)
        """
        if not isinstance(size, tuple):
            size = (size, 1)        
        
        y = np.random.randint(low, high, size=size)     
        columns = []
        length = size[1]
        for i in range(length):
            columns.append('x{}'.format(i+1))        
        df = pd.DataFrame(y, columns=columns)
        return df
        
    def DiscreteNormal(self, low, high, size, scale=1):
        """
        [low, high)
        """
        if not isinstance(size, tuple):
            size = (size, 1)        
        
        x = np.arange(low, high + 1)
        xU, xL = x + 0.5, x - 0.5
        loc=( ( low + high - 1 ) / 2 )
        prob = ss.norm.cdf(xU, loc=loc, scale=scale) - ss.norm.cdf(xL, loc=loc, scale=scale)
        prob = prob / prob.sum()
        y = np.random.choice(x, size=size, p=prob)
        #y = ss.truncnorm(a=low/scale, b=high/scale, scale=scale).rvs(size=size)
        
        columns = []
        length = size[1]
        for i in range(length):
            columns.append('x{}'.format(i+1))        
        df = pd.DataFrame(y, columns=columns)
        
        #plt.hist(df['x1'], bins = len(x))        
        #plt.show()
        
        return df
        
    def DiscreteNormal2(self, mu, sigma, size):
        """
        [low, high)
        """
        if not isinstance(size, tuple):
            size = (size, 1)        
        
        x = np.arange(mu - sigma, mu + sigma + 1)
        xU, xL = x + 0.5, x - 0.5
        prob = ss.norm.cdf(xU, loc=mu, scale=sigma) - ss.norm.cdf(xL, loc=mu, scale=sigma)
        prob = prob / prob.sum()
        y = np.random.choice(x, size=size, p=prob)
        #y = ss.truncnorm(a=low/scale, b=high/scale, scale=scale).rvs(size=size)
        
        columns = []
        length = size[1]
        for i in range(length):
            columns.append('x{}'.format(i+1))        
        df = pd.DataFrame(y, columns=columns)
        
        #plt.hist(df['x1'], bins = len(x))        
        #plt.show()
        
        return df
        
    def Concat(self, *args):
        if args:
            dfs = [df for df in args]
            sizes = [df.shape[0] for df in dfs]
            
            if not all(x == sizes[0] for x in sizes):
                print('Sizes are not identical!')
                return None
                
            df = pd.concat(dfs, axis=1)
            columns = ['x{}'.format(i+1) for i in range(len(df.columns))]
            df.columns = columns        
            return df
        
        print('Missing arguments!')
        return None    

def localize_floats(row):
    return [
        str(el).replace('.', ',') if isinstance(el, float) else el 
        for el in row
    ]