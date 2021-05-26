__author__ = "Szabolcs Szekér  <szeker@dcs.uni-pannon.hu>"
__copyright__ = "Copyright (C) 2019 Szabolcs Szekér"
__license__ = "Public Domain"
__version__ = "0.5"


# %%

import bisect
import numpy as np
import pandas as pd
from tableone import TableOne
    
# %%
  
def _measure_dissimilarity(method, sample, population, ovar, attribute_type, pair_name, w = None):
    n = sample.shape[0]
    N = len(ovar)
    m = population.shape[0]
    
    d = 0
    #d = abs(n-r) * N
    
    if method == 'NNI':
        for index, row in sample.iterrows():
            try:
                sample_element = row[ovar]
                control_element = population.loc[sample.loc[index, pair_name]][ovar]
                
                pd = 0.0        
                
                for attribute in range(N):
                    if attribute_type[attribute] == 'c':
                        s_c1 = sample_element[attribute]
                        c_c1 = control_element[attribute]
                        c1_sorted = np.sort(population, axis=0)[:, attribute]
                        c1_sorted_list = list(c1_sorted)
                        bisect.insort(c1_sorted_list, s_c1)
                        c1_sorted_list_desc = list(reversed(c1_sorted_list))
                        
                        if s_c1 > c_c1:
                            s_c1_idx = c1_sorted_list.index(s_c1)
                            c_c1_idx = c1_sorted_list_desc.index(c_c1)
                            c_c1_idx = m - c_c1_idx
                        else:
                            s_c1_idx = c1_sorted_list_desc.index(s_c1)
                            s_c1_idx = m - s_c1_idx
                            c_c1_idx = c1_sorted_list.index(c_c1)
                        
                        pd = pd + (0.0 if abs(s_c1_idx-c_c1_idx) == 1 else 1.0)
                    elif attribute_type[attribute] == 'b' or attribute_type[attribute] == 'o' or attribute_type[attribute] == 'n':
                        pd = pd + (0.0 if sample_element[attribute] == control_element[attribute] else 1.0)
                    else:
                        print('Unknown attribute type!')
                        return np.nan
                            
                d = d + pd
            except TypeError:
                d = d + N
    elif method == 'LDI':    
        
        limits = np.concatenate((sample[ovar], population[ovar]))
    
        sample_norm = sample[ovar].copy()
        sample_norm = sample_norm * w
        population_norm = population[ovar].copy()
        population_norm = population_norm * w
        
        sample_norm = (sample_norm - limits.min(axis=0)) / (limits.max(axis=0) - limits.min(axis=0))
        population_norm = (population_norm - limits.min(axis=0)) / (limits.max(axis=0) - limits.min(axis=0))
        
        sample_norm[pair_name] = sample[pair_name]
        
        for index, row in sample_norm.iterrows():
            try:
                sample_element = row[ovar]
                control_element = population_norm.loc[sample_norm.loc[index, pair_name]][ovar]
                
                pd = 0.0        
                
                for attribute in range(N):
                    if attribute_type[attribute] == 'b' or attribute_type[attribute] == 'o':
                        pd = pd + abs(sample_element[attribute] - control_element[attribute])
                    elif attribute_type[attribute] == 'c':
                        s_c1 = sample_element[attribute]
                        c1_sorted = np.sort(population, axis=0)[:,attribute]
                        c1_sorted_list = list(c1_sorted)
                        bisect.insort(c1_sorted_list, s_c1)
                        s_c1_idx = c1_sorted_list.index(s_c1)
                        if s_c1_idx == 0:
                            d_nn = abs(sample_element[attribute]-c1_sorted_list[s_c1_idx+1])
                        elif s_c1_idx == m:
                            d_nn = abs(sample_element[attribute]-c1_sorted_list[s_c1_idx-1])
                        else:
                            d_nn = min(abs(sample_element[attribute]-c1_sorted_list[s_c1_idx-1]), abs(sample_element[attribute]-c1_sorted_list[s_c1_idx+1]))
                        d_sel = abs(sample_element[attribute]-control_element[attribute])
                        
                        pd = pd + (0.0 if d_sel == 0 else (1 - (d_nn*1.0)/(d_sel*1.0)))
                    elif attribute_type[attribute] == 'n':
                        pd = pd + (0.0 if sample_element[attribute] == control_element[attribute] else 1.0)
                    else:
                        print('Unknown attribute type!')
                        return np.nan
                            
                d = d + pd
            except TypeError:
                d = d + N
    elif method == 'GDI':
        limits = np.concatenate((sample[ovar], population[ovar]))
    
        sample_norm = sample[ovar].copy()
        sample_norm = sample_norm * w
        population_norm = population[ovar].copy()
        population_norm = population_norm * w
        
        sample_norm = (sample_norm - limits.min(axis=0)) / (limits.max(axis=0) - limits.min(axis=0))
        population_norm = (population_norm - limits.min(axis=0)) / (limits.max(axis=0) - limits.min(axis=0))
        
        sample_norm[pair_name] = sample[pair_name]
        
        i = 0
        for index, row in sample_norm.iterrows():
            try:
                sample_element = row[ovar]
                control_element = population_norm.loc[sample_norm.loc[index, pair_name]][ovar]
                
                pd = 0.0        
                for attribute in range(N): pd = pd + abs(sample_element[attribute] - control_element[attribute])
                d = d + pd
                i += 1
            except TypeError:
                d = d + N
        #print('pairs: {}'.format(i))
    else:
        print('Unknown method!')
        return np.nan
    
    d = d / (n * N * 1.0)
        
    return d


# %%

def NNI(case, control, ovar, attribute_types, **kwargs):
    pair_name = kwargs['pair_name'] if 'pair_name' in kwargs else 'pair'
    return _measure_dissimilarity('NNI', case, control, ovar, attribute_types, pair_name)


# %%
    
def LDI(case, control, ovar, attribute_types, **kwargs):
    pair_name = kwargs['pair_name'] if 'pair_name' in kwargs else 'pair'
    w = kwargs['weights'] if 'weights' in kwargs else None
    return _measure_dissimilarity('LDI', case, control, ovar, attribute_types, pair_name, w)


# %%
    
def GDI(case, control, ovar, attribute_types, **kwargs):
    pair_name = kwargs['pair_name'] if 'pair_name' in kwargs else 'pair'
    w = kwargs['weights'] if 'weights' in kwargs else None
    return _measure_dissimilarity('GDI', case, control, ovar, attribute_types, pair_name, w)


# %%
 
def calculate_SMD(df, treatment_column, ovar):
    return df[ovar].apply(lambda x: SMD(x, df[treatment_column])).round(4)


# %%    

def SMD(feature, treatment):
    """Calculate the standard mean difference (SMD) of a feature between the
    treatment and control groups.

    The definition is available at
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3144483/#s11title

    Args:
        feature (pandas.Series): a column of a feature to calculate SMD for
        treatment (pandas.Series): a column that indicate whether a row is in
                                   the treatment group or not

    Returns:
        (float): The SMD of the feature
    """
    
    t = feature[treatment == 1]
    c = feature[treatment == 0]

    return (t.mean() - c.mean()) / np.sqrt(.5 * (t.var() + c.var()))


# %%
   
def DDI(sample, control, fields, ranges, bins):   
    ns = sample.shape[0];
    nc = control.shape[0]; 
    N = len(fields)
    d = np.float64(0);

    for i in range(len(fields)):
        [s_hist, s_bin_edges] = np.histogram(sample[fields[i]], bins=bins[i], range=ranges[i])
        [c_hist, c_bin_edges] = np.histogram(control[fields[i]], bins=bins[i], range=ranges[i])
        d = d + np.sum(np.absolute(s_hist - c_hist));

    d = d / N;
    d = (d + np.abs(ns - nc)) / ( 2 * ns);
    
    return d


# %%
    
def calculate_individual_balance(df1, df2, columns, categorical, groupby, nonnormal, pval = True):
    df = pd.concat([df1 , df2])

    #print(sorted(df.groupby('treated').groups.keys()))

    bal_tab = TableOne(df, columns, categorical, groupby, nonnormal, pval=pval)
    return bal_tab.tableone

