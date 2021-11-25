import numpy as np
import scipy.spatial.distance as sc
import sys

import pandas as pd

import wnnem

def normalize(sample, population, attributes, w, limits=None):
    if limits is None:
        limits = np.concatenate((sample[attributes], population[attributes]))
    
    sample_norm = sample[attributes].copy()
    sample_norm = sample_norm * w
    population_norm = population[attributes].copy()
    population_norm = population_norm * w
    
    sample_norm = (sample_norm - limits.min(axis=0)) / (limits.max(axis=0) - limits.min(axis=0))
    population_norm = (population_norm - limits.min(axis=0)) / (limits.max(axis=0) - limits.min(axis=0))    
    
    return sample_norm.as_matrix(), population_norm.as_matrix()

def calc_error(row, result_indices, distance_matrix):
    current = distance_matrix[row,:].argsort()[:result_indices[row]+2][-2]
    alter = distance_matrix[row,:].argsort()[:result_indices[row]+2][-1]
    return (distance_matrix[row, alter] - distance_matrix[row, current])# / distance_matrix[row, current]

def error_func(rows, result_indices, distance_matrix):
   "Calculates the error and returns the optimal index"      
   #error = [lambda row: calc_error(row, result_indices, distance_matrix), rows]
   #vmax = np.max(error)
   #indices = [i for i, x in enumerate(error) if x == vmax]
   #return np.random.choice(indices)
   return np.argmax([calc_error(row, result_indices, distance_matrix) for row in rows])

def  update_result(row, vector, indices, distance_matrix):  
    offset = 0
    if row.size == 1:
        offset = 1
    vector[row] = distance_matrix[row,:].argsort()[:indices[row]+2][-2+offset]
    return

def update_result_vector( rows, vector, indices, distance_matrix):
    if isinstance(rows, (int, np.integer)):
        update_result(rows, vector, indices, distance_matrix)
    else:    
        for row in rows:
            update_result(row, vector, indices, distance_matrix)

def _unique_dist( a, b ):   
    w = len(a)
    abs_diff = np.absolute(a - b)
    d = np.sum(abs_diff) / w
    
    return d

def energy( result, e_sam, e_pop, P):   
    res = P.iloc[result]
    res = np.asarray(res)
    en = [_unique_dist(e_sam[i], e_pop[result[i]]) for i in range(len(result))]
    #print(np.sum(en), np.sum(en) + (e_sam.shape[0] - len(result)), e_sam.shape[0], len(result))
    return np.sum(en) + (e_sam.shape[0] - len(result))
    #return np.sum(distance_matrix[:,result])

def WNNSA(distance_matrix, k, sT, e_sam, e_pop, P, max_it):   
    
    with open('simann_LOG_ELES.txt', 'w') as f:
    
        print_info = True
        
        n = k
        full_count = distance_matrix.shape[0]
        
        distance_matrix[distance_matrix == 0] = np.nextafter(np.float64(0), np.float64(1))
        
        reduced_d_val = np.sort(distance_matrix, axis=1)[:,0:n]
        reduced_d_idx = np.argsort(distance_matrix, axis=1)[:,0:n]
        if print_info:
            print(reduced_d_val.dtype, file=f)
          
        be = sys.float_info.max    
        Tmax = sT
        for T in range(Tmax+1, 1, -1):
            
            red = (reduced_d_val**((Tmax-T)*1.0/(Tmax*1.0)))
            if print_info:
                print('### {} ###'.format(T), file=f)
                print('### red ###', file=f)   
                print(red, file=f)
            
            red[red < sys.float_info.min] = sys.float_info.min
            
            p = 1.0/red
            
            if print_info:
                print('### {p} ###', file=f)
                print(p, file=f)
            
            sum_p = np.sum(p, axis=1)
            tmp = sum_p
            sum_p[sum_p > sys.float_info.max] = sys.float_info.max
            tmp = tmp - sum_p
            
            if print_info:
                print('### {sum_p} ###', file=f)
                print(sum_p, file=f)
                
                print('### {diff_p} ###', file=f)
                print(tmp, file=f)
                
            
            p = p / sum_p[:, None]
            
            sp = np.sum(p, axis=1)
            p = p / sp[:, None]
            
            if print_info:
                print('### {p} ###', file=f)
                print(p, file=f)
                
                print('### {sum_p} ###', file=f)
                print(np.sum(p, axis=1), file=f)
                print('~~~~~~~~~~', file=f)
            
       # try:
            tmp_indices = [np.random.choice(range(n), p=p[i,:]) for i in range(full_count)]
            tmp_indices = np.asarray(tmp_indices)
            tmp_vector = [reduced_d_idx[i,tmp_indices[i]] for i in range(len(tmp_indices))]
            tmp_vector = np.asarray(tmp_vector)
        
            res_count = np.unique(tmp_vector).shape[0]
            iteration = 1
            while iteration <= max_it and res_count != full_count:
                vals, inverse, count = np.unique(tmp_vector, return_inverse=True, return_counts=True)
                idx_vals_repeated = np.where(count > 1)[0]
                vals_repeated = vals[idx_vals_repeated]
                
                for val in vals_repeated:
                    conflicting_indices =np.where(tmp_vector == val)
                    conflicting_indices = [item for sublist in conflicting_indices for item in sublist]
                    conflicting_indices.pop(error_func(conflicting_indices, tmp_indices, distance_matrix))
                    if len(conflicting_indices) == 1:
                        conflicting_indices = conflicting_indices[0]
                    tmp_indices[conflicting_indices] = tmp_indices[conflicting_indices] + 1
                    tmp_indices = tmp_indices % n
                    update_result_vector(conflicting_indices, tmp_vector, tmp_indices, distance_matrix)
                    
                res_count = np.unique(tmp_vector).shape[0]
                
                iteration += 1    
                
            te = energy(tmp_vector, e_sam, e_pop, P)
            if te < be:
                result_vector = tmp_vector
                be = te     
             
    return result_vector, iteration - 1

def uni_arr(arr):
    
    for i in range(len(arr)-1):
        for j in range(i+1, len(arr)):
            if arr[j] == arr[i]:
                arr[j] = np.nan
    
    return arr

def match(_to, _from, ovar, w, **kwargs):
    """
    match(_to, _from, ovar, w, **kwargs)

    Pairs individuals of _to and _from.
    
    Parameters
    ----------
    _to : Pandas.DataFrame
        Case group.
    _from : Pandas.DataFrame
        Population.    
    ovar : List
        Observed variables.
    w : List
        Weights.
    **kwargs: {'pair_name'}
        index     : string (default 'index')
        Name of the index column.
        pair_name : string (default 'pair')
        Name of the pair column.
           
    Returns
    -------
    Pandas.DataFrame
        Control group.
    """
    
    pd.options.mode.chained_assignment = None
    
    pair_name = kwargs['pair_name'] if 'pair_name' in kwargs else 'pair'
    _id = kwargs['index'] if 'index' in kwargs else 'index'
    
    if 'dmat' in kwargs:
        d = kwargs['dmat']
    else:
        norm_to, norm_from = wnnem._normalize(_to, _from, ovar, w)
        d = sc.cdist(norm_to, norm_from, _unique_dist)
    
    k = int(est_k(d) * 1.5)
    print('Reduced environment: {}'.format(k))
    sT = 200
    max_it = 500
    
    res, it = WNNSA(d, k, sT, norm_to, norm_from, _from, max_it)
      
    res = uni_arr(res.astype('float64'))
    
    otpt = _from.reset_index().rename(columns={_from.index.name:_id})
      
    _to['tmp'] = res
    
    res = res[~np.isnan(res)]
    control = otpt.loc[res]
      
    for index, row in _to.iterrows():
        try:
            _to.loc[index, pair_name] = control.loc[_to.loc[index, 'tmp'], _id]
        except KeyError:
            _to.loc[index, pair_name] = None
    
    del _to['tmp']
    
    return control[[_id] + ovar + ['treated', 'ps']].dropna().set_index(_id)

def est_k(dist):
    k=1
    utkozes = 1
    while utkozes > 0:
        dist_rank = np.zeros(dist.shape)
        index=np.arange(0, len(dist), 1, dtype=int)
        for idx in index:
            temp = dist[idx].argsort()
            ranks = np.empty(len(dist[idx]), int)
            ranks[temp] = np.arange(len(dist[idx]))
            dist_rank[idx] = np.transpose(ranks)
        
        dist_rank = dist_rank + 1
        messziek=np.argwhere(dist_rank > k)
        dist_rank_kozeli=dist_rank
        dist_rank_kozeli[messziek[:, 0], messziek[:, 1]] = 0
        dist_rank_bool=dist_rank_kozeli
        hol = np.argwhere(dist_rank_kozeli > 0)
        dist_rank_bool[hol[:, 0], hol[:, 1]] = 1
        
        P_igeny = dist_rank_bool.sum(axis=0)*(1.0/k)
        dist_igeny = P_igeny * dist_rank_bool
        dist_igeny[messziek[:, 0], messziek[:, 1]] = 999
        A_minimum = dist_igeny.min(axis=1)
        problemas = np.argwhere (P_igeny > 1)
        problemas=np.transpose(problemas)
        problemas=problemas[0, :]
        
        utkozes = 0
        for prob in problemas:
            osszeg = sum(A_minimum[np.argwhere(dist_rank_bool[:, prob] == 1)])
            elemszam = len(np.argwhere(dist_rank_bool[:, prob] == 1))
            if osszeg/elemszam > 1:
                utkozes = utkozes +1
        k = k+1
        
    return k
