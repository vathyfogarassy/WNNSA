__author__ = "Szabolcs Szekér <szeker@dcs.uni-pannon.hu>"
__copyright__ = "Copyright (C) 2019 Szabolcs Szekér"
__license__ = "Public Domain"
__version__ = "0.5"


# %%

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
import sklearn.metrics.pairwise
import random

# %%

def calculate_ps(df, dependent_variable, ovar=0, C=1.0, **kwargs):
   
    header = list(df)[:-2] if ovar == 0 else ovar

    y = df[dependent_variable]
    X = df[header]

    propensity = LogisticRegression(
        C=C
        #, solver='lbfgs'
        , solver='newton-cg'
        , max_iter=1000
        #, class_weight='balanced'
        #, fit_intercept=False
        )
        
    propensity = propensity.fit(X, y)

    ps_score = propensity.predict_proba(X)[:, 1]
    
    col_name = kwargs['col_name'] if 'col_name' in kwargs else 'ps'
    df[col_name] = ps_score

    orvalues = list(np.exp(propensity.coef_).flatten())

    return df, orvalues, list(propensity.coef_.flatten())


# %%

def match(case, population, PS, **kwargs):
    """
    match(case, population, PS, **kwargs)

    Pairs individuals of case and populatio.
    
    Parameters
    ----------
    case : Pandas.DataFrame
        Case group.
    population : Pandas.DataFrame
        Population.    
    ovar : string
        Name of the column containing the propensity scores.
    **kwargs: {'pair_name'}
        index     : string (default 'index')
        Name of the index column.
        pair_name : string (default 'pair')
        Name of the pair column.
        caliper   : float (default 0.2 * ((var_s**2 + var_p**2)/2)**0.5)
        Caliper size.
        scale     : float (default 1.0)
        Scale multiplier for the caliper size
           
    Returns
    -------
    Pandas.DataFrame
        Control group.
    """ 
    
    
    pd.options.mode.chained_assignment = None
    
    logit = lambda lp: np.log10(lp/(1-lp))
    
    scale = kwargs['scale'] if 'scale' in kwargs else 1.0
    pair_name = kwargs['pair_name'] if 'pair_name' in kwargs else 'pair'
    _id = kwargs['index'] if 'index' in kwargs else 'index'
        
    sample_size = case.shape[0]

    v_sample = logit(case[PS].values.reshape(-1, 1))
    v_population = logit(population[PS].values.reshape(-1, 1))

    var_s = np.var(v_sample)
    var_p = np.var(v_population)
    
    caliper_size = kwargs['caliper'] if 'caliper' in kwargs else 0.2 * ((var_s**2 + var_p**2)/2)**0.5
        
    v_result = np.empty(case.shape[0])
    v_result[:] = np.nan
    v_result = v_result.reshape(-1, 1)

    d = sklearn.metrics.pairwise.pairwise_distances(v_sample, v_population, metric='manhattan')
    d_idx_sorted = np.argsort(d, axis=1)

    matching_order = list(range(sample_size))
    random.shuffle(matching_order)

    for treated in matching_order:
        actual = 0
        while np.any(v_result == d_idx_sorted[treated][actual]):
            actual = actual + 1
        untreated = d_idx_sorted[treated][actual]
        if abs(d[treated][untreated]) < scale*caliper_size:
            v_result[treated] = untreated
        else:
            v_result[treated] = None
    
    otpt = population.reset_index().rename(columns={population.index.name:_id})
    
    v_result = v_result.flatten()
    case['tmp'] = v_result
    
    v_result = v_result[~np.isnan(v_result)]
    control = otpt.loc[v_result]#.set_index('orig_idx')
    
    for index, row in case.iterrows():
        try:
            case.loc[index, pair_name] = control.loc[case.loc[index, 'tmp'], _id]
        except KeyError:
            case.loc[index, pair_name] = None
    
    del case['tmp']
    
    return control.dropna().set_index(_id), scale * caliper_size

