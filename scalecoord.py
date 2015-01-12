# -*- coding: utf-8 -*-
"""
Created on Sun Jan 11 00:46:13 2015

@author: Administrator
"""

import numpy as np, pandas as pd

def scalecoord(fname):
    df = pd.read_csv(fname, header=None)
    for i in range(3):
        df[i+1] *= 1e-3
    df.to_csv('scaled_'+fname, header=None, index=None)
    return


if __name__ == '__main__':
    scalecoord('temp.csv')