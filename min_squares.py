# -*- coding: utf-8 -*-
"""
Turns a rectangle into a minimum number of squares. 

Created on Thu Oct 15 17:24:31 2020

@author: MAW32652
"""

def tiler(m, n, dimList = []):
    """
    Tiler outputs the size of the tiles (squares), that are required to fill 
    a rectangle of arbitrary dimension. 
    
    Note: It does not find the minimum number of squares, but finds a number of squares. 
    
    m = len of horz dimension
    n = len of vert dimension
    """
    
    #if the input lists are of the same length
    if len(m) == len(n):
        #add the length of either dimension to the list. 
        dimList.append(len(m))
        #return the lsit of results. 
        return dimList
    
    #case where the m dimension is larger than the n dimension
    elif len(m) > len(n):
        #append the length of the smaller dimension, n.
        dimList.append(len(n))
        
        #do tiler on the rest of the rectangle
        #start at len(n) to the end and all of n. 
        return tiler(m[len(n):], n, dimList)
    
    #case where the n dimension is larger than the m dimension.
    elif len(n) > len(m):
        #append the length of the smaller dimension, m.
        dimList.append(len(m))
        
        #do tiler on the rest of the rectangle
        #start at all of m and at len(m)
        return tiler(m, n[len(m):], dimList)
    

    