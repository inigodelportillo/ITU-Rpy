# -*- coding: utf-8 -*-
"""
Turns a rectangle into a minimum number of squares. 

Created on Thu Oct 15 17:24:31 2020

@author: MAW32652
"""

def minSq(m, n, r, sqList = []):
    
    # print("this is m:")
    # print(m)
    # print()
    
    # print("this is n:")
    # print(n)
    # print()
    
    # print("this is r:")
    # print(r)
    # print()
    
    #if the input lists are of the same length
    if len(m) == len(n):
        #append the add the results to the list of results.
        sqList.append(r)
        print("Done!")
        print("Your List of square's is: " + str(sqList))
        #return the lsit of results. 
        return sqList
    
    #case where the m dimension is larger than the n dimension
    elif len(m) > len(n):
        #append the results for m - n and all of n
        sqList.append(r[:len(n), :len(n)])
        
        #do minSq on the rest of the rectangle
        #start at len(n) to the end and all of n. 
        return minSq(m[len(n):], n, r[:, len(n):], sqList)
    
    #case where the n dimension is larger than the m dimension.
    elif len(n) > len(m):
        #append the results for all of m and n - m.
        sqList.append(r[:len(m), :len(m)])
        
        #do minSq on the rest of the rectangle
        #start at all of m and at len(m)
        return minSq(m, n[len(m):], r[len(m):, :], sqList)
    

    