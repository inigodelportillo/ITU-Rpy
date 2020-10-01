# -*- coding: utf-8 -*-
"""
Creates 2D lists of the L_red database files for P.840 revisions & and 8.
 
This will allow the user to compare the differences between revisions using
equality commands in the console after running through the script. 

Variables Names:
    Revision 7:
        rev8_01List = L_red database in revision 7 with probability 0.1
        rev8_1List = L_red database in revision 7 with probability 1
        
    Revision 8:
        rev8_01List = L_red database in revision 8 with probability 0.1
        rev8_1List = L_red database in revision 8 with probability 1

Created on Thu Jun 25 15:00:36 2020

@author: MAW32652
"""
import csv

### Creating the revision 7 variables

with open('C:/Users/MAW32652/Desktop/Revision 7/01test7.csv', newline='') as rev7_01:
    reader = csv.reader(rev7_01)
    rev7_01List = list(reader)

with open('C:/Users/MAW32652/Desktop/Revision 7/1test7.csv', newline='') as rev7_1:
    reader = csv.reader(rev7_1)
    rev7_1List = list(reader)

with open('C:/Users/MAW32652/Desktop/Revision 7/02test7.csv', newline='') as rev7_02:
    reader = csv.reader(rev7_02)
    rev7_02List = list(reader)

with open('C:/Users/MAW32652/Desktop/Revision 7/2test7.csv', newline='') as rev7_2:
    reader = csv.reader(rev7_2)
    rev7_2List = list(reader)

with open('C:/Users/MAW32652/Desktop/Revision 7/03test7.csv', newline='') as rev7_03:
    reader = csv.reader(rev7_03)
    rev7_03List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 7/3test7.csv', newline='') as rev7_3:
    reader = csv.reader(rev7_3)
    rev7_3List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 7/05test7.csv', newline='') as rev7_05:
    reader = csv.reader(rev7_05)
    rev7_05List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 7/5test7.csv', newline='') as rev7_5:
    reader = csv.reader(rev7_5)
    rev7_5List = list(reader)

with open('C:/Users/MAW32652/Desktop/Revision 7/10test7.csv', newline='') as rev7_10:
    reader = csv.reader(rev7_10)
    rev7_10List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 7/20test7.csv', newline='') as rev7_20:
    reader = csv.reader(rev7_20)
    rev7_20List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 7/30test7.csv', newline='') as rev7_30:
    reader = csv.reader(rev7_30)
    rev7_30List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 7/50test7.csv', newline='') as rev7_50:
    reader = csv.reader(rev7_50)
    rev7_50List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 7/60test7.csv', newline='') as rev7_60:
    reader = csv.reader(rev7_60)
    rev7_60List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 7/70test7.csv', newline='') as rev7_70:
    reader = csv.reader(rev7_70)
    rev7_70List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 7/80test7.csv', newline='') as rev7_80:
    reader = csv.reader(rev7_80)
    rev7_80List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 7/90test7.csv', newline='') as rev7_90:
    reader = csv.reader(rev7_90)
    rev7_90List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 7/95test7.csv', newline='') as rev7_95:
    reader = csv.reader(rev7_95)
    rev7_95List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 7/99test7.csv', newline='') as rev7_99:
    reader = csv.reader(rev7_99)
    rev7_99List = list(reader)

### Creating the revision 8 variables

with open('C:/Users/MAW32652/Desktop/Revision 8/01test8.csv', newline='') as rev8_01:
    reader = csv.reader(rev8_01)
    rev8_01List = list(reader)

with open('C:/Users/MAW32652/Desktop/Revision 8/1test8.csv', newline='') as rev8_1:
    reader = csv.reader(rev8_1)
    rev8_1List = list(reader)

with open('C:/Users/MAW32652/Desktop/Revision 8/02test8.csv', newline='') as rev8_02:
    reader = csv.reader(rev8_02)
    rev8_02List = list(reader)

with open('C:/Users/MAW32652/Desktop/Revision 8/2test8.csv', newline='') as rev8_2:
    reader = csv.reader(rev8_2)
    rev8_2List = list(reader)

with open('C:/Users/MAW32652/Desktop/Revision 8/03test8.csv', newline='') as rev8_03:
    reader = csv.reader(rev8_03)
    rev8_03List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 8/3test8.csv', newline='') as rev8_3:
    reader = csv.reader(rev8_3)
    rev8_3List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 8/05test8.csv', newline='') as rev8_05:
    reader = csv.reader(rev8_05)
    rev8_05List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 8/5test8.csv', newline='') as rev8_5:
    reader = csv.reader(rev8_5)
    rev8_5List = list(reader)

with open('C:/Users/MAW32652/Desktop/Revision 8/10test8.csv', newline='') as rev8_10:
    reader = csv.reader(rev8_10)
    rev8_10List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 8/20test8.csv', newline='') as rev8_20:
    reader = csv.reader(rev8_20)
    rev8_20List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 8/30test8.csv', newline='') as rev8_30:
    reader = csv.reader(rev8_30)
    rev8_30List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 8/50test8.csv', newline='') as rev8_50:
    reader = csv.reader(rev8_50)
    rev8_50List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 8/60test8.csv', newline='') as rev8_60:
    reader = csv.reader(rev8_60)
    rev8_60List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 8/70test8.csv', newline='') as rev8_70:
    reader = csv.reader(rev8_70)
    rev8_70List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 8/80test8.csv', newline='') as rev8_80:
    reader = csv.reader(rev8_80)
    rev8_80List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 8/90test8.csv', newline='') as rev8_90:
    reader = csv.reader(rev8_90)
    rev8_90List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 8/95test8.csv', newline='') as rev8_95:
    reader = csv.reader(rev8_95)
    rev8_95List = list(reader)
    
with open('C:/Users/MAW32652/Desktop/Revision 8/99test8.csv', newline='') as rev8_99:
    reader = csv.reader(rev8_99)
    rev8_99List = list(reader)


# =============================================================================
# 
# iturFixed = []
# for sublist in itur7List:
#     for item in sublist:
#         iturFixed.append(item)
#         
# rev7Fixed = []
# rev8Fixed = []
# for sublist in rev7List:
#     for item in sublist:
#         rev7Fixed.append(item)
#         
# for sublist in rev8List:
#     for item in sublist:
#         rev8Fixed.append(item)
# =============================================================================
        

# =============================================================================
# for i in range(len(rev7List)):
#     print(rev7List[i] == itur7List[i])
# 
# print()
# print()
# print(rev7List == itur7List)
# =============================================================================
