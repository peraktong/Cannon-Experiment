import numpy as np


"""
ar = np.array([1,2,3,4,5,6])

mask = np.array([0,1,1,0,0,0],dtype=bool)


mask = 2*b > a+c
mask_2 = 3*b>a+c

mask_3 = mask*mask_2

m = np.array(mask,dtype=int)

print(m, mask)



"""


a = np.array([1,2,3,3,2,2])
b = np.array([2,3,4,2,5,4])
c = np.array([3,4,2,2,1,3])

fusion = np.c_[a,b]
fusion = np.c_[fusion,c]

target = list(set(a))

for i in range(0,len(target)):
    index = np.where(a == target[i])
    index = np.array(index)
    index = index.ravel()
    print(index)



