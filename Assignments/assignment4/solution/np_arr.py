import numpy as np


A = np.array([[-0.250000000000, 0.250000000000 ,0.250000000000 ,-0.250000000000],
[0.235576968368 ,-0.232604849813 ,0.267395150187 ,-0.200786667995 ,],
[0.179487213995 ,0.211305395813 ,-0.288694604187 ,-0.256876422369 ,],
[-0.269230820992 ,-0.105594457356 ,-0.105594457356 ,0.058041906280]])

b = np.array([-1.272693304105, -0.272693304105, -0.272693304105, 0.727306695895 ]).reshape((4,1))

print(A@b)