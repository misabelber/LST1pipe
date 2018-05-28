import numpy as np
import Disp

a = np.array([1,2,3,4])
b = np.array([0.5,0.3,1,2])
c = np.array([0.1,0.2,0.3,0.4])
d = np.array([0.2,0.3,0.4,0.5])

result = Disp.calc_CamSourcePos(a,b,c,d,30)
