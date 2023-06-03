import _flute
import numpy as np
import matplotlib.pyplot as plt



a = [2, 0, 3, 1, 4, 6]
b = [0, 1, 2, 3, 5, 3]
_flute.readLUT()
t = _flute.flute(len(a), a, b, 3)

_flute.printtree(t)
arr = _flute.treeToPairArray(t)
print(arr.x)
print(arr.y)
plt.figure(figsize=(7., 7.))

# plotting nodes:
plt.scatter(a, b, s=10, color='r')

for i in range(int(len(arr.x)/2)):
    plt.plot(arr.x[2*i:2*i+2], arr.y[2*i:2*i+2], color='k')    

#plt.plot([2, 2], [0, 1], color='k')
# plotting MST edges:
# plt.plot(arr.x,
#          arr.y,
#          color='k')

plt.xlim(float(min(a)-2), float(max(a)+2))
plt.ylim(float(min(b)-2), float(max(b)+2))
plt.xlabel(r'$X$', size=16)
plt.ylabel(r'$Y$', size=16)
plt.tight_layout()
plt.show()



#0 0, 3 0