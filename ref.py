from timeit import default_timer as clock
import numpy as np

n = 1024*4*2

A = np.empty((n,n), dtype=np.float32)
B = np.empty((n,n), dtype=np.float32)
t1 = clock()
C = np.dot(A, B)
t2 = clock()
print(f"Time: {t2-t1}s")
