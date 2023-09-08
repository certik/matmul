from timeit import default_timer as clock
import numpy as np

n = 1024*4

A = np.empty((n,n), dtype=np.float32)
B = np.empty((n,n), dtype=np.float32)
t1 = clock()
C = np.dot(A, B)
t2 = clock()

GHz = 1e9
t = t2-t1
fma_clock = 0.0625
freq = 3.2*GHz
measured = t * freq / n**3
percent_peak = fma_clock / measured * 100
print(f"Time: {t}s")
print("Clock cycles per element:")
print(f"Theoretical performance peak: {fma_clock} cycles")
print(f"Measured:                     {measured} cycles")
print(f"Percent peak:                 {percent_peak}%")
