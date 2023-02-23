import rainbow as rb
import numpy as np
import numpy.random as rand
import time

q = 2**8
v1 = 6
n = 33

#Y = rand.randint(0, q, (n - v1,))
#s_k = rb.keygen(q, v1, 12, 17, 22, n)

#print(rb.sign(q, Y, s_k))

q = 127

# (1*6796 + 7*9061 - 1) / 12837

#for i in range(0):
#    print("size:", i)
#    t = time.time()
M = np.random.randint(0, q, (26, 26))
#    print(M)
#    Inv = rb.inv_mod(M, q)
#    print(rb.inv_mod(Inv, q))
#    print(time.time() - t)

print("{", end="")
for x in range(len(M)):
    print("{", end="")
    print(",".join([str(M[x][y]) for y in range(len(M))]), end="")
    print("}" + ("," if x < len(M) - 1 else ""), end="")
print("}", end="")
print()



q = 17
v1 = 3
n = 10

Y = rand.randint(1, q, (n - v1,))
s_k = rb.keygen(q, v1, 8, 9, n)
X = rb.sign(q, Y, s_k)
print(rb.verify(q, s_k[2],X,Y))
#print(rb.sign(q, Y, s_k))

#Y = rand.randint(0, q, (n - v1,))
#s_k = rb.keygen(q, v1, 12, 17, 22, n)

#print(rb.sign(q, Y, s_k))