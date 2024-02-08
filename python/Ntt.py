import numpy as np
import copy
from Field import FieldElement, Field

g = Field().g()
mod = Field().p


def ntt(a, inv):
    n = len(a)
    bit = 0
    while (1 << bit) < n:
        bit += 1
    rev = [0] * n
    for i in range(0, n):
        rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (bit-1))
        if (i < rev[i]):
            tmp = a[rev[i]]
            a[rev[i]] = a[i]
            a[i] = tmp
    mid = 1
    while mid < n:
        tmp = g ^ ((mod-1)//(mid*2))
        if (inv == 1):
            tmp = tmp ^ (mod-2)

        for i in range(0, n, mid*2):
            omega = FieldElement(1)
            for j in range(mid):
                x = a[i+j]
                y = omega*a[i+j+mid]
                a[i+j] = (x+y)
                a[i+j+mid] = (x-y)
                omega = omega*tmp

        mid = mid * 2
    if (inv == 1):
        tmp = FieldElement(n, Field()).inverse()
        for i in range(0, n):
            a[i] = a[i] * tmp


n = 1 << 3
a0 = [FieldElement(1) for _ in range(n)]
a1 = [FieldElement(1) if i % 2 == 0 else -FieldElement(1) for i in range(n)]
ntt(a0, 0)
ntt(a1, 0)
print([hex(x.value) for x in a0])
print([hex(x.value) for x in a1])
with open('/scorpio/home/wangcheng/icicle/examples/c++/ntt/ntt_output.txt', 'r') as f:
    for i in range(n):
        out = FieldElement(int(f.readline().strip(), 16))
        print(f"line {i} : {a0[i]} {out}")
        assert a0[i] == out
    for i in range(n):
        out = FieldElement(int(f.readline().strip(), 16))
        print(f"line {i+n} : {a1[i]} {out}")
        assert a1[i] == out
