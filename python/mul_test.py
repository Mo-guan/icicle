modulus = 2**64 - 2**32 + 1


def karatsuba(a, b, m):
    m2 = m//2
    a_lo = a & ((1 << m2) - 1)
    a_hi = a >> m2
    b_lo = b & ((1 << m2) - 1)
    b_hi = b >> m2
    r_lo = a_lo * b_lo
    r_hi = a_hi * b_hi
    diffs_lo = a_hi - a_lo
    diffs_hi = b_lo - b_hi
    carry1 = 0
    carry2 = 0
    if diffs_lo < 0:
        diffs_lo = diffs_lo & ((1 << m2) - 1)
        carry1 = 1
    if diffs_hi < 0:
        diffs_hi = diffs_hi & ((1 << m2) - 1)
        carry2 = 1
    middle_part = diffs_lo * diffs_hi + r_lo + r_hi
    if carry1:
        middle_part -= (diffs_hi << m2)
    if carry2:
        middle_part -= (diffs_lo << m2)
    if carry1 and carry2:
        middle_part = middle_part + 2**m
    middle_part_golden = (a_hi - a_lo) * (b_lo - b_hi) + r_lo + r_hi
    r = r_lo + r_hi * (1 << m) + middle_part * (1 << m2)
    print(f"a = {hex(a)}")
    print(f"b = {hex(b)}")
    print(f"a_lo = {hex(a_lo)}")
    print(f"a_hi = {hex(a_hi)}")
    print(f"b_lo = {hex(b_lo)}")
    print(f"b_hi = {hex(b_hi)}")
    print(f"r_lo = {hex(r_lo)}")
    print(f"r_hi = {hex(r_hi)}")
    print(f"diffs_lo = {hex(diffs_lo)}")
    print(f"diffs_hi = {hex(diffs_hi)}")
    print(f"carry1 = {hex(carry1)}")
    print(f"carry2 = {hex(carry2)}")
    print(f"middle_part = {hex(middle_part)}")
    print(f"r = {hex(r)}\n")

    debug1 = (a_hi - a_lo) * (b_lo - b_hi + 2**32) + a_lo*b_lo + a_hi*b_hi
    debug3 = debug1 - (a_hi-a_lo) * 2**32
    print(
        f"debug1 = {hex((a_hi - a_lo + 2**32) * (b_lo - b_hi + 2**32) + a_lo*b_lo + a_hi*b_hi)}")
    # print(f"debug2 = {hex(debug2)}")
    print(f"debug3 = {hex(debug3)}")
    assert middle_part == middle_part_golden
    assert r == a * b


def reduce(x):
    assert x >= 0 and x < modulus ** 2
    n0 = x & (2**64 - 1)
    n1 = (x >> 64) & (2**32 - 1)
    n2 = x >> 96
    print(f"n0 = {hex(n0)}")
    print(f"n1 = {hex(n1)}")
    print(f"n2 = {hex(n2)}")
    n0minusn2 = n0 - n2
    print(f"n0minusn2 = {hex(n0minusn2)}")
    if n0minusn2 < 0:
        n0minusn2 = n0minusn2 + modulus
    print(f"n0minusn2 = {hex(n0minusn2)}")
    n1new = (2**32-1)*n1
    print(f"(2^32 - 1) n1 = {hex(n1new)}")
    sum = n0minusn2 + n1new
    print(f"sum = {hex(sum)}")
    print(f"carry = {sum < modulus}")
    if sum >= modulus:
        sum = sum - modulus
    print(f"r = {hex(sum)}")
    assert sum == x % modulus


xy_list = [0] * 1024


def check_multiplication(input_file):
    cnt = 0
    with open(input_file, 'r') as file:
        for line in file:
            cnt += 1
            if 'xs: ' in line and 'ys: ' in line and 'xy: ' in line:
                line = line.strip()
                print(f"行{cnt} '{line}'")
                gid = int(line.split('gid: ')[1].split(' ')[0])
                xs = line.split('xs: ')[1].split(' ')[0]
                ys = line.split('ys: ')[1].split(' ')[0]
                xy = line.split('xy: ')[1].split(' ')[0]

                xs_int = int(xs, 16)  # 将xs转换为十进制整数
                ys_int = int(ys, 16)  # 将ys转换为十进制整数
                xy_int = int(xy, 16)  # 将xy转换为十进制整数

                karatsuba(xs_int, ys_int, 64)  # 调用karatsuba函数进行乘法运算

                if xs_int * ys_int == xy_int:
                    # print(f"行 '{line}' 中的 xs * ys = xy 成立")
                    pass
                else:
                    print(f"行{cnt} '{line}' 中的 xs * ys != xy 不成立")
                    print(
                        f"x = {hex(xs_int)}\ny = {hex(ys_int)}\nxy  = {hex(xy_int)}\nx*y = {hex(xs_int * ys_int)}\n")
                xy_list[gid] = xy_int

            if 'reduce_res: ' in line:
                line = line.strip()
                gid = int(line.split('gid: ')[1].split(' ')[0])
                r = line.split('reduce_res: ')[1].split(' ')[0]
                r_int = int(r, 16)
                if r_int == xy_list[gid] % modulus:
                    print(
                        f"pass reduce_res  = {hex(r_int)}\nxy = {hex(xy_list[gid])}\n")
                    pass
                else:
                    print(f"行{cnt} '{line}' 中的 reduce_res == xy mod p 不成立")
                    print(f"xy = {hex(xy_list[gid])}")
                    print(f"reduce_res  = {hex(r_int)}")
                    print(
                        f"reduce_res correct  = {hex(xy_list[gid] % modulus)}")
                reduce(xy_list[gid])


# 用法示例：
# input_file = '/scorpio/home/wangcheng/icicle/icicle/build/Testing/Temporary/LastTest.log'  # 替换为你的输入文件路径
# check_multiplication(input_file)
reduce(0x0285fbbd854d0d89b2c17dad759d91f6)
