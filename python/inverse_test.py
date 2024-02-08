modulus = 2**64 - 2**32 + 1


def split_int(line, tag):
    return int(line.split(tag)[1].split(' ')[0], 16)


def inverse(x, trace):
    idx = 0
    u = x
    v = modulus
    b = 1
    c = 0
    while not (u == 1) and not (v == 1):
        while u % 2 == 0:
            u >>= 1
            if b % 2 == 1:
                b += modulus
            b >>= 1
            print("u in while", hex(u))
            print("b in while", hex(b))
            u_trace = trace[idx]
            b_trace = trace[idx + 1]
            idx += 2
            if u != u_trace:
                print("u: ", hex(u), "u_trace: ", hex(u_trace))
                assert False, "u in while not equal"
            if b != b_trace:
                print("b: ", hex(b), "b_trace: ", hex(b_trace))
                assert False, "b in while not equal"

        while v % 2 == 0:
            v >>= 1
            if c % 2 == 1:
                c += modulus
            c >>= 1
            print("v in while", hex(v))
            print("c in while", hex(c))
            v_trace = trace[idx]
            c_trace = trace[idx + 1]
            idx += 2
            if v != v_trace:
                print("v: ", hex(v), "v_trace: ", hex(v_trace))
                assert False, "v in while not equal"
            if c != c_trace:
                print("c: ", hex(c), "c_trace: ", hex(c_trace))
                assert False, "c in while not equal"
        if v < u:
            u = u - v
            b = b - c
            if b < 0:
                b += modulus
        else:
            v = v - u
            c = c - b
            if c < 0:
                c += modulus
        print("u", hex(u))
        print("v", hex(v))
        print("b", hex(b))
        print("c", hex(c))
        u_trace = trace[idx]
        v_trace = trace[idx + 1]
        b_trace = trace[idx + 2]
        c_trace = trace[idx + 3]
        idx += 4
        if u != u_trace:
            print("u: ", hex(u), "u: ", hex(u_trace))
            assert False, "u not equal"
        if v != v_trace:
            print("v: ", hex(v), "v_trace: ", hex(v_trace))
            assert False, "v not equal"
        if b != b_trace:
            print("b: ", hex(b), "b_trace: ", hex(b_trace))
            assert False, "b not equal"
        if c != c_trace:
            print("c: ", hex(c), "c_trace: ", hex(c_trace))
            assert False, "c not equal"
    if u == 1:
        print("xs_inv: ", hex(b))
        return b
    else:
        print("xs_inv: ", hex(c))
        return c


def quick_pow(a, b):
    if b == 0:
        return 1
    ans = 1
    while b != 0:
        if b & 1:
            ans *= a
        b >>= 1
        a *= a

        ans = ans % modulus
        a = a % modulus
    return ans


cnt = 0
xs_list = []
trace_list = []

tag_list = ['u: ', 'v: ', 'b: ', 'c: ',
            'u in while: ', 'b in while: ', 'v in while: ', 'c in while: ']


def tag_in_line(line):
    for tag in tag_list:
        if tag in line:
            return True, tag
    return False, ''


def check(input_file):
    with open(input_file, 'r') as file:
        tl = []
        for line in file:
            if 'xs: ' in line:
                if tl != []:
                    trace_list.append(tl)
                    tl = []
                xs_int = split_int(line, 'xs: ')
                xs_list.append(xs_int)

            is_in, tag = tag_in_line(line)
            if is_in:
                tl.append(split_int(line, tag))

    if tl != []:
        trace_list.append(tl)
        tl = []
    for xs_int, trace in zip(xs_list, trace_list):
        xs_inv_golden = quick_pow(xs_int, modulus - 2)
        xs_inv = inverse(xs_int, trace)
        assert xs_inv == xs_inv_golden

        assert (xs_inv * xs_int) % modulus == 1


input_file = '/scorpio/home/wangcheng/icicle/icicle/build/Testing/Temporary/LastTest.log'  # 替换为你的输入文件路径
check(input_file)

print("0x380e8b36edca88af * 0xd476feb5f5c899ce = ",
      hex(0x380e8b36edca88af * 0xd476feb5f5c899ce % modulus))
print("0xd476feb5f5c899ce * 0x27892a949ba5d73a = ",
      hex(0xd476feb5f5c899ce * 0x27892a949ba5d73a % modulus))
print("0x1056d5ad8b645f77 * 0x27892a949ba5d73a = ",
      hex(0x1056d5ad8b645f77 * 0x27892a949ba5d73a % modulus))
