# a python version of poseidon
from Constants import *
from Field import FieldElement as Goldilocks


class Round_ctr:
    def __init__(self, round_ctr):
        self.round_ctr = round_ctr

    def inc(self):
        self.round_ctr += 1

    def get(self):
        return self.round_ctr


def constant_layer(state, round_ctr):
    for i in range(12):
        if i < WIDTH:
            round_constant = Goldilocks(
                ALL_ROUND_CONSTANTS[i+WIDTH*round_ctr.get()])
            state[i] += round_constant


def sbox_monomial(x):
    x2 = x * x
    x4 = x2 * x2
    x3 = x2 * x
    return x4 * x3


def sbox_layer(state):
    for i in range(12):
        if i < WIDTH:
            state[i] = sbox_monomial(state[i])


def mds_row_shf(r, v):
    res = Goldilocks(0)
    for i in range(12):
        if i < WIDTH:
            res += v[(i+r) % WIDTH] * Goldilocks(MDS_MATRIX_CIRC[i])
    res += v[r] * Goldilocks(MDS_MATRIX_DIAG[r])
    return res


def mds_layer(state):
    result = [Goldilocks(0) for _ in range(WIDTH)]
    for r in range(12):
        if r < WIDTH:
            sum = mds_row_shf(r, state)
            result[r] = sum
    return result


def partial_first_constant_layer(state):
    for i in range(12):
        if i < WIDTH:
            state[i] += Goldilocks(FAST_PARTIAL_FIRST_ROUND_CONSTANT[i])


def mds_partial_layer_init(state):
    result = [Goldilocks(0) for _ in range(WIDTH)]
    result[0] = state[0]
    for r in range(1, 12):
        if r < WIDTH:
            for c in range(1, 12):
                if c < WIDTH:
                    t = Goldilocks(
                        FAST_PARTIAL_ROUND_INITIAL_MATRIX[r - 1][c - 1])
                    result[c] += state[r] * t
    return result


def full_rounds(state, round_ctr):
    for _ in range(HALF_N_FULL_ROUNDS):
        constant_layer(state, round_ctr)
        sbox_layer(state)
        state_ = mds_layer(state)
        for i in range(WIDTH):
            state[i] = state_[i]
        round_ctr.inc()


def mds_partial_layer_fast(state, r):
    d_sum = Goldilocks(0)
    for i in range(1, 12):
        if i < WIDTH:
            t = Goldilocks(FAST_PARTIAL_ROUND_W_HATS[r][i - 1])
            si = state[i]
            d_sum += si * t
    s0 = state[0]
    mds0to0 = MDS_MATRIX_CIRC[0] + MDS_MATRIX_DIAG[0]
    d_sum += s0*Goldilocks(mds0to0)
    result = [Goldilocks(0) for _ in range(WIDTH)]
    result[0] = d_sum
    for i in range(1, 12):
        if i < WIDTH:
            t = Goldilocks(FAST_PARTIAL_ROUND_VS[r][i - 1])
            result[i] = state[i] + state[0] * t
    return result


def partial_rounds(state, round_ctr):
    partial_first_constant_layer(state)
    state_ = mds_partial_layer_init(state)
    for i in range(WIDTH):
        state[i] = state_[i]

    for i in range(N_PARTIAL_ROUNDS):
        state[0] = sbox_monomial(state[0])
        state[0] = state[0] + Goldilocks(FAST_PARTIAL_ROUND_CONSTANTS[i])
        state_ = mds_partial_layer_fast(state, i)
        for ii in range(WIDTH):
            state[ii] = state_[ii]
        round_ctr.inc()


def poseidon(input):
    state = input[:]
    round_ctr = Round_ctr(0)
    full_rounds(state, round_ctr)
    partial_rounds(state, round_ctr)
    full_rounds(state, round_ctr)
    return state


def hash_n_to_m_no_pad(input, num_outputs):
    input_chunks = [input[i:min(i+8, len(input))] for i in range(0, len(input), 8)]
    state = [Goldilocks(0) for _ in range(12)]
    for input_chunk in input_chunks:
        for i in range(len(input_chunk)):
            state[i] .assign(input_chunk[i])
        state = poseidon(state)

    output = []
    while True:
        for s in state[:8]:
            output.append(s)
            if len(output) == num_outputs:
                return output


def hash_or_noop(input):
    if len(input) < 4:
        output = [Goldilocks(0) for _ in range(4)]
        for i in range(len(input)):
            output[i] = input[i]
        return output
    else:
        return hash_n_to_m_no_pad(input, 4)

def two_to_one(left, right):
    state = [Goldilocks(0) for _ in range(12)]
    for i in range(4):
        state[i].assign(left[i])
        state[i+4].assign(right[i])
    return poseidon(state)[:4]

if __name__ == "__main__":
    # state = [Goldilocks(i) for i in range(12)]
    state = [6099940033764713049, 11297638662466000013, 9055972977164648173, 3296733479173602111, 6154433692288734276, 7577668585637567982,
             9651928966405010883, 3856873209381218373, 12956910039412814074, 8941104527614940572, 13836442418435444548, 7731245241932084924]
    state = [Goldilocks(i) for i in state]
    state = poseidon(state)
    print([str(i) for i in state])
