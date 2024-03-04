from Field import FieldElement as Goldilocks
from Poseidon import poseidon
from PlonkTree import build_plonk_tree

RATE = 8
T = 12


def check_digests(digests, digest_length, name=''):
    if name:
        print(name + ' = ', end='')
    num_digests = len(digests) // digest_length
    lut = [0 for _ in range(num_digests)]
    for i in range(len(digests)):
        if not digests[i].is_zero():
            lut[i // digest_length] = 1

    print('[', end='')
    for i in range(num_digests):
        if lut[i]:
            print(i, end=', ' if i < num_digests - 1 else ']\n')
        else:
            print('X', end=', ' if i < num_digests - 1 else ']\n')


def build_merkle_subtree(state, subtree_idx, height, cap_height, big_tree_digests, start_segment_size, start_segment_offset, digest_length):
    subtree_height = height - cap_height
    leaves_size = 2 ** (subtree_height - 1)
    number_of_blocks = leaves_size // 2
    segment_size = start_segment_size
    segment_offset = start_segment_offset
    number_of_leaves = 2 ** (height - 1)
    num_digests = 2 * (number_of_leaves - (1 << cap_height))

    # print("Subtree idx: ", subtree_idx)
    # print("Subtree height: ", subtree_height)
    # print("Number of blocks: ", number_of_blocks)
    # print("Segment size: ", segment_size)

    state_ = [Goldilocks(0) for _ in range(T * number_of_blocks)]
    for i in range(number_of_blocks):
        for j in range(2 * digest_length):
            state_[i * T + j] = state[i * 2 * digest_length + j]
    # check_digests(state_, T, 'state_')

    level = 1
    while number_of_blocks > 0:
        for i in range(number_of_blocks):
            out = poseidon(state_[i * T: (i + 1) * T])
            if number_of_blocks > 1:
                i_next_layer = i >> 1
                right_offset = digest_length if i & 1 else 0
                for j in range(digest_length):
                    state_[i_next_layer * T + right_offset + j] = out[j]
                if i & 1:
                    for j in range(2 * digest_length, T):
                        state_[i_next_layer * T + j] = Goldilocks(0)
            if number_of_blocks > 1:
                addr = addr_to_index_fast(level, subtree_idx * number_of_blocks + i, height, cap_height)
            else:
                addr = num_digests + subtree_idx
            for j in range(digest_length):
                big_tree_digests[addr * digest_length + j].assign(out[j])

        segment_offset += segment_size
        segment_size >>= 1
        subtree_height -= 1
        number_of_blocks >>= 1
        level += 1


def split_merkle_tree(leaves, digests, height, cap_height, digest_length):
    number_of_leaves = 2 ** (height - 1)
    number_of_subtrees = 2 ** cap_height
    subtree_height = height - cap_height
    subtree_leaves_size = 2 ** (subtree_height - 1)

    for i in range(number_of_subtrees):
        subtree_leaves = leaves[i * subtree_leaves_size *
                                digest_length: (i + 1) * subtree_leaves_size * digest_length]
        start_segment_size = number_of_leaves // 2
        build_merkle_subtree(
            subtree_leaves,      # state
            i,                   # subtree_idx
            height,              # height
            cap_height,          # cap_height
            digests,             # big_tree_digests
            start_segment_size,  # start_segment_size
            0,                   # start_segment_offset
            digest_length        # digest_length
        )
        # check_digests(digests, digest_length)


def absorb_leaves(leaves, digests, height, cap_height, leaf_length, digest_length):
    number_of_leaves = 2 ** (height - 1)
    state = [[Goldilocks(0) for _ in range(T)] for _ in range(number_of_leaves)]
    for leaf_idx in range(0, leaf_length, RATE):
        for i in range(number_of_leaves):
            for j in range(min(RATE, leaf_length - leaf_idx)):
                state[i][j] = leaves[i * leaf_length + leaf_idx + j]
            # print(f"states {i} : ", [format(x.value, '016x') for x in state[i]])
            state[i] = poseidon(state[i])
            if leaf_idx + RATE >= leaf_length:
                addr = addr_to_index_fast(0, i, height, cap_height)
                for j in range(digest_length):
                    digests[addr * digest_length + j].assign(state[i][j])


def build_merkle_tree(leaves, digests, height, cap_height, leaf_length, digest_length):
    assert cap_height < height and cap_height >= 0 and height > 0
    assert len(leaves) == 2 ** (height - 1) * leaf_length
    absorb_leaves(leaves, digests, height, cap_height, leaf_length, digest_length)
    # check_digests(digests, digest_length, 'after absorb_leaves')
    if cap_height == height - 1:
        return

    new_leaves = [Goldilocks(0) for _ in range(2 ** (height - 1) * digest_length)]
    for i in range(2 ** (height - 1)):
        addr = addr_to_index_fast(0, i, height, cap_height)
        for j in range(digest_length):
            new_leaves[i * digest_length + j] = digests[addr * digest_length + j]

    # check_digests(new_leaves, digest_length, 'leaves')
    # check_digests(new_digests, digest_length, 'digests')
    split_merkle_tree(new_leaves, digests, height, cap_height, digest_length)
    # check_digests(new_leaves, digest_length, 'after split_merkle_tree')
    # check_digests(new_digests, digest_length, 'after split_merkle_tree')
    # check_digests(digests, digest_length, 'after split_merkle_tree')


def addr_to_index(level, idx, height, cap_height):
    subtree_height = height - cap_height
    if subtree_height == 1:
        return idx

    number_of_leaves = 2 ** (height - 1)
    number_of_nodes = 2 ** (height - level - 1)
    number_of_nodes_in_subtree = number_of_nodes >> cap_height
    subtree_idx = idx // number_of_nodes_in_subtree

    num_digests = 2 * (number_of_leaves - (1 << cap_height))
    digest_buf_offset = (num_digests >> cap_height) * subtree_idx

    number_of_left_nodes = 0

    if not (idx & 1):  # add self child
        number_of_left_nodes += 2 ** (level + 1) - 2
    else:  # add sibling node and its children
        number_of_left_nodes += 2 ** (level + 1) - 1
    level += 1
    idx >>= 1
    while level < subtree_height - 1:
        if idx & 1:
            number_of_left_nodes += 2 ** (level + 1)
        level += 1
        idx >>= 1

    return digest_buf_offset + number_of_left_nodes


def addr_to_index_fast(level, idx, height, cap_height):
    subtree_height = height - cap_height
    if subtree_height == 1:
        return idx

    number_of_leaves = 2 ** (height - 1)
    number_of_nodes = 2 ** (height - level - 1)
    number_of_nodes_in_subtree = number_of_nodes >> cap_height
    subtree_idx = idx // number_of_nodes_in_subtree

    num_digests = 2 * (number_of_leaves - (1 << cap_height))
    digest_buf_offset = (num_digests >> cap_height) * subtree_idx

    number_of_left_nodes = 2 ** (level + 1) - 1 if idx & 1 else 2 ** (level + 1) - 2
    mask = 2 ** (subtree_height - 1) - 1
    number_of_left_nodes += (((idx >> 1) << level + 1) & mask) << 1

    return digest_buf_offset + number_of_left_nodes


if __name__ == "__main__":
    height = 4
    cap_height = 1
    leaf_length = 20
    digest_length = 4

    number_of_leaves = 2 ** (height - 1)
    num_digests = 2 * (number_of_leaves - (1 << cap_height))
    len_cap = 2 ** cap_height
    print("Number of leaves: ", number_of_leaves)
    print("Number of digests: ", num_digests)
    leaves = [Goldilocks(i) for i in range(number_of_leaves * leaf_length)]
    digests = [Goldilocks(0) for _ in range((num_digests + len_cap) * digest_length)]
    build_merkle_tree(leaves, digests, height, cap_height, leaf_length, digest_length)
    print(digests)

    plonk_digests, plonk_cap = build_plonk_tree([leaves[i*leaf_length: (i+1) * leaf_length] for i in range(number_of_leaves)], cap_height)
    for x in [_.value for _ in digests]:
        print(f"{{0x{format(x&0xFFFFFFFF, '08x')}, 0x{format((x>>32)&0xFFFFFFFF, '08x')}}}, ")


    cap = digests[-(2**cap_height * digest_length):]
    digests = digests[:-(2**cap_height * digest_length)]
    for i in range(len(cap)):
        assert cap[i] == plonk_cap[i]

    for idx in range(num_digests):
            for i in range(digest_length):
                assert digests[idx * digest_length + i] == plonk_digests[idx * digest_length + i]

    print("Done")
