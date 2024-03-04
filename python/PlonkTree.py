from Field import FieldElement as Goldilocks
from Poseidon import poseidon, hash_or_noop, two_to_one
import math

digest_length = 4


def split_at_mut(buf: list, length: int):
    assert length <= len(buf)
    return buf[:length], buf[length:]


def split_last_mut(buf: list):
    assert len(buf) > 0
    return buf[-digest_length:], buf[:-digest_length]


def split_first_mut(buf: list):
    assert len(buf) > 0
    return buf[:digest_length], buf[digest_length:]


def fill_subtree(digests: list[Goldilocks], leaves: list[list[Goldilocks]]):
    assert len(leaves) == len(digests) // digest_length // 2 + 1
    if len(digests) == 0:
        hashout = hash_or_noop(leaves[0])
        return hashout

    left_digests_buf,    right_digests_buf = split_at_mut(digests, len(digests) >> 1)
    left_digest_mem,    left_digests_buf = split_last_mut(left_digests_buf)
    right_digest_mem,    right_digests_buf = split_first_mut(right_digests_buf)

    left_leaves, right_leaves = split_at_mut(leaves, len(leaves) >> 1)

    left_digest = fill_subtree(left_digests_buf, left_leaves)
    right_digest = fill_subtree(right_digests_buf, right_leaves)

    for i in range(digest_length):
        left_digest_mem[i].assign(left_digest[i])
        right_digest_mem[i].assign(right_digest[i])

    return two_to_one(left_digest, right_digest)


def fill_digest_buf(digests: list[Goldilocks], cap: list[Goldilocks], leaves: list[list[Goldilocks]], cap_height: int):
    if len(digests) == 0:
        assert len(cap) // digest_length == len(leaves)
        for i in range(len(cap) // digest_length):
            hashout = hash_or_noop(leaves[i])
            for j in range(digest_length):
                cap[i * digest_length + j].assign(hashout[j])
        return

    subtree_digests_len = len(digests) >> cap_height
    subtree_leaves_len = len(leaves) >> cap_height

    digests_chunks = [digests[i * subtree_digests_len: (i + 1) * subtree_digests_len] for i in range(2 ** cap_height)]
    leaves_chunks = [leaves[i * subtree_leaves_len: (i + 1) * subtree_leaves_len] for i in range(2 ** cap_height)]

    for i in range(2 ** cap_height):
        hash_out = fill_subtree(digests_chunks[i], leaves_chunks[i])
        for j in range(digest_length):
            cap[i * digest_length + j].assign(hash_out[j])


def build_plonk_tree(leaves: list[list[Goldilocks]], cap_height: int):
    log2_leaves_len = int(math.log2(len(leaves)))
    assert cap_height <= log2_leaves_len

    num_digests = 2 * (len(leaves) - (1 << cap_height))
    digests = [Goldilocks(0) for _ in range(num_digests * digest_length)]

    len_cap = 1 << cap_height
    cap = [Goldilocks(0) for _ in range(len_cap * digest_length)]

    fill_digest_buf(digests, cap, leaves, cap_height)

    return digests, cap


if __name__ == "__main__":
    height = 8
    len_leaf = 20
    cap_height = 4

    n = 2 ** (height - 1)
    leaves = [[Goldilocks(i + j * len_leaf) for i in range(len_leaf)] for j in range(n)]

    digests, cap = build_plonk_tree(leaves, cap_height)

    print("Digests: ", digests)
    print("Cap: ", cap)
