import numpy as np

def vecNum(vec):
    return int("".join(str(int(n)) for n in vec), base=2)

def printBinary(mat, width):
    for vec in mat:
        print("0b{:0{}b},".format(vecNum(vec), width))

def genToParityCheck(genParity):
    xpose = genParity.transpose()
    return np.hstack((xpose, np.eye(xpose.shape[0])))

gen = np.array([
    [1, 0, 0, 0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1, 1, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 1, 1, 1, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 1, 1, 1],
    [0, 0, 0, 1, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1, 1, 1],
    [0, 0, 0, 0, 1, 0, 0, 0, 0,   1, 1, 0, 1, 1, 0, 1, 1],
    [0, 0, 0, 0, 0, 1, 0, 0, 0,   1, 1, 1, 1, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 1, 0, 0,   1, 1, 1, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0,   0, 1, 1, 1, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1,   0, 0, 1, 1, 1, 0, 0, 1],
])

parity = gen[:, -8:]

print("generator:")
printBinary(parity.transpose(), 9)

print("parity check:")
parityCheck = genToParityCheck(parity)
printBinary(parityCheck, 17)

error = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
syndromes = {}

for r in range(17):
    w = np.roll(error, r)
    w[-1] = 1

    s = (parityCheck @ w) % 2
    syndromes[vecNum(s)] = vecNum(w)

for x in syndromes.keys():
    assert x < 256

print("syndrome/pattern mappings:")

for syn in range(256):
    try:
        print("\n0b{:017b},".format(syndromes[syn]))
    except KeyError:
        print("0, ", end="")

print()
