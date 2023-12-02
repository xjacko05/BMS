import numpy as np
from pyldpc import make_ldpc, encode, decode, get_message, coding_matrix


def text_to_binary(text):
    binary_array = []

    for char in text:
        # Convert each character to its binary representation
        binary_representation = bin(ord(char))[2:]
        
        # Pad the binary representation to 8 bits (if needed)
        padded_binary = binary_representation.zfill(8)
        
        # Extend the binary array with the bits of the character
        binary_array.extend(map(int, padded_binary))

    return binary_array


length = 6
n = 2*length
d_v = length -1
d_c = length
snr = 100
H, G = make_ldpc(n, d_v, d_c, systematic=False, sparse=True)
#print(H)
#print(G)
k = G.shape[1]
#v = np.random.randint(2, size=k)
#v = [1,0,1]
#G = np.array([[1,0,0],[0,1,0],[0,0,1],[1,0,1],[0,0,1],[0,1,1]])
v = [0,1,1,1,1,0,0,1]
#G = np.array([[0, 0, 1, 1, 0, 1, 0, 0],[1, 0, 0, 0, 0, 0, 0, 0],[0, 1, 1, 0, 0, 0, 0, 1],[0, 0, 1, 1, 1, 1, 0, 1],[1, 0, 0, 1, 0, 0, 1, 0],[1, 1, 1, 1, 0, 0, 0, 1],[0, 0, 0, 0, 0, 1, 0, 0],[1, 0, 0, 0, 1, 1, 1, 1],[0, 1, 1, 1, 1, 1, 1, 1],[0, 1, 0, 0, 0, 0, 0, 0],[0, 0, 1, 0, 0, 0, 0, 0],[0, 0, 0, 1, 0, 0, 0, 0],[0, 0, 0, 0, 1, 0, 0, 0],[0, 0, 0, 0, 0, 1, 0, 0],[0, 0, 0, 0, 0, 0, 1, 0],[0, 0, 0, 0, 0, 0, 0, 1]])
v = np.array(text_to_binary('hello'))
H = np.loadtxt('matica.csv', delimiter=',', dtype=int)
G = coding_matrix(H)
res = encode(G, v, snr)
enc = np.array([0 if elem == 1 else 1 for elem in res])
#print(*[0 if elem == 1 else 1 for elem in res], sep='')
#d = decode(H, res, snr)
x = get_message(G, enc)
assert abs(x - v).sum() == 0