# Hill-Cipher algorithm
import numpy as np
import math

# Convert text to number
def int_conversion(cipher, n):
    num = np.zeros(n, dtype=int)
    for i in range(n):
        if cipher[i] == " ":
            num[i] = 26
        else:
            num[i] = ord(cipher[i]) - ord('A')
    return num

# Convert number to text
def text_conversion(num):
    text = ""
    for val in num:
        val = int(val) % 27  
        if val == 26:
            text += " "
        else:
            text += chr(ord('A') + val)
    return text

# Finding inverse of integer a modulo m
def modinv(a, m):
    for i in range(1, m):
        if (a * i) % m == 1:
            return i
    
# Modular inverse of  matrix A mod 27
def mod_matrix_inverse(A, mod):
    det = int(round(np.linalg.det(A))) % mod
    det_inv = modinv(det, mod)

   # Finding adjugate/adjoint matrix of A 
    cofactors = np.zeros_like(A) # First caluclate cofactor matrix
    for row in range(3):
        for col in range(3):
            minor = np.delete(np.delete(A, row, axis=0), col, axis=1)
            cofactors[row, col] = ((-1) ** (row + col)) * int(round(np.linalg.det(minor)))

    adj = cofactors.T % mod 
    inverse = (det_inv * adj) % mod # This is inverse matrix of A modulo 27
    return inverse

#We could have used the following to directly calculate adjugate/adjoint:
# adj = np.round(np.linalg.det(A) * np.linalg.inv(A)).astype(int) % 27
# Then, inverse= (det_inv*adj)
# But it is not recommended it might give error with ill-conditoned matrices


# Encryption
def encryption(A, P):
    n = len(P)
    num = int_conversion(P, n)

    
    while len(num) % 3 != 0:
        num = np.append(num, 26)  # implement padding with space character in multiples of 3

    cipher_nums = []
    for i in range(0, len(num), 3):
        block = num[i:i+3]
        c = np.dot(A, block) % 27 # C=A.P
        cipher_nums.extend(c)

    encrypted_text = text_conversion(cipher_nums) # Obtaining ciphertext
    return encrypted_text

# Decryption
def decryption(A, ciphertext):
    n = len(ciphertext)
    num = int_conversion(ciphertext, n) # converting ciphertext into blocks of integers

    while len(num) % 3 != 0:
        num = np.append(num, 26) # implement padding with space character in multiples of 3

    A_inv = mod_matrix_inverse(A, 27) # finding inverse matrix of A modulo 27

    plain_nums = []
    for i in range(0, len(num), 3):
        block = num[i:i+3] # obtaining (3,1) vector as one block
        p = np.dot(A_inv, block) % 27 # P= A^-1 .C
        plain_nums.extend(p)

    decrypted_text = text_conversion(plain_nums) # retrieving our original message
    return decrypted_text

# Implementing these functions to encrypt messages
P = str(input("Enter a message:\n")).upper()
A = np.array([
    [1, 2, 5],
    [6, 16, 15],
    [17, 19, 3]
],dtype=int)

Cipher = encryption(A, P)
print("Encrypted Text:", Cipher)

if(math.gcd((int(np.linalg.det(A))),27)==1): # check if gcd(det(A),m)=1
    Decrypted = decryption(A, Cipher)
    print("Decrypted Text:", Decrypted)
else:
    print("Decryption  is not possible")