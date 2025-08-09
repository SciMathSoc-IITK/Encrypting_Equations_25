#Importing Libraries
import numpy as np
import math
import random









# Function for converting string to integer
def int_conversion(cipher, n):

    #Initializing the intger
    num = np.zeros(n, dtype=int)

    #Converting char to int
    for i in range(n):
        if cipher[i] == " ":
            num[i] = 26
        else:
            num[i] = ord(cipher[i]) - ord('A')
    return num





# Function for converting integer to string
def text_conversion(num):

    #Intializing the string
    text = ""

    #converting each int to char
    for val in num:
        val = int(val) % 27  
        if val == 26:
            text += " "
        else:
            text += chr(ord('A') + val)
    return text





# Finding inverse of integer a modulo m using itirative method
def modinv(a, m):

    for i in range(1, m):

        if (a * i) % m == 1:
            return i
    




# Modular inverse of  matrix A mod 27
def mod_matrix_inverse(A, mod):

    det = int(round(np.linalg.det(A))) % mod
    det_inv = modinv(det, mod)

   # Finding adjugate/adjoint matrix of A 
    cofactors = np.zeros_like(A) 

    # First caluclate cofactor matrix
    for row in range(3):

        for col in range(3):

            minor = np.delete(np.delete(A, row, axis=0), col, axis=1)
            cofactors[row, col] = ((-1) ** (row + col)) * int(round(np.linalg.det(minor)))

    adj = cofactors.T % mod 
    inverse = (det_inv * adj) % mod # This is inverse matrix of A modulo 27
    return inverse

# We could have used the following to directly calculate adjugate/adjoint:

    # adj = np.round(np.linalg.det(A) * np.linalg.inv(A)).astype(int) % 27
    # Then, inverse= (det_inv*adj)
    # But it is not recommended it might give error with ill-conditoned matrices




# Encryption
def encryption(A, P):

    n = len(P)
    num = int_conversion(P, n)

    # implement padding with space character in multiples of 3
    while len(num) % 3 != 0:
        num = np.append(num, 26)  


    cipher_nums = []

    for i in range(0, len(num), 3):

        block = num[i:i+3]

        # C=A.P
        c = np.dot(A, block) % 27 

        cipher_nums.extend(c)

    encrypted_text = text_conversion(cipher_nums) 
    
    # Obtaining ciphertext
    return encrypted_text









# Decryption
def decryption(A, ciphertext):
    
    n = len(ciphertext)
    num = int_conversion(ciphertext, n) 
    
    # implement padding with space character in multiples of 3
    while len(num) % 3 != 0:
        num = np.append(num, 26) 
        

    # finding inverse matrix of A modulo 27
    A_inv = mod_matrix_inverse(A, 27) 
    

    plain_nums = []

    for i in range(0, len(num), 3):

        # obtaining (3,1) vector as one block
        block = num[i:i+3] 

        # P= A^-1 .C
        p = np.dot(A_inv, block) % 27 

        plain_nums.extend(p)

    decrypted_text = text_conversion(plain_nums) 

    return decrypted_text














# Implementing these functions to encrypt messages
P = str(input("Enter a message:\n")).upper()

A = np.array([
    [14, 0, 11],
    [1, 2, 8],
    [3, 21, 5]
],dtype=int)

Cipher = encryption(A, P)

print("Encrypted Text:", Cipher)

if math.gcd(int(round(np.linalg.det(A))), 27) == 1: 

    # check if gcd(det(A),m)=1
    Decrypted = decryption(A, Cipher)

    print("Decrypted Text:", Decrypted)

else:
    
    print("Decryption  is not possible")

