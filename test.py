def salam(a, b):
    a = b
    a.remove(2)
    return b

a = [1, 2]
b = [2, 3]
a = b
a.remove(2)
print(b)