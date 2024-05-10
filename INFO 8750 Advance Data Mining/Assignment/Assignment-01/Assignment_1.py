#factorial function
def factorial(n):
    result = 1
    for i in range(1, n + 1):
        result *= i
    return result
 
#Closing condition
x = 1
while factorial(x) < 1000000:
    x += 1
 
print("The least integer x which satisfies x! >= 1,000,000 is:", x)
print("x! =", factorial(x))