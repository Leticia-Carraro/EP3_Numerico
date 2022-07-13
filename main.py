import numpy as np
from math import e, pi, cos, sin, log, log10, log2, tan, sqrt, tanh, cosh, acos, asin, atan, sinh
import functions as f

print("EP 3 - Numérico\nModelagem de um Sistema de Resfriamento de Chips\n")
print("Indique o método que gostaria de utilizar para resolver seu problema, apenas digite o número entre colchetes:")
chosen_mode = input(" [1] - MÉTODO DE RITZ RALEIGH CÚBICO\n [2] - MÉTODO DAS DIFERENÇAS FINITAS LINEAR\n [3] - MÉTODO DOS ELEMENTOS FINITOS\n")
print("Agora, serão pedidos os dados do problema, responda na ordem solicitada, de acordo com o enunciado do problema:")
n = int(input("n = "))
k = input("k(x) = ")
function = input("q(x) = ")
a = 0
b = int(input("L = "))
u_0 = float(input("u(0) = "))
u_L = float(input("u(L) = "))
if chosen_mode == "1":
    x, y = f.cubic_spline_rayleigh_ritz(n, function)
    print("Os valores de x foram:\n",x)
    print("Os valores de y foram:\n",y)
if chosen_mode == "2":
    x, y = f.linear_finite_diference(a, b, str(2300*750), k, function, a, b, n)
    print("Os valores de x foram:\n",x)
    print("Os valores de y foram:\n",y)
if chosen_mode == "3":
    x, y = f.finite_elements(n, function)
    print("Os valores de x foram:\n",x)
    print("Os valores de y foram:\n",y)