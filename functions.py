import numpy as np
from math import e, pi, cos, sin, log


## GAUSS DE DOIS PONTOS
GAUSS_TABLE_X = [0.9061798459, 0.5384693101, 0.0000000000, -0.5384693101, -0.9061798459]
GAUSS_TABLE_C = [0.2369268850, 0.4786286705, 0.5688888889,  0.4786286705,  0.2369268850]

def table_lookup_gauss(target_n):
  return GAUSS_TABLE_X, GAUSS_TABLE_C

def transformacao(a, b):
  x_aux = [0.5*(b-a), 0.5*(a+b)]
  dx_function = "*("+str(0.5*(b-a))+")"
  x_function = str(x_aux[0])+"*t+"+str(x_aux[1])
  return x_function, dx_function

def gauss(function, target_n, a=-1, b=1):
  result = 0
  x, c = table_lookup_gauss(target_n)
  if a!=-1 or b!=1:
    t, dt = transformacao(a, b)
    function = function.replace("x", "("+t+")")
    f = lambda t: eval(str(function+dt))
  else:
    f = lambda x: eval(function)
  for i in range(len(x)):
    result+=c[i]*f(x[i])
  return result

#------------------------------------------------------------

# MATRIZES TRIDIAGONAIS - RESOLUÇÃO

def solving_system_LU(alpha, beta, b, n):
  a = np.ones(n+1)
  a[1] = alpha [1]
  l=np.ones(n) #iniciando o vetor l com 1s por conveniência
  u=np.zeros(n+1) #iniciando o vetor u com 0s por conveniência
  c=np.zeros(n+1)
  l[1] = beta[1]/alpha[1]
  u[1] = b[1]/alpha[1]
  for i in range(n):
    if i > 1:
      a[i] = alpha[i] - beta[i-1]*l[i-1]
      l[i] = beta[i]/a[i]
      u[i] = b[i]-(beta[i-1]*u[i-1])/a[i]
  a[n] = alpha[n] - beta[n-1]*l[n-1]
  u[n] = b[n]-(beta[n-1]*u[n-1])/a[n]
  c[n] = u[n]

  for i in range(n, 0, -1):
    c[i] = u[i] - (l[i]*c[i+1])
  return l, u, c


def decomposicao_abcd(a, b, q_function, p_function, r_function, bound_a, bound_b, n):
  h = (b-a)/(n+1)
  x = a+h
  a_iterations = np.zeros(n+1)
  q = lambda x: eval(str(q_function))

  b_iterations = np.zeros(n+1)
  p = lambda x: eval(str(p_function))

  c_iterations = np.zeros(n+1)

  d_iterations = np.zeros(n+1)
  r = lambda x: eval(str(r_function))
  
  a_iterations[1] = 2+((h**2)*q(x))
  b_iterations[1] = -1+((h/2)*p(x))

  for i in range(2, n, 1):
    x = a+(i*h)
    a_iterations[i] = 2+((h**2)*q(x))
    b_iterations[i] = -1+((h/2)*p(x))
    c_iterations[i] = -1-((h/2)*p(x))
    d_iterations[i] = -(h**2)*r(x)
  
  x = b - h
  a_iterations[n] = 2+((h**2)*q(x))
  c_iterations[n] = -1-((h/2)*p(x))
  d_iterations[n] = (-(h**2)*r(x))+((1-((h/2)*p(x)))*bound_b)

  return a_iterations, b_iterations, c_iterations, d_iterations

def decomposicao_LU_nova(a, b, c, d, n):
  
  l = np.ones(n+1) #iniciando o vetor l com 1s por conveniência
  u = np.zeros(n) #iniciando o vetor u com 0s por conveniência
  z = np.zeros(n+1) #iniciando o vetor auxiliar z com 0s por conveniência
  l[1]=a[1]
  u[1]=b[1]/a[1]
  z[1]=d[1]/l[1]

  for i in range(2, n, 1):
    l[i] = a[i]-(u[i-1]*c[i])
    u[i] = b[i]/l[i]
    z[i] = (d[i]-(c[i]*z[i-1]))/l[i]

  l[n] = a[n]-(u[n-1]*c[n])
  z[n] = (d[n]-(c[n]*z[n-1]))/l[n]

  return l, u, z

def encontrando_y(n, bound_a, bound_b, u, z):
  '''Nessa parte, estou encontrando os valores de x e y a partir das fórmulas fo
     rnecidas no enunciado. '''
  y = np.zeros(n+2)
  y[0] = bound_a
  y[n+1] = bound_b
  y[n] = z[n]

  for i in range(n-1, 0, -1):
    y[i] = z[i]-(u[i]*y[i+1])
  return y


def matrix_a(n, h):
  matrix_a = np.zeros((n+2, n+2))
  for i in range(n+2):
    for j in range(n+2):
      if i == j:
        matrix_a[i][j] = 2
      elif i == j-1 or i == j+1:
        matrix_a[i][j] = -1
      else:
        pass
  return np.dot((1/h),matrix_a)

def finding_d(f, h, n, x_iteration):
  phi = phi_calc(n, x_iteration, h)
  d = np.zeros(n+2)
  for i in range(n+1):
    d[i] = produto_interno(f, phi[i])
  return d

def decomposicao_abc(n, matriz_a):
  '''Essa função faz a decomposição da matriz dada nos vetores a, b
     e c dos elementos da matriz'''
  #O primeiro passo é iniciar os vetores a, b e c com zeros
  a=np.zeros(n) 
  b=np.zeros(n)
  c=np.zeros(n)
  for i in range(n):
    for j in range(n):
      if i == j+1: #Verifique na matriz tridiagonal que os elementos de a aparecem nessas posições
        a[i] = matriz_a[i][j]
      if i == j-1: #Verifique na matriz tridiagonal que os elementos de c aparecem nessas posições
        c[i]=matriz_a[i][j]

  aux=np.diag_indices(n) #pegando os indices diagonais, pois é onde o b se localiza na matriz A
  b=matriz_a[aux] #inserindo os valores em b
  return a, b, c

def decomposicao_LU(a, b, c, n):
  '''Essa função faz a decomposição da matriz dada em L e U segundo
     a fórmula fornecida. Ela recebe os vetores a, b e c, que devem 
     ser np arrays e a dimensão da matriz. Devolve os vetores l e u
     que são também np arrays'''
  l=np.ones(n) #iniciando o vetor l com 1s por conveniência
  u=np.zeros(n) #iniciando o vetor u com 0s por conveniência
  u[0]=b[0] #Dado
  aux_n = np.arange(1, n)
  for i in aux_n:
    l[i]=a[i]/u[i-1] #Implementacao l_i=a_i/u_(i-1)
    u[i]=b[i]-l[i]*c[i-1] #Implementacao u_i=b_i-l_i*c_(i-1)

  return l, u

def encontrando_xy(n, c, d, l, u):
  '''Nessa parte, estou encontrando os valores de x e y a partir das fórmulas fo
     rnecidas no enunciado. n é int, c, d , l e u são np arrays. O retorno x e y
     também são arrays'''
  y=np.empty(n) #Iniciando com empty para rapidez
  x=np.empty(n)
  y[0]=d[0]
  aux_n=np.arange(1,n)
  for i in aux_n:
    y[i]=d[i]-l[i]*y[i-1]
  x[-1]=y[-1]/u[-1]
  for i in range(n-2, -1, -1):
    x[i]=(y[i]-c[i]*x[i+1])/u[i]
  return x, y
#------------------------------------------------------------

# MÉTODO DE RITZ RALEIGH

def get_Q_values(n_dim, h, x_list, x, q = lambda x: eval("0"), p = lambda x: eval("1"), f = lambda x: eval("12*x*(1-x)")):
  Q = np.zeros(7, n_dim+1)
  for i in range(n_dim):
    if i>0:
      function_1 = str("("+x_list[i+1]+"-x)*(x-"+x_list[i]+")*"+q(x))
      Q[1][i] = ((1/(h[i]))**2)*gauss(function_1, 5, x_list[i], x_list[i+1])
      function_2 = str("((x-"+x_list[i-1]+")**2)*"+q(x))
      Q[2][i] = ((1/(h[i-1]))**2)*gauss(function_2, 5, x_list[i-1], x_list[i])
      function_3 = str("(("+x_list[i+1]+"-x)**2)*"+q(x))
      Q[3][i] = ((1/(h[i]))**2)*gauss(function_3, 5, x_list[i], x_list[i+1])
      function_4 = str(p(x))
      Q[4][i] = ((1/(h[i-1]))**2)*gauss(function_4, 5, x_list[i-1], x_list[i])
      function_5 = str("(x-"+x_list[i-1]+"*"+f(x))
      Q[5][i] = (1/(h[i-1]))*gauss(function_5, 5, x_list[i-1], x_list[i])
      function_6 = str("("+x_list[i+1]+"-x)*"+f(x))
      Q[6][i] = (1/(h[i]))*gauss(function_6, 5, x_list[i], x_list[i+1])
  return Q

def get_function_s(x):
  for i in range(len(x)):
    if x[i]<=-2:
      S = lambda x: eval("0")
    elif x[i] > -2 and x[i] <=-1:
      S = lambda x: eval("(1/4)*((2+x)**3)")
    elif x[i] > -1 and x[i] <= 0:
      S = lambda x: eval("(1/4)*(((2+x)**3)-4*((1+x)**3))")
    elif x[i] > 0 and x[i] <= 1:
      S = lambda x: eval("(1/4)*(((2-x)**3)-4*((1-x)**3))")
    elif x[i] > 1 and x[i] <= 2:
      S = lambda x: eval("(1/4)*((2+x)**3)")
    elif x[i] > 2:
      S = lambda x: eval("0")
  return S

def ritz_raleigh(x, x_target, n_given, q = lambda x: eval("0"), p = lambda x: eval("1"), f = lambda x: eval("12*x*(1-x)")):
  h = np.zeros(n_given+1)
  phi = np.zeros(n_given+1)
  Q = np.zeros(7, n_given+1)
  alpha = np.zeros(n_given+1)
  beta = np.zeros(n_given+1)
  b = np.zeros(n_given)

  for i in range(n_given+1):
    h[i] = x[i+1] - x[i]
    if i > 0:
      for i in range(n_given+1):
        if x_target >= 0 and x_target<=x[i-1]:
          phi[i] = 0
        elif x_target> x[i-1] and x_target <=x[i]:
          phi[i] = (x_target - x[i-1])/h[i-1]
        elif x[i] < x_target and x_target<=x[i+1]:
          phi[i] = (x[i+1]-x_target )/h[i]
        elif x[i+1] < x_target and x_target<=1:
          phi[i] = 0
  
  Q = get_Q_values(n_given, h, x, x_target, q, p, f)    
  for i in range(n_given):
    if i > 0:
      alpha[i] = Q[4][i]+Q[4][i+1]+Q[2][i]+Q[3][i]
      beta[i] = Q[1][i]+Q[4][i+1]
      b[i] = Q[5][i]+Q[6][i]
  alpha[n_given] = Q[4][n_given-1]+Q[4][n_given]+Q[2][n_given]+Q[3][n_given]
  b[n_given] = Q[5][n_given]+Q[6][n_given]

  return alpha, beta, b

#------------------------------------------------------------

# MÉTODO DE RITZ RALEIGH CÚBICO

def get_cubic_spline_basis(S, n, x, x_seq, h):
  phi = np.zeros(n+2)
  for i in range(n):
    if i == 0:
      phi[i] = S((x-x_seq[i])/h) - 4*S((x+h)/h)
    if i == 1:
      phi[i] = S((x-x_seq[i])/h) - S((x+h)/h)
    else:
      phi[i] = S((x-x_seq[i])/h)
  phi[n] = S((x-x_seq[n])/h) - S((x+((n+2)*h))/h)
  phi[n+1] = S((x-x_seq[n+1])/h) - 4*S((x-((n+2)*h))/h)
  return phi

def cubic_spline_rayleigh_ritz(n, f = "12*x*(1-x)"):
  h = 1/(n+1)
  x = np.zeros(n+6)
  matriz_a = np.zeros((n+2, n+2))
  d = np.zeros(n+2)
  for i in range(n+2):
    x[i] = i*h
  x[-2] = x[-1] = 0
  x[n+2] = x[n+3] = 1
  S = get_function_s(x)
  phi = get_cubic_spline_basis(S, n, x[0], x, h)

  for i in range(n+2):
    for j in range(i, min(i+3, n+1), 1):
      L = max(x[j-2], 0)
      U = min(x[i+2], 1)
      matriz_a[i][j] = gauss(str(phi[i]*phi[j]), 5, L, U)
      if i!=j: 
        matriz_a[j][i] = matriz_a[i][j]
      if i>=4:
        for j in range(i-3):
          matriz_a[i][j] = 0
      if i<=n-3:
        for j in range(i+4, n+2, 1):
          matriz_a[i][j] = 0
      L = max(x[i-2], 0)
      U = min(x[i+2], 1)
      d[i] = gauss((str(f)+"*"+str(phi[i])), 5, L, U)
  a, b, c = decomposicao_abc(n, matriz_a)
  l, u = decomposicao_LU(a, b, c, n)
  x, y = encontrando_xy(n, c, d, l, u)

  return x, y

#------------------------------------------------------------

# MÉTODO DAS DIFERENÇAS FINITAS LINEAR

def linear_finite_diference(a, b, p_function, q_function, r_function, bound_a, bound_b, n):
  x = np.zeros(n+2)
  h = (b-a)/(n+1)
  a_iterations, b_iterations, c_iterations, d_iterations = decomposicao_abcd(a, b, q_function, p_function, r_function, bound_a, bound_b, n)
  l, u, z = decomposicao_LU_nova(a_iterations, b_iterations, c_iterations, d_iterations, n)
  y = encontrando_y(n, bound_a, bound_b, u, z)
  for i in range(n+2):
    x[i] = a+(i*h)
  
  return x, y

#------------------------------------------------------------

# MÉTODO DOS ELEMENTOS FINITOS

def phi_calc(n, x_iteration, h):
  phi = np.zeros(n+2, dtype=object)
  aux = np.arange(0.0, x_iteration[-1], h/2)
  for i in range(-1, n+1, 1):
    for j in aux:
      if j > x_iteration[i-1] and j < x_iteration[i]:
        phi[i] = str(1/h)+"*(x-"+str(x_iteration[i-1])+")"
      elif j > x_iteration[i] and j < x_iteration[i+1]:
        phi[i] = str(1/h)+"*(x-"+str(x_iteration[i+1])+")"
  return phi

def mesh(n):
  h = 1/(n+1)
  x = np.zeros(n+2)
  for i in range(n+2):
    x[i] = (i*h)
  return h, x

def matrix_a(n, h):
  matrix_a = np.zeros((n+2, n+2))
  for i in range(n+2):
    for j in range(n+2):
      if i == j:
        matrix_a[i][j] = 2
      elif i == j-1 or i == j+1:
        matrix_a[i][j] = -1
      else:
        pass
  return np.dot((1/h),matrix_a)

def finding_d(f, h, n, x_iteration):
  phi = phi_calc(n, x_iteration, h)
  d = np.zeros(n+2)
  for i in range(n+1):
    d[i] = produto_interno(f, phi[i])
  return d


def produto_interno(function_v, function_w, a = 0, b = 1):
  return gauss(str(function_v)+"*"+str(function_w), 5, a, b)


def finite_elements(n, f):
    h, x_iteration = mesh(n)
    matriz_a = matrix_a(n, h)
    a, b, c = decomposicao_abc(n, matriz_a)
    d = finding_d(f, h, n, x_iteration)
    l, u = decomposicao_LU(a, b, c, n)
    x, y = encontrando_xy(n, c, d, l, u)
    return x, y