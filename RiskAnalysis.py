import math
from _ast import alias

from scipy.stats import norm
from scipy.stats import binom
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

############################################### Init. Parameters:


########## Paremeters Description #############
# N          The number of the individuals in the beacon service
# a_prime    Database property
# b_prime    Database property
# alpha      False positive probability
# sigma      Sequencing Error
###############################################

alpha = 0.05
sigma = 1e-06

# 1000 Genomes Phase 1 Dataset
N, a_prime, b_prime, = 1092, 0.0735, 1.0096

# SSMP Dataset
# N, a_prime, b_prime = 100 , 0.1848, 0.85

# GoNL Dataset
#N, a_prime, b_prime= 498 , 0.1131, 0.8574

# 1000 Genomes Phase 1 Affymetrix array Dataset
# N, a_prime, b_prime = 1074 , 0.6483, 1.2876


###############################################

a = a_prime + 1
b = b_prime + 1
n = math.ceil(N ** a)


def D(N, a, b):
    return (math.gamma(a + b)) / (math.gamma(b) * (2 * N + a + b) ** a)


teta_0 = n * (1 - D(N, a, b))
variance_0 = n * D(N, a, b) * (1 - D(N, a, b))
teta_1 = n * (1 - sigma * D(N - 1, a, b))
variance_1 = n * sigma * D(N, a, b) * (1 - sigma * D(N, a, b))
# t_alpha=norm.ppf((1-alpha),teta_0,math.sqrt(variance_0))
t_alpha = binom.ppf((1 - alpha), n, (1 - D(N, a, b)))

############### Calculating the number of yes

print("The number of needed queries: ", n)
message = "The number of needed yes to say that the query individual is in the database: "
print(message, t_alpha)


############### Emprical Power computing form the responds of the beacon servic

def emperical_power(numberOfyes):
    return 1 - binom.cdf(numberOfyes, n, 1 - sigma * D(N - 1, a, b))


############### Ploting the two different distibution

t = [i for i in range(int(teta_0) - 5, int(teta_0) + 5)]
Dist_H0 = binom(n, 1 - D(N, a, b))
Dist_H1 = binom(n, 1 - sigma * D(N - 1, a, b))
Y_H0 = [Dist_H0.pmf(i) for i in t]
Y_H1 = [Dist_H1.pmf(i) for i in t]
#
plt.plot(t, Y_H0, t, Y_H1,np.full((10), t_alpha),np.arange(0,1.0,0.1), 'r--', label='t_alpha')
plt.xticks(t, t)
# plt.text(t_alpha-1,.1)
plt.ylabel('Probability')
plt.xlabel('# of Yes in respond to queries')
# plt.annotate('t_alpha', xy=(t_alpha,0), xytext=(teta_0-2, 0.5),arrowprops=dict(arrowstyle="->"))
plt.grid(True)
# tmp = np.arange(0,1,0.1)
print(t)
plt.show()
