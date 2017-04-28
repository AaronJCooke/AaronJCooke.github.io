__author__ = 'Aaron'
import numpy as np
import math
import scipy as sp
def tauchen(s,rho): #approximates an AR process with s state s by s dimensional Markov Chain
    sigma_eps = .37# variance of epsilon in y(t+1)=a + rho*y(t) + epsilon
    m = 3 #allowable standard deviation of the state variable from the mean
    sig_y = sigma_eps * (1-(rho**2))**(-0.5) # Standard deviation of the the state variable
    y_state_0 = -1*sig_y*m # The state with the smallest value
    y_state_n = sig_y*m # The state with the largest value
    y_state = np.linspace(y_state_0, y_state_n, num=s) # Equally spaced state space
    w = y_state[s-1] - y_state[s-2] # The difference between consecutive state values
    P_state = sp.zeros((s,s))#Matrix of transition probability
    phi = lambda x:(1.0 + math.erf(x / math.sqrt(2.0))) / 2.0 #'Cumulative distribution function for the standard normal distribution'
    #print(phi(-3), phi(0), phi(3), 'should be approximately zero, 0.5, and 1')
    print('w',w)
    print('sigma_eps',sigma_eps)
    print('rho',rho)
    print('y_state_0',y_state_0)
    print('y_state_n',y_state_n)
    print('y_state',y_state)
    print('phi',phi(.5))
    for i in range(s):
        for k in range(s):
            if k==0:
                P_state[i,k]=phi((y_state_0 - rho*y_state[i] + 0.5*w)/sigma_eps) # see equation 3b of Tauchen 1986
            elif k == s:
                P_state[i,k]= 1.0 - phi((y_state_n - rho*y_state[i] - 0.5* w)/sigma_eps) # see equation 3b of Tauchen 1986
            else:
                P_state[i,k]=phi((y_state[k] - rho*y_state[i] + 0.5*w)/sigma_eps) - phi((y_state[k] - rho*y_state[i] - 0.5*w)/sigma_eps)
    return P_state #, y_state # " Note: the TM is the [0] and the state is [1] element of the result"
a = tauchen(7, .4)
#with open('C:/Users/Aaron/Desktop/tauchen.txt', 'w') as f:
 #   f.write('\n'.join([str(x) for x in a]))
#f.close()
q = sp.zeros((1, 7))
q = np.sum(a, 1)
#print(q)
#print(a[1,1])
