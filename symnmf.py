import math
import sys
import copy
import numpy as np

np.random.seed(0)

def calculate_average(W, num_of_elements):
    sum = 0.0
    for i in range(num_of_elements):
        for j in range(num_of_elements):
            sum += W[i][j]
    average = sum/(num_of_elements**2)
    return average

def initialize_H(W, num_of_elements, k):
    average = calculate_average(W, num_of_elements)
    H = np.random.uniform(low=0, high=2*math.sqrt(average), size=(num_of_elements, k))
    return H