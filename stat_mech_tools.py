import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# define function
def func(w):
	num = np.power(w,2)
	den = np.exp(w) - 1
	return np.divide(num,den)

def monte_carlo_integration(function, lower_x, upper_x, upper_y, num_iter):
	# set grid limits from observation
	# can be modified to take in linspace and calculate max

	# upper y limit
	LY = upper_y

	# lower x limit
	a = lower_x # will give error
	b = upper_x # is ~0 outside

	# area of grid
	grid_area = abs(b - a)*abs(LY - 0)

	# set parameters
	N = num_iter

	# make random array of size N
	x_range = a + abs(b-a)*np.random.random(N)
	y_range = LY*np.random.random(N) # for completeness

	# MC iteration
	counts = 0 # initialize
	for i in range(N):
		x = x_range[i]
		y = y_range[i]
		check = y - function(x)
		if check <= 0:
			# this means that the point is inside
			counts += 1
		elif check > 0:
			# this means that the point is outside
			pass

	# add up all of the counts
	M = counts

	# compute the "area" 
	fraction = M/N
	area = fraction*grid_area

	return area

def average_integrals(areas, n):
	# simply computes an average over n of the integrals
	avg = np.sum(areas)/n
	return avg

def std_integrals(areas, n):
	# calculates sqrt(<(I_1 - <I_1>)^2>)
	mean_area = average_integrals(areas, n)
	diff = areas - mean_area
	diff_sqd = np.power(diff,2)
	mean_diff_sqd = average_integrals(diff_sqd, n)
	final_rms = np.sqrt(mean_diff_sqd)
	return final_rms












