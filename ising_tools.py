from numpy import *
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("white")
from random import *
from scipy.optimize import curve_fit as opt


def Energy_all_neighbours(spins, i, j, L):
    S = spins[i, j]
    neighbours = spins[(i + 1) % L, j] + spins[(i - 1) %
                                               L, j] + spins[i, (j + 1) % L] + spins[i, (j - 1) % L]
    nn_neighbours = (spins[(i + 1) % L, (j + 1) % L] + spins[(i + 1) % L, (j - 1) % L] + spins[(i - 1) % L, (j + 1) % L]
                     + spins[(i - 1) % L, (j - 1) % L])
    e_ij = - S * (neighbours + nn_neighbours)
    return e_ij

# define energy function


def Energy(spins, L):
    E = 0.0
    for i in range(len(spins)):
        for j in range(len(spins)):
            e_ij = Energy_all_neighbours(spins, i, j, L)
            E += e_ij
    return E / 2.  # to avoid double counting


def Initialize(L, state = 0):
    """
    If spin == 2, fill with random array
    """
    # args is either 1 or -1
    spins = zeros((L, L), dtype=int)

    if state == 0:
        for i in range(L):
            for j in range(L):
                spins[i, j] = choice([-1, 1])

    elif (state == -1 or state == 1):
        for i in range(L):
            for j in range(L):
                spins[i, j] = state

    else:
        print("State (spin value) must be \pm 1 or 0 (random)")
        from sys import exit
        exit()

    return spins


def deltaU(spins, i, j, L):
    e_ij = Energy_all_neighbours(spins, i, j, L)
    Ediff = -2 * e_ij
    return Ediff


def McMove(spins, L, T, cases, retene=False):
    E = Energy(spins, L)
    N = L * L
    for k in range(N):
        itest = int(random() * L)
        jtest = int(random() * L)
        (i, j) = (itest, jtest)
        S = spins[i, j]
        Ediff = deltaU(spins, i, j, L)

        # if del_E < 0, we flip. If not, pick a 0 < random() < 1 and it flips if prob > random()

        # prob > 1 so always flip
        if Ediff <= 0:
            spins[i, j] = -S
            # adjust energy to reflect flip
            E += Ediff

        # otherwise get the prob from cases
        else:
            idx = int((16 - Ediff) / 4)
            prob = cases[idx]
            #prob = exp(-Ediff/T)
            if prob > random():
                spins[i, j] = -S
                # adjust energy to reflect flip
                E += Ediff

    if retene == True:
        return spins, E
    else:
        return spins


def ising_main(MCs, L, T, state = 0, retene=False):
    """
    takes in steps, gridsize, temperature and initial state
    and advances the system by {steps}

    returns: all_spins: (L,L,MCs) array containing the progression of system
             ene: (MCs) array containing progression of energy
    """
    # set up spins
    spins = Initialize(L, state)

    # set up cases
    possible_ediffs = asarray([16, 12, 8, 4])
    cases = exp(-(possible_ediffs) / T)

    # run it
    if retene == True:
        for k in range(MCs):
                # record last state
                # move it forward
            spins, ene = McMove(spins, L, T, cases, retene=True)
        return spins, ene

    elif retene == False:
        for k in range(MCs):
            # record last state
            # move it forward
            spins = McMove(spins, L, T, cases)

        return spins


def ising_main_flipbook(MCs, L, T, state = 0, retene=False):
    """
    takes in steps, gridsize, temperature and initial state
    and advances the system by {steps}

    returns: all_spins: (L,L,MCs) array containing the progression of system
             ene: (MCs) array containing progression of energy
    """
    # set up spins
    spins = Initialize(L, state)

    # set up cases
    possible_ediffs = asarray([16, 12, 8, 4])
    cases = exp(-(possible_ediffs) / T)

    # run it
    all_spins = zeros((L, L, MCs), dtype=int)

    if retene == True:

        ene = zeros((MCs), dtype=float)

        for k in range(MCs):
            # record last state
            all_spins[:, :, k] = spins
            # move it forward
            spins, ene[k] = McMove(spins, L, T, cases, retene=True)

        return all_spins, ene

    elif retene == False:

        for k in range(MCs):
            # record last state
            all_spins[:, :, k] = spins
            # move it forward
            spins = McMove(spins, L, T, cases)

        return all_spins


def check_mc_steps(spins_array, checks = 5):
    """
    Lets you look at the evolution of the monte carlo program at 
    evenly logarithmically spaced steps
    """
    # get steps
    steps = spins_array.shape[2]

    # set up time checks
    timechecks = geomspace(1, steps, num = checks)
    timechecks = timechecks.astype(int)

    # few adjustments
    timechecks[0] = 0
    timechecks[-1] = timechecks[-1] - 1

    # plot the results
    # since T = 1.0, they should be magnetized after 300 steps
    fig, axs = plt.subplots(ncols = checks, figsize = (16,16))
    plt.subplots_adjust(left = 0.125, right = 2.9, bottom = 0.1, top = 2.9, wspace = 0.4, hspace = 0.2)
    for i,ax in enumerate(axs):
        now = timechecks[i]
        ax.set_title("Steps = {}".format(now + 1), fontsize = 36)
        ax.imshow(spins_array[:,:,now], cmap = 'viridis')


def make_temperature_spins(MCs, L, start_temp, end_temp, num_temp, state = 0):
    """
    Makes an array of spins for various temperatures
    """
    # temp linspace
    temps = linspace(start_temp, end_temp, num = num_temp)

    # this 3d-array will contain the final configuration (2d) for a given T (+1d)
    t_spins = zeros((L, L, num_temp), dtype = int)
    for k,temp in enumerate(temps):
        t_spins[:,:,k] = ising_main(MCs, L, temp, state = state)

    return t_spins


def check_mc_temps(spin_array, start_temp, end_temp, num_temp, checks = 10):
    """
    Lets you look at how the spins vary as a function of temperature
    """
    # reproduce temps array
    temps = linspace(start_temp, end_temp, num = num_temp)

    # get idxs
    indexchecks = linspace(0, num_temp, num = checks)
    indexchecks = indexchecks.astype(int)
    indexchecks[-1] = indexchecks[-1] - 1

    # plot
    fig, axs = plt.subplots(ncols = checks, figsize = (16,16))
    plt.subplots_adjust(left = 0.125, right = 2.9, bottom = 0.1, top = 2.9, wspace = 0.4, hspace = 0.2)
    for i,ax in enumerate(axs):
        now = indexchecks[i]
        tempnow = temps[now]
        ax.set_title("T = {:3.2f}".format(tempnow), fontsize = 36)
        ax.imshow(spin_array[:,:,now], cmap = 'viridis')


def make_total_spins_energy(MCs, L, start_temp, end_temp, num_temp, state = 0):
    """
    Makes a 4d array corresponding to spins (2d), mc steps (1d) and temperature (1d)
    also returns energy as it evolves through mc steps
    """
    # make temps array
    temps = linspace(start_temp, end_temp, num = num_temp)

    # since we have to average over MC time, make 4d array
    # the last axis will have the indexes corresponding to the temp
    total_spins = zeros((L, L, MCs, num_temp), dtype = int)
    energies = zeros((MCs, num_temp), dtype = float)

    # fill it up
    for k,temp in enumerate(temps):
        total_spins[:,:,:,k], energies[:,k] = ising_main_flipbook(MCs, L, temp, retene = True, state = state)

    return total_spins, energies


def make_m_T(total_spins, start_temp, end_temp, num_temp):
    """
    Takes in 4d spins array and returns magnetication averaged over MC time (steps)
    as a function of the temperature 
    """
    # make temps array
    temps = linspace(start_temp, end_temp, num = num_temp)

    # to get m, we average over MC time and over the grid
    m_avg_over_grid = average(total_spins, axis = (0,1))

    # set up m space
    m_T = zeros((num_temp), dtype = float)

    # fill it up by averaging over steps for a given T
    for i in range(num_temp):
        # average over MCs
        m_T[i] = average(m_avg_over_grid[:,i], axis = 0)

    return m_T


def make_energy_per_spin(energies, L):
    """
    Takes in energy array and temperature array that it is based on
    Returns 1d energy per spin array with dimension of temperature array
    """
    # we will  calculate E
    # we need to average it over the time steps
    temp_size = energies.shape[1]
    E_T = zeros((temp_size), dtype = float)

    # fill it up by averaging over steps for a given T
    for i in range(temp_size):
        # average over MCs
        E_T[i] = average(energies[:,i], axis = 0)

    # finally, divide by the total number of spins
    N_spins = L ** 2
    e_T = E_T / N_spins

    return e_T


def plot_vs_temps(array, name, MCs, L, start_temp, end_temp, num_temp):
    """
    Plots array vs temperatures
    name must be in string format
    """
    # make temps
    temperatures = linspace(start_temp, end_temp, num = num_temp)

    # plot
    fig = plt.figure(figsize = (16, 8))
    plt.title(name + " as a function of temperature for {0} Monte Carlo Steps and {1} $x$ {1} grid size "
        .format(MCs, L), fontsize = 20)
    plt.scatter(temperatures, array, c = 'k', label = name)
    plt.ylim(min(array) - 0.1, max(array) + 0.1)
    plt.axhline(y = 0, ls = '--', c = 'darkcyan')
    plt.ylabel(name)
    plt.xlabel("Temperature in units of $J/k_b$")
    #plt.axvline(x = Tc, c = 'orange')
    plt.legend()
    plt.show()


def make_heat_capacity(energy, start_temp, end_temp, num_temp):
    """
    Takes energy per spin and T params and makes Cv array
    """
    # delta T
    delta_T = abs(end_temp - start_temp)
    delta_T /= num_temp

    # gradient
    Cv_T = gradient(energy, delta_T)

    return Cv_T


def critical_temperature(energies, L, start_temp, end_temp, num_temp):
    """
    Takes in energies as a function of T and MCs (2d) and temperature (1d) arrays
    Returns scalar critical temperature
    Uses the maximization of Cv as representation of critical temperature
    """
    # temps
    temps = linspace(start_temp, end_temp, num = num_temp)

    # find Cv
    eps = make_energy_per_spin(energies, L)
    heat_cap = make_heat_capacity(eps, start_temp, end_temp, num_temp)
    
    # to ignore early transients, chop off first 1/5th
    one_fifth_length = int(len(heat_cap)/5)
    heat_cap_fixed = heat_cap[one_fifth_length:]
    temps_fixed = temps[one_fifth_length:]
    
    crit_idx = argmax(heat_cap_fixed)
    crit_temp = temps_fixed[crit_idx]
    
    return crit_temp


def get_cv_exponent(Cv, Tc, start_temp, end_temp, num_temp):
    """
    Takes in Cv array and underlying temps params and returns critical exponent and error
    """
    # temps
    temps = linspace(start_temp, end_temp, num = num_temp)

    # fit the function

    # split it in two
    cv_best_params, cv_cov = opt(lambda t, a, c: tc_fit_abs(t, Tc, a ,c), xdata = temps, ydata = Cv, p0 = [0.1, 0])

    # get best params and errors
    (best_alpha, best_c) = cv_best_params
    err_params = sqrt(diag(cv_cov))
    best_alpha_err = err_params[0]

    # print out best params estimated with this method
    print("Estimated critical exponent = {:.3} +/- {:.1}".format(best_alpha, best_alpha_err))

    # make best fit array with
    best_fit = tc_fit_abs(temps, Tc, best_alpha, best_c)

    # plot results
    fig = plt.figure(figsize = (16, 8))
    plt.title("Heat capacity fit with " + "$|T_c - T|^{-alpha}$" + " and $alpha =$ {:.3} $+/-$ {:.1}".format(best_alpha, best_alpha_err))
    plt.scatter(temps, Cv, c = 'k')
    plt.plot(temps, best_fit, c = 'darkcyan', label = "Best fit")
    plt.ylabel("Heat capacity")
    plt.xlabel("Temperature in units of $J/k_b$")
    plt.legend()
    plt.show()

    return best_alpha, best_alpha_err, best_fit


def get_m_exponent(m, Tc, start_temp, end_temp, num_temp):
    """
    Takes in Cv array and underlying temps params and returns critical exponent and error
    """
    # temps
    temps = linspace(start_temp, end_temp, num = num_temp)

    # trim the function to only T < Tc
    crit_temp_idx = argmin(abs(temps - Tc))
    m_trimmed = m[:crit_temp_idx]
    temps_trimmed = temps[:crit_temp_idx]

    # fit the function
    m_best_params, m_cov = opt(lambda t, b, c, d: tc_fit(t, Tc, b, c, d), xdata = temps_trimmed, ydata = m_trimmed, p0 = [1/8, 0, 0])

    # get best params and errors
    (best_beta, best_c, best_d) = m_best_params
    err_params = sqrt(diag(m_cov))
    best_beta_err = err_params[0]

    # print out best params estimated with this method
    print("Estimated critical exponent = {:.3} +/- {:.1}". format(best_beta, best_beta_err))

    # make best fit array with
    best_fit = tc_fit(temps_trimmed, Tc, best_beta, best_c, best_d)

    # plot results
    fig = plt.figure(figsize = (16, 8))
    plt.title("Magnetization fit with " + "$(T_c - T)^{beta}$" + " and $beta =$ {:.3} $+/-$ {:.1}".format(abs(best_beta), best_beta_err))
    plt.scatter(temps_trimmed, m_trimmed, c = 'k')
    plt.plot(temps_trimmed, best_fit, c = 'darkcyan', label = "Best fit")
    plt.ylabel("Magnetization")
    plt.xlabel("Temperature in units of $J/k_b$")
    plt.legend()
    plt.show()

    return best_beta, best_beta_err, best_fit


def tc_fit_abs(temps, Tc, alpha, C):
    """
    Takes in temps array, critical temp, and critical exponent and returns (Tc - T)^-alpha
    C helps with the fit

    Unfortunately, we had to make a concession here. Namely, I manually entered the critical
    temperature so that the fit wouldn't crash by dividing by zero (because of the negative powers)
    """
    fixed = zeros(temps.shape, dtype = float)
    for i,temp in enumerate(temps):
        if temp == Tc:
            fixed[i] = 1e10 # big number
        else:
            fixed[i] = absolute(temp - Tc) ** -alpha + C
    return fixed


def tc_fit(temps, Tc, beta, C, D):
    """
    Takes in temps array, critical temp, and critical exponent and returns (Tc - T)^beta
    C and D helps with the fit
    """
    return (Tc - temps) ** beta + C * temps + D

def plot_m_and_cv(m, Cv, Tc, start_temp, end_temp, num_temp):
    """
    Plots magnetization per spin and heat capacity on the same plot to compare
    """
    # make temps space
    temps = linspace(start_temp, end_temp, num = num_temp) 

    # plot
    fig = plt.figure(figsize = (16, 8))
    plt.title("Magnetization and Heat Capacity as a function of Temperature")
    plt.scatter(temps, m, c = 'maroon', label = "$m$ data")
    plt.scatter(temps, Cv, c = 'darkcyan', label = "$C_V$ data")
    plt.axvline(x = Tc, c = 'orange', ls = '--')
    plt.legend()
    plt.show()






