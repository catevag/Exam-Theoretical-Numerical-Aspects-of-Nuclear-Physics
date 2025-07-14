

#%%starting



import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from scipy.special import erf
from scipy.optimize import minimize
import matplotlib.cm as cm  # for colormap


#create ouput directories if not present
for folder in ["potentials", "eigenfunctions", "energy", "phase_shift", "relative_errors"]:os.makedirs(folder, exist_ok=True)
  

#%% choice of the global parameters
r_max = 2000 #must be >100
x_min = -8. #corresponds to r_min == 3 * 1E-4 Bohr radii
dx = 0.01 #grid spacing 
n_max = 20 #maximum energy level (number of nodes limited by range)
a = 1 #cutoff for renormalization theory
a_values = [0.7, 1, 2, 5, 10] #to compare results varying the cutoff, in the last part of the program

#%%setting the coulomb potential

def coulomb_initialization(r):

    #definition of the potential
    v_potential = -1/r #Z=1, m=1

    #saving the potential to CSV using pandas df
    df_pot = pd.DataFrame({"r": r, "V(r)": v_potential})
    df_pot.to_csv("potentials/coulomb_potential.csv", index=False)
    print("potential saved in 'potentials/coulomb_potential.csv'\n")

    
    #plotting the potential
 
    plt.figure()
    plt.title("Coulomb potential")
    plt.xlabel("r")
    plt.ylabel("V(r)")
    plt.xlim(0, 0.01)  # radius range
    plt.grid(True)

    plt.plot(r, v_potential, color='black', label=r'$V(r)$, atomic units')
    plt.legend()

    plt.savefig("potentials/coulomb_potential_plot.pdf")
    plt.close()
    
    return v_potential


#%% setting the real potential

def short_range_initialization(r):

    #definition of the potential
    #v_potential = -1/r - ( np.exp(- r)) / r 
    v_potential = -1/r - 1 / (1 + np.exp(r))

    #saving the potential to CSV
    df_pot = pd.DataFrame({"r": r, "V(r)": v_potential})
    df_pot.to_csv("potentials/real_potential.csv", index=False)
    print("potential saved in 'potentials/real_potential.csv'\n")

    
    #plotting the potential
    plt.figure()
    plt.title("Real potential with short range term")
    plt.xlabel("r")
    plt.ylabel("V(r)")
    plt.xlim(0, 0.01)
    plt.grid(True)

    plt.plot(r, v_potential, color='green', label=r'$V(r)$, atomic units')
    plt.legend()
    
    plt.savefig("potentials/real_potential_plot.pdf")
    plt.close()

    return v_potential

#%%setting effective potential

def effective_initialization(r, a, c, d_1):

    #definition of the potential 
    v_potential = (
    -1 / r * erf(r / (np.sqrt(2) * a))
    + c * np.power(a, 2) * np.exp(-np.power(r, 2) / (2 * np.power(a, 2)))
      / (np.power(2 * np.pi, 1.5) * np.power(a, 3))
    + d_1 * np.power(a, 4) * ((np.power(r, 2) / np.power(a, 2)) - 3)
      * np.exp(-np.power(r, 2) / (2 * np.power(a, 2)))
      / (np.power(2 * np.pi, 1.5) * np.power(a, 5))
    )
    

    #saving the potential to CSV
    df_pot = pd.DataFrame({"r": r, "V(r)": v_potential})
    if d_1 == 0:
        df_pot.to_csv("potentials/effective_potential_a2.csv", index=False) #
    else:
        df_pot.to_csv("potentials/effective_potential_a4.csv", index=False) #
    

    #plotting the potential
    plt.figure()

    plt.title("Effective Potential")
    plt.xlabel("r")
    plt.ylabel("V(r)")
    plt.grid(True)
    plt.xlim(0, 10)
    #plt.show()
    plt.plot(r, v_potential, color='green', label=r'V(r), atomic units')
    plt.legend()

    if d_1 == 0:
        plt.savefig("potentials/effective_potential_a2_plot.pdf")
    else:
        plt.savefig("potentials/effective_potential_a4_plot.pdf")
    

    plt.close()


    return v_potential

#%%schrodinger


def solve_schrodinger(n_grid, dx, v_potential, r2, r, sqrt_r, n, l, i_100):
    #solve the radial schrodinger equation on a logarithmic
    #grid by Numerov method

    eps = 1E-10 #tolerance for eigenvalue
    n_iter = 200 #max iterations

    #useful quantities
    ddx12 = (dx * dx) / 12
    lpar = (l + 0.5) * (l + 0.5) 
    
    #initial lower and upper bounds to the eigenvalue
    e_up = v_potential[n_grid] #bound state has negative energy
    e_low = e_up

    #min energy is greater than the min of the effective potential
    for j in range(0, n_grid + 1):
            e_low = np.minimum(e_low, lpar / (2 * r2[j] ) + v_potential[j])

    if (e_up - e_low < eps):
        print("error in solving schrodinger: e_up and e_low coincide", file = sys.stderr)
        sys.exit(1)
            
    e = (e_low + e_up) * 0.5 #first rough estimate
  

    g = np.zeros(n_grid + 1) # Numerov g
    
    de = 1E10 #any number greater than epsilon to start the cycle
    
    inv_point = -1 #index of classical inversion point
    
    #start loop to find energy eigenvalue
    i = 0
    while i  < n_iter and np.absolute(de) > eps:
        #set up the f-function (in a way to determine the position of its last change of sign)

        g[0] = ddx12 * ((2 * r2[0] * (v_potential[0] - e)) + lpar)
        for j in range(1, n_grid + 1):
            g[j] = ddx12 * ((2 * r2[j] * (v_potential[j] - e)) + lpar)

            #if g[j] is zero the change of sign is not observed 
            #trick to prevent missing change of sign
            if (g[j] == 0.):
                g[j] = 1E-10
                
            # g > 0 approximately means classically forbidden region
            # g < 0  allowed
            #take the index of classical inversion
            if np.sign(g[j]) != np.sign(g[j - 1]):
                inv_point = j
                
            
                

        if inv_point < 0 or inv_point >= n_grid - 2:
           
            print(f"inversion = {inv_point:4d}, {n_grid:4d}, n = {n}, e ={e}")
            print("error in solving schrodinger: last change of sign too far", file = sys.stderr)
            plt.figure()
            
            sys.exit(1)

        #rewrite the f-function how required by numerov method
        f = 1 - g

        y = np.zeros(n_grid + 1) #wavefunction initialization
        nodes = n - l - 1 #analytical
        
        #wavefunction in the first two points
        y[0] = 1E-12
        y[1] = 1E-5

        #outward integration with node counting
        n_cross = 0
        for j in range(1, inv_point):
            #numerov formula
            y[j + 1] = ((12. - f[j] * 10.) * y[j] - (f[j - 1] * y[j - 1])) / f[j + 1]
            if np.sign(y[j]) != np.sign(y[j + 1]):
                n_cross += 1
            if (y[j + 1] > 1E10):
                for m in range(1, j + 2, +1): #prevent overflow
                    y[m] = y[m] / y[j + 1]

        scale_factor = y[inv_point] #value of the wavefunction at classical turning point, to match outward and inward
        
        #check the number of crossings
        if (n_cross != nodes):
            #incorrect number of nodes, adjusting eigenvalue
            if (n_cross > nodes): #means my eigenvalue is too high for n
                e_up = e
            else:
                e_low = e

            e = (e_up + e_low) * 0.5 
            

        else:
            #correct number of nodes, we can perform inward integration.

            #determination of the wavefunction in last two points
            #assuming
            
            y[n_grid] = 0
            y[n_grid - 1] = dx

            #inward integration
            for j in range(n_grid - 1, inv_point, -1):
                y[j - 1] = ((12. - f[j] * 10.) * y[j] - (f[j + 1] * y[j + 1])) / f[j - 1]
                if (y[j - 1] > 1E10):
                    for m in range(n_grid, j - 2, -1): #prevent overflow
                        y[m] = y[m] / y[j - 1]

            #rescale the function to match at the classical turning point
            scale_factor /= y[inv_point] #ratio between outward and inward
            for j in range(inv_point, n_grid + 1):
                y[j] *= scale_factor


            #normalize wavefunction
            norm = 0.

            for j in range(0, n_grid + 1):
                norm += y[j] * y[j] * r2[j] * dx #approx integration

            norm = np.sqrt(norm)

            for j in range(0, n_grid + 1):
                y[j] /= norm
         



            

            #improving convergence with perturbation theory
            
            #find the value of the cusp at the matching point
            j = inv_point
            y_cusp = (y[j - 1] * f[j - 1] + f[j + 1] * y[j + 1] + f[j] * 10. * y[j]) / 12.
            #ycusp is the value predicted by the Numerovs method using xcl as central point
            #fcusp is value of f if consider delta term in potential
            #df = fcusp - f(inv_point)
            df_cusp = f[j] * ((y[j] / y_cusp) - 1.)

            # eigenvalue update using perturbation theory
            de = df_cusp / ddx12 * y[j] * y[j] * dx / 2
            if (de > 0.):
                e_low = e
            if (de < 0.):
                e_up = e
            
            if np.abs(de) > (e_up - e_low) :
                # bisection step
                e = 0.5 * (e_low + e_up)
            else:
                # perturbative step
                e = e + de
           
            #prevent e to go out of bounds ( e > e_up or e < e_low)
            #could happen far from convergence
            e = np.minimum(e, e_up)
            e = np.maximum(e, e_low)
         

        i += 1
    
    #convergence not achived
    if (np.abs(de) > eps): #i reached n_iter
        if n_cross != nodes:
            print(f"n_cross={n_cross:4d} nodes={nodes:4d} inv_point={inv_point:4d} " 
                f"e={e:16.8e} e_low={e_low:16.8e} e_up={e_up:18.8e}", file=sys.stderr)
        else:
            print(f"e={e:16.8e} de={de:16.8e}", file=sys.stderr)

        #print (n)
        
        print(f"solve_schrodinger not converged after {n_iter} iterations", file=sys.stderr)
       
        
        sys.exit(1)

    
    #compute phase shift ar r = 100
    
    phase_shift = np.arcsin(y[i_100]) - np.sqrt(np.absolute(2 * e)) * r[i_100] + (l * np.pi)/2
    #phase_shift = np.mod(phase_shift, 2 * np.pi)

    return e, y, phase_shift

#%%analysis


def analysis(n_max, pot_name, e_n, p_n, y_n, n_grid, dx, v_potential, r2 , r, sqrt_r, i_100):
    #quantic numbers
    n=1 #increases inside the cycle
    l=0 #fixed, considering only s-wave

    for i in range(1, n_max +1):
        e_n[i-1], y_n[i-1], p_n[i-1] = solve_schrodinger(n_grid, dx, v_potential, r2, r, sqrt_r, n ,l, i_100)

        #saving and plotting eigenfunctions
        y_df = pd.DataFrame({"r": r, "y": y_n[i-1]})
        y_df.to_csv(f"eigenfunctions/eigenfunction_{pot_name}_{i}.csv", index=False)


        plt.figure()

        plt.title(f"{pot_name} wavefunction n = {i}")
        plt.xlabel("r")
        plt.ylabel(r'$\psi(r)$')
        plt.grid(True)
        
        plt.plot(r, sqrt_r * y_n[i-1], color='red', label=r'$\psi$') #go back to phys wf
        plt.legend()
        
        plt.savefig(f"eigenfunctions/{pot_name}_wavefunction_{i}.pdf")
        plt.close()

        
        n += 1 #next energy level and the cycle restarts (i++)

    n_array = np.arange(1, n_max + 1) #1 2 ... nmax.

    #saving eigenvalues, phase shifts
    e_df = pd.DataFrame({"n": n_array, "E": e_n})
    #print("-----------------------------------------------\n")
    #print(f"{pot_name} energy eigenvalues\n")
    #print(e_df.to_markdown(index=False)) #table wo indexes
    e_df.to_csv(f"energy/energy_eigenvalues_{pot_name}.csv", index=False)

    p_df = pd.DataFrame({"n": n_array, "phase shift": p_n})
    #print("-----------------------------------------------\n")
    #print(f"{pot_name} phase shifts at r = 100\n")
    #print(p_df.to_markdown(index=False))
    p_df.to_csv(f"phase_shift/phase_shift_{pot_name}.csv", index=False)
    
    df_combined = pd.DataFrame({
    "n": n_array,
    "E": e_n,
    "phase shift": p_n
    })

    print("-----------------------------------------------\n")
    print(f"{pot_name} energy eigenvalues and phase shifts \n")
    print(df_combined.to_markdown(index=False))

    return e_n, y_n, p_n







#%%grid preparation

def grid_preparation(n_grid, x_min, dx):

    #preparing x-array with constant step
    x = np.linspace(x_min, x_min + ((n_grid - 1) * dx), n_grid)
    x = np.append(x, x[n_grid - 1] + dx) #another point for convenient labelling

    #generate r, sqrt_r, and r^2 (logarithmic grid)
    r = np.exp(x)
    sqrt_r = np.sqrt(r)
    r2 = np.power(r, 2)
    for i in range(n_grid):
        if r[i]>100:
            i_100 = i  #index corresponding to r = 100
            break

    #print grid information
    print("Radial grid:\n")
    print("dx = ", dx)
    print("x_min = ", x_min)
    print("number of points = ", n_grid)
    print("r_min = ", r[0])
    print("r_max = ", r[n_grid])
    print("-----------------------------------------------\n")

    return r, sqrt_r, r2, i_100


#%%cost functions for effective theory tuning

#cost function for a2-theory
def cost_function_1(c, a, n_grid, dx, r2, r, sqrt_r, n_max):
    v_potential = effective_initialization(r, a, c, 0) #imposing d_1 = 0
    e, y, phase_shift = solve_schrodinger(n_grid, dx, v_potential, r2, r, sqrt_r, n_max, 0, i_100)
    return np.absolute(ph_shift_real[-1] - phase_shift) #index -1 for last element
    
#cost function for a4-theory
def cost_function_2(parameters, a, n_grid, dx, r2, r, sqrt_r, n_max):
    c = parameters[0]
    d_1 = parameters[1]
    v_potential = effective_initialization(r, a, c, d_1)
    e, y, phase_shift = solve_schrodinger(n_grid, dx, v_potential, r2, r, sqrt_r, n_max, 0, i_100)
    
    return np.absolute(ph_shift_real[-1] - phase_shift) #phase_shift at n_max(=10)
#%%make grid


#number of points of the grid
n_grid = int((np.log(r_max) - x_min) / dx)

#initialize logarithmic grid (globally)
r, sqrt_r, r2, i_100 = grid_preparation(n_grid, x_min, dx)


#%%coulomb potential


#initialize potential
v_potential_1 = coulomb_initialization(r)
v_potential_1_name = "coulomb"

#initialize arrays for eigenvalues, eigenfunctions and phase shifts
e_coulomb = np.zeros(n_max)
y_coulomb = np.zeros((n_max, n_grid + 1))
ph_shift_coulomb = np.zeros(n_max)

e_coulomb, y_coulomb, ph_shift_coulomb = analysis(
    n_max, v_potential_1_name, e_coulomb, ph_shift_coulomb, y_coulomb,
    n_grid, dx, v_potential_1, r2, r, sqrt_r, i_100
)

#comparison with analytical values -1/n^2
e_theory = np.arange(1, n_max + 1)
e_theory = -1 / (2 * e_theory ** 2)

plt.figure()

plt.title("Relative Error (vs Analytical Values) Coulomb Energies")
plt.xlabel(r'$|E|$')
plt.ylabel(r'$|\Delta E / E|$')
plt.grid(True)

plt.plot(
    np.absolute(e_theory),
    np.absolute(e_coulomb - e_theory) / np.absolute(e_theory),
    marker='o',
    linestyle='--',
    color='black',
    label="relative error coulomb"
)

plt.xscale("log")
plt.yscale("log")
plt.legend()

plt.savefig("relative_errors/relative_error_coulomb_eigenvalues.pdf")
plt.show()
plt.close()



#%%short range potential


v_potential_2 = short_range_initialization(r)
v_potential_2_name = "real"

#initialize arrays for eigenvalues, eigenfunctions and phase shifts
e_real = np.zeros(n_max)
y_real = np.zeros((n_max, n_grid + 1))
ph_shift_real = np.zeros(n_max)

e_real, y_real, ph_shift_real = analysis(
    n_max, v_potential_2_name, e_real, ph_shift_real, y_real,
    n_grid, dx, v_potential_2, r2, r, sqrt_r, i_100
)

c = np.sqrt(np.pi) *  np.power(n_max, 3) * ( e_real[n_max - 1] + 1 / (2 * np.power(n_max, 2) ) )
print("\nc = ", c)


#%%relative errors perturbative delta


#coulomb eigenvalues in e_theo

#approximate eigenvalues
e_pert = np.arange(1,n_max + 1) #1, 2, 3, 4...
#formula 5 Lepage
e_pert = - 1/ ( 2 * np.power(e_pert, 2) ) + c / (np.sqrt(np.pi) * np.power(e_pert, 3))

#plotting relative errors
plt.figure()

plt.title("Relative Error for Energies - Comparison Coulomb vs $\delta$")
plt.xlabel(r'$|E|$')
plt.ylabel(r'$|\Delta E / E|$')
plt.grid(True)
plt.xscale("log")
plt.yscale("log")

plt.plot(
    np.absolute(e_real),
    np.absolute(e_real - e_theory) / np.absolute(e_real),
    marker='o',
    linestyle='--',
    color='black',
    label="Coulomb"
)

plt.plot(
    np.absolute(e_real),
    np.absolute(e_real - e_pert) / np.absolute(e_real),
    marker='o',
    linestyle='--',
    color='blue',
    label="$\delta$ - Perturbative"
)

plt.legend()
plt.savefig("relative_errors/relative_error_coulomb_vs_delta.pdf")
plt.show()
plt.close()







#%% a^2 theory


#tuning of the parameter (lowest energy data - phase shift)
c_a2 = -40 # initial guess of the parameter a2



#compute parameter minimizing the phase shift

result1 = minimize(
    cost_function_1, 
    c_a2, 
    method='nelder-mead', #accurate
    args=(a, n_grid, dx, r2, r, sqrt_r, n_max),
    options={'maxiter': 10000, 'maxfev': 10000, 'xatol': 1E-8, 'disp': False}
)



c_a2 = result1.x
    

    
#print(f"-------------------------------------------------\n")
print(f"(a^2 theory):    c = {c_a2} \n")


#initialize potential
v_potential_3 = effective_initialization(r, a, c_a2, 0)
v_potential_3_name = "effective_a^2"

plt.figure()
plt.title("Effective Potential")
plt.xlabel("r")
plt.ylabel("V(r)")
plt.grid(True)
plt.xlim(0, 10)

plt.plot(r, v_potential_3, color='green', label=r'V(r), atomic units')
plt.legend()
plt.show()

plt.close()
#initialize arrays for eigenvalues, eigenfunctions and phase shifts
e_eff_a2 = np.zeros(n_max)
y_eff_a2 = np.zeros((n_max, n_grid + 1))
ph_shift_eff_a2 = np.zeros(n_max)



#analysis
e_eff_a2, y_eff_a2, ph_shift_eff_a2 = analysis(
n_max, v_potential_3_name, e_eff_a2, ph_shift_eff_a2, y_eff_a2,
n_grid, dx, v_potential_3, r2, r, sqrt_r, i_100
)

#%% a^4 theory

#initial guess of parameters
init_parameters = [-40, -1] #[c, d_1] a4

result2 = minimize(
    cost_function_2, 
    init_parameters, 
    method='nelder-mead',
    args=(a, n_grid, dx, r2, r, sqrt_r, n_max),
    options={'maxiter': 10000, 'maxfev': 10000, 'xatol': 1E-8, 'disp': False}
)

parameters_a4 = result2.x
c_a4 = parameters_a4[0]
d_1 = parameters_a4[1]

print(f"(a^4 theory): c = {c_a4}, d_1 = {d_1}\n")

v_potential_4 = effective_initialization(r, a, c_a4, d_1)
v_potential_4_name = "effective_a^4"

plt.figure()
plt.title("Effective Potential")
plt.xlabel("r")
plt.ylabel("V(r)")
plt.grid(True)
plt.xlim(0, 10)
#plt.show()
plt.plot(r, v_potential_4, color='green', label=r'V(r), atomic units')
plt.legend()

plt.close()


e_eff_a4 = np.zeros(n_max)
y_eff_a4 = np.zeros((n_max, n_grid + 1))
ph_shift_eff_a4 = np.zeros(n_max)

e_eff_a4, y_eff_a4, ph_shift_eff_a4 = analysis(
n_max, v_potential_4_name, e_eff_a4, ph_shift_eff_a4, y_eff_a4,
n_grid, dx, v_potential_4, r2, r, sqrt_r, i_100
)

    
    
#%%final plots

# --- Plot 0 : a^2 vs a^4
plt.figure()
plt.plot(
    np.abs(e_real),
    np.abs(e_real - e_eff_a4) / np.abs(e_real),
    marker='o',
    linestyle='--',
    label=r'$a^4$ effective theory',
    color='orange'
    )
plt.plot(
    np.abs(e_real),
    np.abs(e_real - e_eff_a2) / np.abs(e_real),
    marker='o',
    linestyle='--',
    label=r'$a^2$ effective theory',
    color='purple'
    )

plt.xlabel(r'$|E|$')
plt.ylabel(r'$|\Delta E / E|$')
plt.xscale("log")
plt.yscale("log")
plt.title("effective theory Comparison")
plt.legend()
plt.grid(True)
plt.savefig("relative_errors/relative_errors_comparison_a_orders.pdf")
plt.show()
plt.close()


# --- Plot 1: Potentials Comparison ---
plt.figure()
#plt.plot(r, v_potential_1, label=r'$V_{coulomb}(r)$', color='black')
plt.plot(r, v_potential_2, label=r'$V_{real}(r)$', color='green')
plt.plot(r, v_potential_4, label=r'$V_{a^4}(r)$', color='orange')

plt.xlabel("r")
plt.ylabel("V(r)")
plt.title("Potentials - Comparison")
plt.grid(True)
plt.legend()
plt.xlim(0, 0.01)
plt.savefig("potentials/potential_comparison.pdf")
#plt.show()
plt.close()


# --- Plot 2: Relative Error on Energies ---
plt.figure()
plt.plot(
    np.abs(e_real),
    np.abs(e_real - e_theory) / np.abs(e_real),
    marker='o',
    linestyle='--',
    label='Coulomb',
    color='black'
)
plt.plot(
    np.abs(e_real),
    np.abs(e_real - e_pert) / np.abs(e_real),
    marker='o',
    linestyle='--',
    label=r'$\delta$ - Perturbation',
    color='blue'
)
plt.plot(
    np.abs(e_real),
    np.abs(e_real - e_eff_a4) / np.abs(e_real),
    marker='o',
    linestyle='--',
    label=r'$a^4$ Theory',
    color='orange'
)

plt.xlabel(r'$|E|$')
plt.ylabel(r'$|\Delta E / E|$')
plt.xscale("log")
plt.yscale("log")
plt.title("Relative Error on Energy Eigenvalues - Comparison")
plt.legend()
plt.grid(True)
plt.savefig("relative_errors/relative_error_energy_comparison.pdf")
plt.show()
plt.close()


# --- Plot 3: Relative Error on Phase Shift ---
plt.figure()
plt.plot(
    np.abs(e_real),
    np.abs(ph_shift_real - ph_shift_coulomb) / np.abs(ph_shift_real),
    marker='o',
    linestyle='--',
    label='Coulomb',
    color='black'
)
plt.plot(
    np.abs(e_real),
    np.abs(ph_shift_real - ph_shift_eff_a4) / np.abs(ph_shift_real),
    marker='o',
    linestyle='--',
    label=r'$a^4$ Theory',
    color='orange'
)

plt.xlabel(r'$|E|$')
plt.ylabel(r'$|\Delta \varphi / \varphi|$')
plt.xscale("log")
plt.yscale("log")
plt.title("Relative Error on Phase Shift - Comparison")
plt.legend()
plt.grid(True)
plt.savefig("relative_errors/relative_error_phase_shift_comparison.pdf")
plt.show()
plt.close()


#%% comparison of different a values



n_checks_a = len(a_values)

# Tuning of the parameter (lowest energy data - phase shift)
c_a2_comp = -40  # initial guess

# Initialize arrays outside loop
e_eff_a2_comp = np.zeros((n_max, n_checks_a))
y_eff_a2_comp = np.zeros((n_max, n_grid + 1, n_checks_a))
ph_shift_eff_a2_comp = np.zeros((n_max, n_checks_a))

# Temporary arrays for results
e_eff_a2_temp = np.zeros(n_max)
y_eff_a2_temp = np.zeros((n_max, n_grid + 1))
ph_shift_eff_a2_temp = np.zeros(n_max)

for j in range(n_checks_a):
    
    # use data already computed
    if (a_values[j] == a):
        v_potential_3_temp = v_potential_3
        v_potential_3_name_temp = f"effective_a^2_{a_values[j]}"
    else:
    # Compute parameter minimizing the phase shift
        result1 = minimize(
            cost_function_1,
            c_a2_comp,
            method='nelder-mead',
            args=(a_values[j], n_grid, dx, r2, r, sqrt_r, n_max),
            options={'maxiter': 10000, 'maxfev': 10000, 'xatol': 1E-8, 'disp': False}
            )
        c_a2_temp = result1.x
        
        print(f"(a^2 theory), a = {a_values[j]}:    c = {c_a2_temp} \n")
        
        # Initialize potential
        v_potential_3_temp = effective_initialization(r, a_values[j], c_a2, 0)
        v_potential_3_name_temp = f"effective_a^2_{a_values[j]}"
        
        
        
        # Analysis
        e_eff_a2_temp, y_eff_a2_temp, ph_shift_eff_a2_temp = analysis(
            n_max, v_potential_3_name_temp, e_eff_a2_temp, ph_shift_eff_a2_temp, y_eff_a2_temp,
            n_grid, dx, v_potential_3_temp, r2, r, sqrt_r, i_100
            )
        
    print("\n ----------------------------------------------- \n")
        
    # Store results
    e_eff_a2_comp[:, j] = e_eff_a2_temp
    y_eff_a2_comp[:, :, j] = y_eff_a2_temp
    ph_shift_eff_a2_comp[:, j] = ph_shift_eff_a2_temp


#%% a values plot

# Plot all relative errors in one plot 
plt.figure()
colors = cm.Purples(np.linspace(0.3, 0.9, n_checks_a))  # generate visually distinct colors

for j in range(n_checks_a):
    plt.plot(
        np.abs(e_real),
        np.abs(e_real - e_eff_a2_comp[:, j]) / np.abs(e_real),
        marker='o',
        linestyle='--',
        label=fr'$a={a_values[j]}$',
        color=colors[j]
    )

plt.xlabel(r'$|E|$')
plt.ylabel(r'$|\Delta E / E|$')
plt.xscale("log")
plt.yscale("log")
plt.title("Relative Error for Different Values of a, a^2 Theory")
plt.legend()
plt.grid(True)
plt.savefig("relative_errors/relative_error_energy_comparison_a.pdf")
plt.show()
plt.close()

