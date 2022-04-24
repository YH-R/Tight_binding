import os
import numpy as np
from math import exp
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


# configure color map to be transparent when white
# get colormap
ncolors = 2
color_array = plt.get_cmap('gray')(range(ncolors))
color_array[:,-1] = np.linspace(1.0,0.0,ncolors) # change alpha values
map_object = LinearSegmentedColormap.from_list(name='gray_alpha',colors=color_array) # create a colormap object
plt.register_cmap(cmap=map_object) # register this new colormap with matplotlib


# physical constants (in SI units)
e_over_kb = 11603
e2_over_hbar = 0.00024349


# check MATLAB for values used in calculations

# t: hopping coefficient (eV)
t = 3;

# mu: chemical potential (eV)
# mu_list: range of mu to investigate
mu_steps = 120
mu_step = 2*4*t / mu_steps # max = 4t, min = -4t
# mu_list = [x * mu_step - 4*t for x in range(0, mu_steps+1)]
mu_list = [-3]

# beta: inverse of temperature times k_b(eV^-1)
# beta_list: range of beta to investigate
T_steps = 60
T_low = 1 # K
T_high = 80000
T_step = (T_high - T_low)/T_steps
#beta_list = [e_over_kb/(T_low + x*T_step) for x in range(T_steps)] # eV^-1
beta_list=[39]
#beta_list = [3**(1-x/15) for x in range(T_steps)]

# B: magnetic field (T), update manually every run
B=1

# list of csv files, must be in the same directory
# usually use one file at a time
# use multiple files when investigating size-dependent effects
csv_list=['current_30x30_0.00011577_hexagon.csv']


def main(beta_list=beta_list, mu_list=mu_list, t=t, B=B, csv_list=csv_list, triangle_mode=False):
    out_counter = 0 # to establish ordering of output files, to make gif
    for current_csv in csv_list:
        print("Using " + current_csv)
        
        # file input
        (Nx, Ny, energies, current_data) = read_currents(current_csv)
        
        for beta in beta_list:
            print("Inverse temperature: ", beta)
            
            for mu in mu_list:
                # vectorized fermi function
                vectorized_f = np.vectorize(lambda E: f(E, beta, mu), otypes=[np.float])

                (current_list, start_indices_list, end_indices_list) = calculate_average_current(vectorized_f, energies, current_data)
                

                # for scaling when plotting
                maximum_I = max(map(abs,current_list)) # default choice
                # maximum_I = 2.4e-07 # manual normalization
                print("Current normalization is " + str(maximum_I) + "A.")

                # array to store currents
                edges = np.zeros((Ny * 2 - 1, Nx * 2 - 1))
                    
                atoms = np.fromfunction(np.vectorize(is_atom), (2*Ny-1,2*Nx-1)) # black squares for lattice sites

                for count in range(len(current_list)):
                    i, j = map_edge_to_edges(start_indices_list[count], end_indices_list[count], Ny)
                    edges[j, i] = current_list[count]


                out_counter += 1
                # plots the current and optionally saves it
                plotting(beta, t, mu, B, Nx, Ny, edges, atoms, maximum_I, out_counter, savefig=False, showfig=True)


# E and mu: eV
# beta: eV^-1
def f(E, beta, mu):
    '''Fermi function'''
    # print(E)
    if beta*(E-mu) < 99: # to avoid overflow
        return 1/(exp(beta*(E-mu))+1)
    else:
        return 0


def read_currents(filename, verbose=0):
    '''read in spreadsheet and store data in numpy arrays'''

    # 2D numpy array   
    data = np.genfromtxt(filename, delimiter=',')

    Nx = int(max(data[:,-2]))
    Ny = int(max(data[:,-1]))

    # first row is energy eigenvalues
    energies = data[0][:-4]
    #print("energies")
    #print(energies)

    # remove first row, the rest are expectation values in different eigenstates
    current_data = data[1:]

    if verbose > 0:
        print("energies")
        print(energies)
    if verbose > 1:
        print("current_data")
        print(current_data)

    return (Nx, Ny, energies, current_data)


# current expectation using Fermi-Dirac distribution as weights
def calculate_average_current(vectorized_f, energies, current_data, verbose=0):
    fermi_weights = vectorized_f(energies) # 1D array
    #print("fermi_weights")
    #print(fermi_weights)

    current_list = []
    start_indices_list = []
    end_indices_list = []

    for row in current_data:
        currents = row[:-4] * e2_over_hbar # convert to Ampere
        start_indices = row[-4:-2].astype(int) # (x, y)
        end_indices = row[-2:].astype(int)

        # sum over expectation values of currents in energy eigenstates,
        # weighted by fermi function of energy eigenvalue
        average_current = np.sum(currents*fermi_weights)

        current_list.append(average_current)
        start_indices_list.append(start_indices)
        end_indices_list.append(end_indices)

    if verbose > 0:
        print("current_list")
        print(current_list)
    if verbose > 1:
        print("start_indices_list")
        print(start_indices_list)
        print("end_indices_list")
        print(end_indices_list)

    return (current_list, start_indices_list, end_indices_list)
    

def is_atom(i, j):
    '''returns 1 if the coordinates (i,j) represents an atom, otherwise 0.'''
    # the returned value is used for plotting subsequently
    
    if i%2==0 and j%2==0:
        return 0 # atoms is black
    else:
        return 1 # otherwise white (transparent)


def edge_type(start_indices, end_indices):
    '''returns string description of edge type, can be either HORIZONTAL, VERTICAL, DIAGONAL or ERROR'''
    # indices are of the form (x,y)
    # assume start_indices are smaller than end indices (currents going upwards or to the right)
    if start_indices[1] == end_indices[1] and end_indices[0] - start_indices[0]==1:
        # same y coordinates
        return "HORIZONTAL"
    elif start_indices[0] == end_indices[0] and end_indices[1] - start_indices[1]==1:
        # same x coordinates
        return "VERTICAL"
    elif end_indices[0] - start_indices[0] == 1 and end_indices[1] - start_indices[1]==-1:
        # going down and right
        return "DIAGONAL"
    else:
        print("Oh no! Error in edge_type.")
        print("start_indices: " + str(start_indices))
        print("end_indices: " + str(end_indices))
        
        return "ERROR"
    

def map_edge_to_edges(start_indices, end_indices, Ny):
    '''returns the coordinate of the edges on the map'''
    orientation = edge_type(start_indices, end_indices)

    i = 2*start_indices[0] - 2
    j = 2*Ny - 2*start_indices[1] # python convention is going down

    if orientation == "HORIZONTAL":
        i += 1
    if orientation == "VERTICAL":
        j -= 1
    if orientation == "DIAGONAL":
        i += 1
        j += 1

    return (i, j)

def plotting(beta, t, mu, B, Nx, Ny, edges, atoms, maximum_I, out_counter, savefig=False, showfig=True):    
    plt.imshow(edges, cmap='bwr') # blue positive, red negative

    if maximum_I>0:
        plt.clim(-maximum_I, maximum_I)
    else:
        plt.clim(-1, 1)


    temp = e_over_kb / beta # in K, for labelling plot
    # print("Temperature is " + "{:.2f}".format(temp) + "K. ")


    plt.subplots_adjust(bottom=0.35) # make space at bottom for text

    # text labels
    plt.text(x=0.37, y=0.3, s=r"$N_x$: " + str(Nx) + r"$ ,N_y$: " + str(Ny), fontsize=10, transform=plt.gcf().transFigure)
    plt.text(x=0.37, y=0.25, s=r"$T$: " + "{:.2f}".format(temp) + "K", fontsize=10, transform=plt.gcf().transFigure)
    plt.text(x=0.37, y=0.2, s=r"$B$: " + "{:.4g}".format(B) + "T", fontsize=10, transform=plt.gcf().transFigure)
    plt.text(x=0.37, y=0.15, s=r"$\mu$: " + "{:.4f}".format(mu) + "eV", fontsize=10, transform=plt.gcf().transFigure)
    plt.text(x=0.37, y=0.1, s=r"$I$ normalization: " + "{0:.4g}".format(maximum_I) + "A", fontsize=10, transform=plt.gcf().transFigure)

    plt.axis('off')

    img_atoms = plt.imshow(atoms, cmap='gray_alpha')

    if savefig:
        out_dir = "D:/Cloud/OneDrive/OneDrive - National University of Singapore/Notes/NUS/Physics/FYP/Matlab/PlotCurrent_out/" # check of directory exists
        out_name = str(out_counter) + "_beta="+"{0:.4g}".format(beta)+"_t="+str(t)+"_mu="+"{:+.4f}".format(mu)+"_B="+"{:+.4g}".format(B)+"_Nx="+str(Nx)+"_Ny="+str(Ny)+".png"
        out_path = os.path.join(out_dir, out_name)
        plt.savefig(out_path, bbox_inches='tight')
        print("Successfully saved " + out_path)

    if showfig:
        plt.show()

    plt.clf()


# main gaurd
if __name__ == '__main__':
    main()
