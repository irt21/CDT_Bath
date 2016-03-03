import math
#A file of constants and useful functions

lx=60
ly=30
lz=30
lattice_points=lx*ly*lz
lattice_size=[lx,ly,lz]
lattice_spacing=2e-9 #2 nanometres
lattice_spacing_nm=2 #2 nanometres
cross_section=ly*lz*(lattice_spacing**2) #Cross sectional area of the terminal

fundamental_charge=1.602e-19
eV=fundamental_charge
kb_J=1.380649e-23 #in J per Kelvin
kb=kb_J/eV #kb in eV
pi=3.14159
epsilon_0=8.854e-12 #Farads per metre
Polaric_energy=0.187#eV 2.99574e-20 J
localisation=2
Forster_radius=2.1e-9 #m
coulomb_constant=4*pi*epsilon_0*3.5

HOMO_energy_n=-6.1
LUMO_energy_n=-3.7
HOMO_energy_p=-5.0
LUMO_energy_p=-2.8

sigma_n_type=0.01 #The standard deviation of the Gaussian disorder in the n-type material in eV
sigma_p_type=0.01 #The standard deviation of the Gaussian disorder in the p-type material in eV

#We are going to model traps as a site that has the same energy for HOMO and LUMO
traps_on=True #Toggle to True to put traps in
#The density of defects in the n-type material--This is the number of defects per lattice site so <1
trap_density_n_type=0.01
trap_energy_n_type=-5.5 #The energy of defects in the n-type material in eV
trap_density_p_type=0.01
trap_energy_p_type=-5.0 #The energy of defects in the p-type material in eV

mobility_e=1e-7
mobility_h=1e-7

#All rates and times are in picoseconds - see Marcus rate prefactors
recombination_rate=1e-4
dissociation_rate=10 #1 over lifetime in picoseconds
singlet_decay_rate=2e-3# 1 over lifetime = 54ns
#Lets assume all photons are 460nm and 1 sun=1367W per m^2
#generation_rate=3.157e9 #photons per ps per m^2
photon_rate=3.157e-9 #for nm^2
absorption_factor=(4*pi*0.5)/460 #absorption factor per nm

#A function that returns the index of particle with given ID, only one particle has each ID so the loop quits once found
def find_particle_by_ID(list,ID):
    for i in range(len(list)):
        if(list[i].ID==ID):
            return i
    print "Particle ID "+ID+" not found"
    return math.inf

#A function that finds the index of a particle by its location
def find_particle_by_location(list,location):
    for i in range(len(list)):
        if(list[i].location==location):
            return i
    print "Particle at location ",location," not found"
    return math.inf

#A functio nthat returns the distance between two lattice site, periodic in y and z
def get_distance(site_1,site_2):
    coordinates_1=coordinates_from_index(site_1)
    coordinates_2=coordinates_from_index(site_2)
    separation_sq=0.0
    for dim in xrange(3):
        delta=coordinates_1[dim]-coordinates_2[dim]
        if(dim>0 and delta>0.5*lattice_size[dim]):
            delta=delta-lattice_size[dim]
        if(dim>0 and delta<-0.5*lattice_size[dim]):
            delta+=lattice_size[dim]
        separation_sq+=delta*delta
    return math.sqrt(separation_sq)*lattice_spacing

#A function that returns a list of coordinates from a lattice site index
def coordinates_from_index(index):
    coordinates=[None]*3
    coordinates[2]=index//(lattice_size[0]*lattice_size[1])
    coordinates[1]=(index-(coordinates[2]*lattice_size[0]*lattice_size[1]))/lattice_size[0]
    coordinates[0]=(index-(coordinates[2]*lattice_size[0]*lattice_size[1]))%lattice_size[0]
    return coordinates

#A function that returns a lattice site index from a list of coordinates
def index_from_coordinates(coordinates):
    if(coordinates[1]<0):
        coordinates[1]+=ly
    if(coordinates[1]>=ly):
        coordinates[1]-=ly
    if(coordinates[2]<0):
        coordinates[2]+=lz
    if(coordinates[2]>=lz):
        coordinates[2]-=lz
    return ((coordinates[2]*lattice_size[0]*lattice_size[1])+(coordinates[1]*lattice_size[0])+coordinates[0])


