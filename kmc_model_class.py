#A class that holds the entire lattice, the list of particles and the interaction handler - essentially the entire model
import kmc_site_class as SITE
import kmc_particle_class as PARTICLE
import kmc_event_class as EVENT
import kmc_constants as const
import kmc_interaction_class as INTERACTIONS
import math
import random

class KMC:
#The constructor for the KMC class, it takes three extra arguments
    def __init__(self,temperature,bias,morphology_file_name):
        self.lattice=[]
        self.interactions=INTERACTIONS.Interaction_Handler(temperature) #An object that will calculate all possible events
        self.particlelist=[] #An empty list that will hold the particle information
        self.simulation_time=0.0
        self.recombinations=0
        self.generations=0
        self.dissociations=0
        self.singlet_dissociation_time=[]
        self.electrons_out=0
        self.electron_collection_time=[]
        self.holes_out=0
        self.hole_collection_time=[]
        self.num_particles=0
        self.newest_ID=0
        self.field=bias/(const.lx) #The field gradient per lattice site in the x direction
        self.current_density=0
        self.load_morphology(morphology_file_name) #Load the morphology for the system

#Run a simulation limited by the numebr of Monte Carlo steps
    def run_kmc_steps(self,num_steps):
        for i in xrange(num_steps):     #Run the simulation until we have performed num-steps Monte Carlo steps
            self.kmc_step()
            if(i%500==0):   #Periodically write information to the terminal for the user to see
                print i, self.simulation_time/1e12, self.current_density
        self.output_results()
        return self.current_density

#Run a simulation limited by the amount of simulation time in seconds
    def run_kmc_time(self,target_time):
        steps=0
        while (self.simulation_time<target_time*1e12): #Run the simultion until the simulation time exceeds the target time (N.B. s->ps)
            self.kmc_step()
            steps+=1
            if(steps%500==0):   #Periodically write information to the terminal for the user to see
                print steps, self.simulation_time/1e12, self.current_density
        self.output_results()
        return self.current_density

#Run a simulation limited by the number of electrons out
    def run_kmc_carriers(self,carriers_out):
        steps=0
        while (self.electrons_out<carriers_out):    #Run the simulation until 
            self.kmc_step()
            steps+=1
            if(steps%500==0):   #Periodically write information to the terminal for the user to see
                print steps, self.simulation_time/1e12, self.current_density
        self.output_results()
        return self.current_density

#load morphology and apply the bias to the energy levels
    def load_morphology(self,morphology_file_name):
    #Load a simple bilayer
        if(morphology_file_name=="DEFAULT"):
            for z in xrange(const.lz):  #Loop through all lattice sites
                for y in xrange(const.ly):
                    for x in xrange(const.lx):
                        #Add a site to the lattice with the properties of n or p type -- Bilayer here is n from x=0 -> x=L/2
                        if(x<const.lx/2):
                            (self.lattice.append(SITE.Site(const.HOMO_energy_n,const.LUMO_energy_n,
                                                           "n",const.index_from_coordinates([x,y,z])))) #Add a n-type site
                        else:
                            (self.lattice.append(SITE.Site(const.HOMO_energy_p,const.LUMO_energy_p,
                                                           "p",const.index_from_coordinates([x,y,z])))) #Add a p-type site
        else:  #Load a file with the given file name
            counter=0
            with open(morphology_file_name,'r') as file:
                for line in file: #The file format has to be one line per site and in columns HOMO LUMO type ID
                    line=line.split()
                    list.append(SITE.Site(list[0],list[1],list[3],list[4])) #Add the site to the lattice
                    counter+=1
            file.closed
            if(counter!=lattice_points):
                print("Error in the number of input lines")     #Error check that we have the right number of sites
                exit(35) #Exit with error code 35

        for i in xrange(const.lattice_points): #Loop through all sites
            coordinates=const.coordinates_from_index(i)
            self.lattice[i].HOMO_energy+=(self.field*coordinates[0])  #Apply the electric field to the site energy levels
            self.lattice[i].LUMO_energy+=(self.field*coordinates[0])
   #Apply Gaussian disorder to the energy levels, there is the magnitude of the disorder and the variation
            if(self.lattice[i].site_type=="n"): #n-type material
                self.lattice[i].HOMO_energy+=random.gauss(0,const.sigma_n_type)
                self.lattice[i].LUMO_energy+=random.gauss(0,const.sigma_n_type)
            elif(self.lattice[i].site_type=="p"): #p-type material
                self.lattice[i].HOMO_energy+=random.gauss(0,const.sigma_p_type)
                self.lattice[i].LUMO_energy+=random.gauss(0,const.sigma_p_type)
        if(const.traps_on==True):
            self.apply_traps()

#Add traps to the system at random
    def apply_traps(self):
        for i in xrange(const.lattice_points):
            if(self.lattice[i].site_type=="n"): #n-type material
                if(random.random()<const.trap_density_n_type):
                    self.lattice[i].HOMO_energy=const.trap_energy_n_type
                    self.lattice[i].LUMO_energy=const.trap_energy_n_type
            if(self.lattice[i].site_type=="p"): #n-type material
                if(random.random()<const.trap_density_p_type):
                    self.lattice[i].HOMO_energy=const.trap_energy_p_type
                    self.lattice[i].LUMO_energy=const.trap_energy_p_type

#KMC_Step function
    def kmc_step(self):
        Fastest_event=self.interactions.get_fastest_event(self.lattice,self.particlelist) #Get the fastest event
        self.simulation_time=self.simulation_time+Fastest_event.waiting_time #Increment the simulation time by the event's waiting time
        particle_collection_time=self.event_handler(Fastest_event) #Execute the fastest event
        net_charge=(self.electrons_out+self.holes_out)*const.fundamental_charge
        self.current_density=net_charge/((self.simulation_time/1e12)*const.cross_section)#Calcualte current density in SI units

#Output results to files 
    def output_results(self):
        with open("output_V%1.1f.csv"%self.field,'w') as output:
            s='Simulation time \t %f ps\t%f s\n'%(self.simulation_time,self.simulation_time/1e12)
            output.write(s)
            s='Internal field of %f per lattice site\n'%(self.field)
            output.write(s)
            s='Net bias of %f\n'%(self.field*const.lx)
            output.write(s)
            s='electrons out \t %f\n'%(self.electrons_out)
            output.write(s)
            s='holes out \t %f \n'%(self.holes_out)
            output.write(s)
            s='dissociations \t %f\n'%(self.dissociations)
            output.write(s)
            s='recombinations \t%f\n'%(self.recombinations)
            output.write(s)
            s='generations \t %f\n'%(self.generations)
            output.write(s)
            output.write('Current density\t%f\n'%(self.current_density))
        output.closed
        with open("HOMO-LUMO.csv",'w') as HL:
            for i in xrange(const.lx):
                HL.write('%g\t%g\t%g\n'%(i,self.lattice[i].HOMO_energy,self.lattice[i].LUMO_energy))
        HL.closed
        with open("electron_collection_time.csv",'w') as e_coll:
            for i in xrange(len(self.electron_collection_time)):
                e_coll.write('%f\n'%(self.electron_collection_time[i]))
        e_coll.closed
        with open("hole_collection_time.csv",'w') as h_coll:
            for i in xrange(len(self.hole_collection_time)):
                h_coll.write('%f\n'%(self.hole_collection_time[i]))
        h_coll.closed
        with open("singlet_dissociation_time.csv",'w') as s_coll:
            for i in xrange(len(self.singlet_dissociation_time)):
                s_coll.write('%f\n'%(self.singlet_dissociation_time[i]))
        s_coll.closed

#Perform the event - this consists of updating the particle's location, the particle list and the lattice occupancy
    def event_handler(self,Fastest_event):
        particle_collection_time=-1 #A variable to hold particle lifetimes if any are removed
#Move a particle
        if(Fastest_event.event_type=="move"):
            self.lattice[Fastest_event.first_site].occupied=False #Where the particle used to be is now empty
            self.lattice[Fastest_event.second_site].occupied=True #New site is occupied
            particle_index=const.find_particle_by_ID(self.particlelist,Fastest_event.particle_ID[0]) #Find the location of the particle in the list
            self.particlelist[particle_index].location=Fastest_event.second_site #Go to the correct particle in the list and update its location
            coordinates=const.coordinates_from_index(Fastest_event.second_site) #Find the 3D coordinates of the particle from the site ID
            #Check for extraction event
            if(self.particlelist[particle_index].particle_type=="e" and coordinates[0]==0): #If an electron is on the electrode
                particle_collection_time=self.simulation_time-self.particlelist[particle_index].creation_time
                self.electron_collection_time.append(particle_collection_time)
                self.electrons_out+=1
                self.num_particles-=1
                self.lattice[Fastest_event.second_site].occupied=False #Remove the particle from the lattice site
                del self.particlelist[particle_index] #Delete the particle from the list of particles as its been extracted
            elif(self.particlelist[particle_index].particle_type=="h" and coordinates[0]==const.lx-1): #If a hole is on the anode
                particle_collection_time=self.simulation_time-self.particlelist[particle_index].creation_time
                self.hole_collection_time.append(particle_collection_time)
                self.holes_out+=1
                self.num_particles-=1
                self.lattice[Fastest_event.second_site].occupied=False #Remove the particle from the lattice site
                del self.particlelist[particle_index] #Delete the particle from the list of particles as its been extracted
#Generate a singlet
        if(Fastest_event.event_type=="generation"):
            self.generations+=1
            self.lattice[Fastest_event.first_site].occupied=True #Add the singlet to the lattice
            self.particlelist.append(PARTICLE.Particle("s",Fastest_event.first_site,self.newest_ID,self.simulation_time)) #Add the singlet to the list of particles
            self.newest_ID+=1
            self.num_particles+=1
#Do a singlet dissociation on n-type side of boundary
        if(Fastest_event.event_type=="dissociate" and self.lattice[Fastest_event.first_site].site_type=="n"):
            singlet_index=const.find_particle_by_ID(self.particlelist,Fastest_event.particle_ID[0]) #Find where the singlet is in the particle list
            del self.particlelist[singlet_index] #Delete the singlet as it has dissociated - N.B. we could technically combine this line with the one above
            self.dissociations+=1
            #self.lattice[Fastest_event.first_site].occupied=True #Technically unnecessary as it was already true from the singlet
            self.lattice[Fastest_event.second_site].occupied=True #The second lattice site must be updated with the hole
            self.particlelist.append(PARTICLE.Particle("e",Fastest_event.first_site,self.newest_ID,self.simulation_time)) #Add the electron to the particle list
            self.newest_ID+=1
            self.particlelist.append(PARTICLE.Particle("h",Fastest_event.second_site,self.newest_ID,self.simulation_time)) #Add the hole to the particle list
            self.newest_ID+=1
            self.num_particles+=2
#Do a singlet dissociation on p-type side of boundary
        if(Fastest_event.event_type=="dissociate" and self.lattice[Fastest_event.first_site].site_type=="p"):
            singlet_ID=const.find_particle_by_ID(self.particlelist,Fastest_event.particle_ID[0])
            particle_diss_time=self.simulation_time-self.particlelist[singlet_ID].creation_time
            self.singlet_dissociation_time.append(particle_diss_time)
            del self.particlelist[singlet_ID]
            self.dissociations+=1
            self.lattice[Fastest_event.first_site].occupied=True
            self.lattice[Fastest_event.second_site].occupied=True
            self.particlelist.append(PARTICLE.Particle("h",Fastest_event.first_site,self.newest_ID,self.simulation_time))
            self.newest_ID=self.newest_ID+1
            self.particlelist.append(PARTICLE.Particle("e",Fastest_event.second_site,self.newest_ID,self.simulation_time))
            self.newest_ID=self.newest_ID+1
            self.num_particles+=2
#Do a singlet decay, it recombines and emits a photon
        if(Fastest_event.event_type=="decay"):
            self.lattice[Fastest_event.first_site].occupied=False #Nothing remains on the lattice site
            singlet_index=const.find_particle_by_ID(self.particlelist,Fastest_event.particle_ID[0]) #Find where the singlet is in the list of particles
            del self.particlelist[singlet_index] #Delete the singlet from the list of particles
            self.num_particles-=1
#Do a recombination event, an electron and a hole recombine and emit a photon
        if(Fastest_event.event_type=="recombine"):
            self.recombinations+=1
            self.lattice[Fastest_event.first_site].occupied=False #Nothing remains on the lattice site
            self.lattice[Fastest_event.second_site].occupied=False #Ditto from above
            particle_ID_1=const.find_particle_by_ID(self.particlelist,Fastest_event.particle_ID[0]) #Find the location of the first particle in the list of particles
            particle_1_lifetime=self.simulation_time-self.particlelist[particle_ID_1].creation_time
            del self.particlelist[particle_ID_1] #Delete the particles from the list
            particle_ID_2=const.find_particle_by_ID(self.particlelist,Fastest_event.particle_ID[1]) #Find the location of the second particle in the list of particles
            particle_2_lifetime=self.simulation_time-self.particlelist[particle_ID_2].creation_time
            del self.particlelist[particle_ID_2] #Delete the particle from the list
            self.num_particles-=2
        return particle_collection_time
