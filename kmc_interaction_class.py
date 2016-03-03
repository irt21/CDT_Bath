#A class that calculates the possible events, rates and times
import kmc_constants as const
import kmc_event_class as EVENT
import kmc_particle_class as PARTICLE
import kmc_site_class as SITE
import math
import random

class Interaction_Handler:
#The constructor for the interaction handler class it needs to know the temperature
    def __init__(self,temperature):
        self.T=temperature
        mobility_prefactor=(6*const.kb_J*self.T*const.mobility_e)/(const.fundamental_charge*const.lattice_spacing*const.lattice_spacing) #This is the mobility prefctor for e-
        #The e-12 is to got picosends in the rate prefactor
        self.hop_e=(mobility_prefactor/(math.exp(-const.Polaric_energy/(4*const.kb*self.T))*math.exp(-2*const.localisation*const.lattice_spacing)))*1e-12 #Marcus rate prefactor for e- hops
        mobility_prefactor=(6*const.kb_J*self.T*const.mobility_h)/(const.fundamental_charge*const.lattice_spacing*const.lattice_spacing) #Mobility prefactor for h+
        self.hop_h=(mobility_prefactor/(math.exp(-const.Polaric_energy/(4*const.kb*self.T))*math.exp(-2*const.localisation*const.lattice_spacing)))*1e-12 #Marcus rate prefactor for h+ hops

#This function is long but quite simple, it checks each particle, where it can hop and if there is an allowed event associated with that hop
    def get_fastest_event(self,lattice,particlelist):
        fastest_event=EVENT.Event() #Create an empty event to store our fastest event
    #Make a generation event
        new_site=random.randint(0,const.lattice_points-1)
    #Uniform absorption moving down along x
        generation_rate=const.photon_rate*const.lx*(const.lattice_spacing_nm**3)*const.absorption_factor #Absorbtion in a column of the device - the coefficients are per nm
        waiting_time=(-1.0*math.log(random.random()))/(generation_rate)
        if(waiting_time<fastest_event.waiting_time and lattice[new_site].occupied==False): #If this is the fastest event
            fastest_event.first_site=new_site
            fastest_event.waiting_time=waiting_time
            fastest_event.event_type="generation"
        for particle in particlelist: #Loop through all mobile particles in the system
            current_site=particle.location #The current site
            old_electrostatic_energy=self.get_coulomb_energy(particle,current_site,particlelist) #The current electrostatic energy
            if(particle.particle_type=="e"):
                old_site_energy=lattice[current_site].LUMO_energy #Site energy for e- that live in LUMO
            if(particle.particle_type=="h"):
                old_site_energy=lattice[current_site].HOMO_energy #h+ live in HOMO
            #A 3 long list to hold x,y,z positions in lattice
            coordinates=const.coordinates_from_index(current_site)
            #loop through the nearby sites and check new possible lattice sites
            for dz in range(-particle.hopping_range,particle.hopping_range+1): #dz is the change in the z direction
                for dy in range(-particle.hopping_range,particle.hopping_range+1):
                    for dx in range(-particle.hopping_range,particle.hopping_range+1):
                        distance=math.sqrt((dx*dx)+(dy*dy)+(dz*dz))*const.lattice_spacing_nm #The distance between the suggested and current sites in nm
                        new_site=const.index_from_coordinates([coordinates[0]+dx,coordinates[1]+dy,coordinates[2]+dz]) #The site ID of the suggested site
                        if(new_site!=current_site and (coordinates[0]+dx)>=0 and (coordinates[0]+dx)<const.lx):
                            if(particle.particle_type=="e"):
                                new_site_energy=lattice[new_site].LUMO_energy
                            if(particle.particle_type=="h"):
                                new_site_energy=lattice[new_site].HOMO_energy
                            #Hopping event for a charged particle
                            if(lattice[new_site].occupied==False #Particles can only hop to an empty site
                            and ((lattice[new_site].site_type=="n" and particle.particle_type=="e") #e- can only hop to n-type
                            or (lattice[new_site].site_type=="p" and particle.particle_type=="h"))):#h+ can only hop to p-type
                                new_electrostatic_energy=self.get_coulomb_energy(particle,new_site,particlelist)
                                delta_energy=(new_electrostatic_energy+new_site_energy)-(old_electrostatic_energy+old_site_energy)
                                #Calculate hopping rate  - Marcus rates (see lecture notes)
                                if(particle.particle_type=="e"):
                                    rate=self.hop_e*math.exp(-2*const.localisation*distance)
                                    rate=rate*math.exp(-1*( ((delta_energy+const.Polaric_energy)**2)/(4*const.Polaric_energy*const.kb*self.T)))
                                elif(particle.particle_type=="h"):
                                    rate=self.hop_h*math.exp(-2*const.localisation*distance)
                                    rate=rate*math.exp(-1*(((-delta_energy+const.Polaric_energy)**2)/(4*const.Polaric_energy*const.kb*self.T)))
                                waiting_time=(-1.0*math.log(random.random()))/rate
                                if(waiting_time<fastest_event.waiting_time): #If this possible event is the fastest save it
                                    fastest_event.first_site=current_site #The old site
                                    fastest_event.second_site=new_site #The new location of the particle
                                    fastest_event.waiting_time=waiting_time #Event waiting time
                                    fastest_event.event_type="move" #The type of move
                                    fastest_event.particle_ID[0]=particle.ID #Which particle is involved
                            #Hopping event for a singlet Forster process -- compare to the charge hop, same code but a different rate
                            if(lattice[new_site].occupied==False and particle.particle_type=="s"):
                                rate=const.singlet_decay_rate*((const.Forster_radius/distance)**6)
                                waiting_time=(-1.0*math.log(random.random()))/const.dissociation_rate
                                if(waiting_time<fastest_event.waiting_time):
                                    fastest_event.first_site=current_site
                                    fastest_event.second_site=new_site
                                    fastest_event.waiting_time=waiting_time
                                    fastest_event.event_type="move"
                                    fastest_event.particle_ID[0]=particle.ID
                            #Singlet decay event -- compare to the hops, again similar idea but a different rate and there is no second site involved
                            if(particle.particle_type=="s"):
                                waiting_time=(-1.0*math.log(random.random()))/const.singlet_decay_rate
                                if(waiting_time<fastest_event.waiting_time):
                                    fastest_event.first_site=current_site
                                    fastest_event.event_type="decay"
                                    fastest_event.particle_ID[0]=particle.ID
                                    fastest_event.waiting_time=waiting_time
                            #Dissociation event on p type at boundary, a singlet in the p-type material near a boundary splits into a hole at the same location and an electron on the n-type side
                            if(lattice[new_site].occupied==False and particle.particle_type=="s"
                               and lattice[new_site].site_type=="n" and lattice[current_site].site_type=="p"):
                                waiting_time=(-1.0*math.log(random.random()))/const.dissociation_rate
                                if(waiting_time<fastest_event.waiting_time):
                                    fastest_event.particle_ID[0]=particle.ID
                                    fastest_event.first_site=current_site
                                    fastest_event.second_site=new_site
                                    fastest_event.event_type="dissociate"
                                    fastest_event.waiting_time=waiting_time
                            #Dissociation event on n type at boundary, same as the p-type dissociation but the singlet is in n-type
                            if(lattice[new_site].occupied==False and particle.particle_type=="s"
                               and lattice[new_site].site_type=="p" and lattice[current_site].site_type=="n"):
                                waiting_time=(-1.0*math.log(random.random()))/const.dissociation_rate
                                if(waiting_time<fastest_event.waiting_time):
                                    fastest_event.particle_ID[0]=particle.ID
                                    fastest_event.first_site=current_site
                                    fastest_event.second_site=new_site
                                    fastest_event.event_type="dissociate"
                                    fastest_event.waiting_time=waiting_time
                            #Recombination event, if a charge carrier can hop onto another of the opposite species they can recombine and emit a photon
                            if(lattice[new_site].occupied==True and (particle.particle_type=="e" or particle.particle_type=="h")):
                                particle_2=particlelist[const.find_particle_by_location(particlelist,new_site)]
                                if((particle.particle_type=="e" and particle_2.particle_type=="h") or
                                (particle.particle_type=="h" and particle_2.particle_type=="e")):
                                    waiting_time=(-1.0*math.log(random.random()))/const.recombination_rate
                                    if(waiting_time<fastest_event.waiting_time):
                                        fastest_event.first_site=current_site
                                        fastest_event.second_site=new_site
                                        fastest_event.event_type="recombine"
                                        fastest_event.waiting_time=waiting_time
                                        fastest_event.particle_ID[0]=particle.ID
                                        fastest_event.particle_ID[1]=particle_2.ID
        return fastest_event #Give the fastest event back to the Model class to execute the event

    def get_coulomb_energy(self,particle,site,particlelist): #Calculate the Coulomb energy of the given particle on the given site 
        coulomb_energy=0.0
        for particle_2 in particlelist: #Loop through all the particles
            if(particle.ID!=particle_2.ID and site!=particle_2.location):#The particle cannot interact with itself at the old site
                separation=const.get_distance(site,particle_2.location) #Distance between particles in m
                coulomb_energy+=(particle.charge*particle_2.charge)/(const.coulomb_constant*separation) #The electrostatic energy between the given particle and the looped particle
        return coulomb_energy/const.eV #We are using eV to describe the energies in the exponentials
