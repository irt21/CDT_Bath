#A class that defines a mobile particle
class Particle:

    def __init__(self,particle_type,location,ID,timestamp): #The constructor that creates a particle object of given type at given location
        self.particle_type=particle_type
        if(particle_type=="e"): #If we make an electron
            self.charge=-1.602e-19
            self.hopping_range=1
        elif(particle_type=="h"): #If we make a hole
            self.charge=1.602e-19
            self.hopping_range=1
        elif(particle_type=="s"): #If we make a singlet
            self.charge=0
            self.hopping_range=5
        else:
            print "CREATING A PARTICLE OF UNKNOWN TYPE" #Error message 
            exit(45) #Quit the program with error code 45
        self.location=location #N.B. the particle does not update the lattice, that is up to the user
        self.ID=ID #The unique ID of the particle - not its location in the particle list because particles are added and deleted
        self.creation_time=timestamp #The creation time of the particle
