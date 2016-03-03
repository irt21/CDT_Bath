#A class that describes a single site of the lattice
class Site:

    def __init__(self,HOMO_energy,LUMO_energy,site_type,ID):
        self.HOMO_energy=HOMO_energy #The HOMO energy of the site (in eV)
        self.LUMO_energy=LUMO_energy #The LUMO energy of the site iin eV)
        self.site_type=site_type #The type of site - at the moment it will be n or p but could be more
        self.ID=ID #The ID of the lattice, also its poisiotn in the lattice list
        self.occupied=False #True or False depending on if the site is occupied by a particle
