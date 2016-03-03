#A class that describes an event - this is simply a data object, it has no functions of its own
class Event:

    def __init__(self):
        self.waiting_time=float('inf') #The default waiting time is infinity so any real event will be faster
        self.first_site=-1 #The first site that the event occurs on
        self.second_site=-1 #The second site that the event will occur on - for hops or 2 body interactions we need both sites
        self.event_type="o" #This will describe the type of event
        self.particle_ID=[None]*2 #The ID(s) of the particle(s) involved.

