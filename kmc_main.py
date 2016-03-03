import kmc_model_class as KMC
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("Mode",help="Mode to run in: (S)teps, T(ime) or (C)arriers")
parser.add_argument("Duration",help="The duration of the simulaiton in the previous mode",type=int)
parser.add_argument("Net_bias",help="The total reversed bias across the device (V)",type=float)
args=parser.parse_args()
mode=args.Mode
duration=args.Duration
bias=args.Net_bias

Model=KMC.KMC(300,bias,"DEFAULT")
if(mode=="S"):
    current_density=Model.run_kmc_steps(duration)
elif(mode=="T"):
    current_density=Model.run_kmc_time(duration)
elif(mode=="C"):
    current_density=Model.run_kmc_carriers(duration)
else:
    print "Mode not recognised"
    exit(101)
print 'Simulation complete J=',current_density,'A/m2'
