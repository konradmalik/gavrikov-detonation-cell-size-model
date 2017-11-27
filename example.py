import SDToolbox as sd
import numpy as np
from effective_activation_energy import ea_r
from znd_model import znd_CJ
from gavrikov import gavrikov_model

##############################################################################

Tinit = 300.0
Pinit = sd.one_atm
NITROGEN_OXYGEN_RATIO = 79.0/21.0

# phases
mech = 'gri30_highT.cti'

# gaseous fuel species
fuel_specie = 'H2'

# oxidant
oxidant = 'air' # oxygen, air

# fuel fraction
fuel_vol = 0.3

##############################################################################

if (oxidant == 'air'):
    o2_vol = (1 - fuel_vol)/(NITROGEN_OXYGEN_RATIO + 1.0)
    n2_vol = 1 - fuel_vol - o2_vol
else:
    o2_vol = 1 - fuel_vol
    n2_vol = 0.0

q = '{0}:{1:0.9g}, O2:{2:0.4g}, N2:{3:0.4g}'.format(fuel_specie, fuel_vol, o2_vol, n2_vol)
print q

# calculate Ea/R, postshock T
[Ea_R, Tps] = ea_r(Pinit, Tinit, q, mech)
print 'Ea_R is: {}, Tps is: {} K'.format(Ea_R, Tps)

# calculate ZND detonation to get 0.75 Ma reaction length [m]
[_, [maxM, _, Tvn, _, reaction_length, _]] = znd_CJ(Pinit, Tinit, q, mech)
print 'Maximum Mach was {}, Reaction length is: {} m, Tvn is {}'.format(maxM,reaction_length,Tvn)

# gavrikov
cell = 1000.0*gavrikov_model(reaction_length, Ea_R, Tps, Tvn, Tinit) # in mm
print 'Cell size was {} mm'.format(cell)
