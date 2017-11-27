import SDToolbox as sd
import numpy as np
from scipy.integrate import odeint
import time
import matplotlib.pyplot as plt

def znd_shk(U1, P1, T1, q, mech):
    
    gas1 = sd.Solution(mech);
    
    gas1.TPX = T1, P1, q
    
    gas = sd.PostShock.PostShock_fr(U1, P1, T1, q, mech);
    
    out = znd_detonation(gas,gas1,U1);
    
    return out
    
def znd_CJ(P1, T1, q, mech):

    plt_num = 1
    cj_speed, _ = sd.PostShock.CJspeed(P1, T1, q, mech, plt_num);
        
    outCJ = znd_shk(cj_speed, P1, T1, q, mech)
    
    return [cj_speed, outCJ]

'''
%Set of ODE's to solve ZND Detonation Problem
%
% FUNCTION
% SYNTAX
% dydt = ZNDReactor(t,y,gas,U1,r1,PSC)
% ALWAYS called by and ODE solver, in this case: out = ode15s(@ZNDReactor,tel,y0,options,gas,U1,r1,PSC)
%       tel = time span
%       y0 = initial conditions array
%       options = specified in znd_detonation.m
%       gas = Cantera gas object
%       U1 = Shock Velocity
%       r1 = initial density
%       PSC = Post-shock pressure
%
% INPUT
% t = time
% y = solution array
% gas = Cantera gas object
% U1 = Shock Velocity
% r1 = initial density
% PSC = Post-shock pressure
%
% OUTPUT
% dydt = Array of ODEs to be solved by ode15s
%
% SUBFUNCTION CALLS
% Cantera Functions: set.m, meanMolecularWeight.m, gasconstant.m,
%       density.m, nSpecies.m, netProdRates.m, enthalpies_RT.m,
%       molecularWeights.m, cp_mass.m, soundspeed.m,  
'''

def ZNDReactor(y,t,gas,U1,r1,PSC):

    gas.TDY = gas.T, y[1], y[3:]
    
    wt = gas.mean_molecular_weight;

    rho = gas.density;
    T = (y[0]*PSC/y[1])*(wt/sd.gas_constant);

    gas.TDY = T, rho, y[3:]
    nsp = gas.n_species;
            
    # Vectors
    wdot = gas.net_production_rates;

    hs = gas.standard_enthalpies_RT*sd.gas_constant*T;
    mw = gas.molecular_weights;
        
    # Scalars
    cp = gas.cp_mass;

    c = soundSpeedM(gas);
    
    U = U1*r1/rho;
    M = U/c;                       #Mach Number
    eta = 1 - M**2;                 #Sonic Parameter

    # % Loop through all of the species,calculate thermicity for this time step
    dykdt = np.array([0]*nsp, dtype=np.float);
    sum = 0;
    for z in range(0,nsp):
        dykdt[z] = mw[z]*wdot[z]/rho; #Net production rate of species z (kg/sec)
        drdy = -wt/mw[z];             # mean molecular weight / molecular weight of species z
        a = hs[z]/(cp*T*mw[z]);       # factor that accounts for chemical energy change
        sum = sum - (drdy + a)*dykdt[z];
        
    sigmadot = sum;                 #Thermicity

    Pdot = -rho*U**2*sum/eta/PSC;   #Pressure Derivative (Normalized by Post-Shock pressure to reduce size)
    rdot = -rho*sum/eta;           #Density Derivative
    #xdot = U;                      %Distance Derivative
    '''
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % SET UP COLUMN VECTOR FOR DYDT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    '''
    list3 = np.array([Pdot, rdot, U]);
    list4 = np.array([0]*nsp, dtype=np.float)
    dydt = np.concatenate((list3, list4));
        
    # species equations
        
    for i in range(0,nsp):
        dydt[i+3]=mw[i]*wdot[i]/rho;
        
    return dydt

def znd_detonation(gas,gas1,U1):

    r1 = gas1.density;
    
    P1 = gas1.P;
    T1 = gas1.T;
    PSC = gas.P;
    rnew = gas.density;
    nsp = gas.n_species;

    list1 = np.array([1, rnew, 0]);
    list2 = np.array(gas.Y);

    y0 = np.concatenate((list1, list2));


    # ujemna temperatura/gestosc znaczy zmniejsz 2 liczbe, bo za dlugo sie liczy  
    # blad w int(cjInd) znaczy zwieksz 2 liczbe bo za krotko sie liczy i nie dochodzimy do M=1  
    tel = np.linspace(0, 0.000005, num=100000);

    out = odeint(ZNDReactor, y0, tel, args=(gas,U1,r1,PSC),
                                full_output = 0, col_deriv = False,
                                printmessg = 1, rtol = 1.e-5, atol = 1.e-8)
    #print out
    b,a = out.shape

    i = 1;
    j = 1;
    
    y = np.array([0]*nsp, dtype=np.float) 

    #Initalize ouput matrices
    outputInd_len_ZND = 0;
    outputInd_time_ZND = 0;
    outputExo_len_ZND = 0;
    outputExo_time_ZND = 0;
    outputTime = np.array([0]*b, dtype=np.float) ;
    outputDistance = np.array([0]*b, dtype=np.float) ;
    outputP = np.array([0]*b, dtype=np.float) ;
    outputT = np.array([0]*b, dtype=np.float) ;
    outputU = np.array([0]*b, dtype=np.float) ;
    outputRho = np.array([0]*b, dtype=np.float) ;
    outputThermicity = np.array([0]*b, dtype=np.float) ;
    outputM = np.array([0]*b, dtype=np.float) ;
    outputAf = np.array([0]*b, dtype=np.float) ;
    outputG = np.array([0]*b, dtype=np.float) ;
    outputWt = np.array([0]*b, dtype=np.float) ;
    outputSonic = np.array([0]*b, dtype=np.float) ;
    outputSpecies = np.zeros(shape=(b,nsp), dtype=np.float);
    temp_grad = np.array([0]*b, dtype=np.float);

    # Extract TIME, TEMPERATURE, PRESSURE, and DENSITY arrays from integrator output
    for i in range(0,b):
        
        outputTime[i] = tel[i];

        for m in range(3,a):
            y[m-3] = abs(out[i,m]);
            
        gas.TDY = gas.T, out[i,1], y;

        wt = gas.mean_molecular_weight;

        Tout = (out[i,0]*PSC/out[i,1])*(wt/sd.gas_constant);
        
        gas.TPY = Tout, out[i,0]*PSC, y;
            
        outputP[i] = out[i,0]*PSC/sd.one_atm;
        outputRho[i] = out[i,1];
        outputDistance[i] = out[i,2];
        outputT[i] = Tout;
        for k in range(0,nsp):
            outputSpecies[i,k] = y[k];

        # Extract WEIGHT, GAMMA, SOUND SPEED, VELOCITY, MACH NUMBER, c^2-U^2,
        # THERMICITY, and TEMPERATURE GRADIENT
 
        den = outputRho[i];
        nsp = gas.n_species;

        #Vectors
        wdot = gas.net_production_rates;
        hs = gas.standard_enthalpies_RT*sd.gas_constant*Tout;
        mw = gas.molecular_weights;

        #Scalars
        cp = gas.cp_mass;
        cv = gas.cv_mass;
        g = cp/cv;
        af = np.sqrt(g*sd.gas_constant/wt*Tout); #% soundspeed(gas);
        r = gas.density;
            
        U = U1*r1/r;
        M = U/af;                               #Mach Number
        eta = 1 - M**2;                          #Sonic Parameter
        sonic = eta*af**2;

        sum = 0;
        for n in range(0,nsp):
            h = hs[n]/mw[n];
            wd = wdot[n];
            w = mw[n];
            dykdt = w*wd/den;
            drdy = -den*wt/w;
            term = den*h/(cp*Tout);
            sum = sum - (drdy + term)*dykdt;
            
        sigma = sum/den;                        #Thermicity
        Pdot = -U**2*sum/eta/PSC;               #Pressure Derivative
        rdot = -sum/eta;                        #Density Derivative

        # FIND TEMPERATURE GRADIENT
        temp_grad[i] = Tout*(Pdot/out[i,0] - rdot/r);

            
        # Assign ouptput structure
        outputU[i] = U;
        outputThermicity[i] = sigma;
        outputM[i] = M;
        outputAf[i] = af;
        outputG[i] = g;
        outputWt[i] = wt;
        outputSonic[i] = sonic;

    #%Find INDUCTION TIME and LENGTH based on MAXIMUM TEMPERATURE GRADIENT and Maximum
    #%Thermicity
        
    # Procdure to find max value with 
    maxTgrad = np.max(temp_grad)
    maxThermicity = np.max(outputThermicity)
    k1 = np.where(temp_grad == maxTgrad)
    n = np.where(outputThermicity == maxThermicity)

    #print 'reaction zone according to MaxGradTemp'
    rZoneGradTemp = outputDistance[k1]
    #print rZoneGradTemp
    
    outputInd_time_ZND = outputTime[n];
    #print 'reaction zone according to MaxThermicity'
    outputInd_len_ZND = outputDistance[n];
    #print outputInd_len_ZND
    
    maxT = np.max(outputT)
    k = np.where(outputT == maxT)
    # maxT,k = outputT[:].max(0),outputT[:].argmax(0)

    #print 'reaction zone according to MaxTemp'
    rZoneMaxTemp = outputDistance[k]
    #print rZoneMaxTemp

    #############################################################
    #%Find reaction LENGTH based on CJ point M = 0.9
    Ind = np.where((outputM <= 0.9001) & (outputM >= 0.8999))
    Ind = int(np.mean(Ind))
    
    reactionLen09 = outputDistance[Ind]
    #%Find reaction LENGTH based on CJ point M = 0.75
    Ind = np.where((outputM <= 0.7501) & (outputM >= 0.7499))
    Ind = int(np.mean(Ind))
    
    reactionLen075 = outputDistance[Ind]
    #############################################################
    
    if (k == 1):
        maxx = outputInd_len_ZND*5;
    else:
        maxx = outputDistance[k];

    minx = 0;

    maxT = np.max(outputT)+0.1*np.min(outputT);
    minT = np.min(outputT)-0.1*np.min(outputT);
    maxP = np.max(outputP)+0.1*np.min(outputP);
    minP = np.min(outputP)-0.1*np.min(outputP);

    # %Check for Eigenvalue Detonation

    if(n == b):
        print 'Error: Maximum thermicity occurs at the end of the reaction zone'
        print 'You may have an eigenvalue detonation, your final integration length may be too short,'
        print' your mixture may be too rich/lean, or something else may be wrong'
        print ' '
        print 'Mach Number (end of reaction): ' +str(outputM[b])+ ' - if close to 1, check for eigenvalue detonation'
        outputIind_time_ZND = outputTime[b]; 
        outputInd_len_ZND = outputDistance[b]; 
        outputExo_time_ZND = 0; 
        outputExo_len_ZND = 0;  
        print 'Induction Time: ' +str(outputInd_time_ZND)+'';
        print 'Exothermic Pulse Time: ' +str(outputExo_time_ZND)+'';

    elif (n == 1):
        print 'Error: Maximum thermicity occurs at the beginning of the reaction zone'
        print 'You may have an eigenvalue detonation, your final integration length may be too short,'
        print 'your mixture may be too rich/lean, or something else may be wrong'
        print ' '
        print 'Mach Number (end of reaction): ' +str(outputM[b])+ ' - if close to 1, check for eigenvalue detonation'
                                                     
        outputInd_time_ZND = outputTime[0]; 
        outputInd_len_ZND = outputDistance[0]; 
        outputExo_time_ZND = 0; 
        outputExo_len_ZND = 0;  
        print 'Induction Time: '+str(outputInd_time_ZND)+'';
        print 'Exothermic Pulse Time: '+str(outputExo_time_ZND)+'';
    else:
        print 'There is no Error'
        max_sigmadot = max(outputThermicity); # max thermicity
        half_sigmadot_flag1 = 0;
        half_sigmadot_flag2 = 0;
        # Go into a loop to find two times when sigma_dot is half its maximum
        for j in range(0,b):
            if(half_sigmadot_flag1 == 0):
                
                if(outputThermicity[j] > 0.5* max_sigmadot):
                    
                    half_sigmadot_flag1 = 1;
                    tstep1 = j;
                
            elif(half_sigmadot_flag2 == 0):
                if(outputThermicity[j] < 0.5* max_sigmadot):
                    half_sigmadot_flag2 = 1;
                    tstep2 = j;
                else:
                    tstep2 = 0;

    #Exothermic time for ZND detonation
    if(tstep2 == 0):
        print 'Error: No pulse in the thermicity'
        print '       You may have an eigenvalue detonation, your final integration length may be too short,'
        print '       your mixture may be too rich/lean, or something else may be wrong'
        outputExo_time_ZND = 0; 
        outputExo_len_ZND = 0;  
    else:
        outputExo_time_ZND = outputTime[tstep2] - outputTime[tstep1]; 
        outputExo_len_ZND = outputDistance[tstep2] - outputDistance[tstep1];
    
    return np.array([reactionLen09,reactionLen075,rZoneGradTemp]) # in m

def soundSpeedM(a):
    # Matlab  SOUNDSPEED  - Speed of sound (m/s).
    gamma = a.cp_mass/a.cv_mass;
    wtm = a.mean_molecular_weight;
    r = sd.gas_constant/wtm;
    b = np.sqrt(gamma * r * a.T);
    return b

##########################################################################

###########
# EXAMPLE #
###########
P1 = sd.one_atm;
T1 = 300;
q = 'H2:2 O2:1 N2:3.76'
mech = 'h2air_highT.cti'

out = znd_CJ(P1, T1, q, mech)
print 'Lengths in m:'
print out[1]
