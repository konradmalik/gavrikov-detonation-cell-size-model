from SDToolbox import *
import math
import numpy as np

def explosion(gas,fname,out,b,j):
    """

    explosion.m
    Computes the time evolution of a constant volume explosion
    
    FUNCTION
    SYNTAX
    [ind_time,ind_time_10,ind_time_90,time,temp,press,species] = explosion(gas,fig_num)
    
    INPUT
    gas = working gas object
    fig_num = figure number (0 for no plot)
    
    OUTPUT
    exo_time = pulse width (in secs) of temperature gradient (using 1/2 max)
    ind_time = time to maximum temperature gradient
    ind_len = distance to maximum temperature gradient
    ind_time_10 = time to 10% of maximum temperature gradient
    ind_time_90 = time to 90% of maximum temperature gradient
    time = array of time
    temp = temperature profile array
    press = pressure profile array
    species = matrix of species profiles

    """
    #b = 200001; j = gas.nSpecies();
    rho = gas.density
    
    r = Reactor(gas)
    sim = ReactorNet([r])
    t = 0.0
    
    # optional settigns:
    sim.set_initial_time(1.0e-8)
    sim.set_max_time_step(1.0e-5)
    
    temp_grad = zeros(b,float)
    temp_grad_max=0.0
    temp_max = 0.0
    y = zeros(j,float)

    #EXTRACT OUTPUT INFORMATION
    for n in range(b):
        out.time[n] = t
        out.T[n] = r.T
        out.P[n] = r.thermo.P
        #print 't = ', t, ' T = ', r.temperature(), ' P= ', r.pressure()
        for i in range(j):
            out.species[n,i] = r.Y[i]
            y[i] = r.Y[i]
        
        gas.TDY = out.T[n], r.density, y
        P = gas.P/one_atm

        # FIND TEMPERATURE GRADIENT
        # Conservation of Energy implies that e0 = e(T)
        # e = cv*T; dedt = 0; cv*dTdt + sum(ei(T)*dyidt) = 0; dyidt = wdoti*wti/rho
        # dTdt = -sum(ei(T)*wdoti*wti)/(cv*rho)
        cv = gas.cv_mass;
        wdot = gas.net_production_rates;
        mw = gas.molecular_weights;
        hs = gas.standard_enthalpies_RT
        R = gas_constant;
        wt = gas.mean_molecular_weight
        
        sumT = 0.0
        for z in range(j):
            w = mw[z]; e = R*out.T[n]*(hs[z]/w - 1/wt)
            wd = wdot[z]; sumT = sumT + e*wd*w;
        temp_grad[n] = -sumT/(rho*cv)
        sim.step()
        t = sim.time
        
        #RUDY - braking simulation if decreasing max grad is achieved
        # calc
        
        if temp_grad[n] > temp_grad_max or temp_grad_max == 0.0:
            temp_grad_max = temp_grad[n]
        if r.T > temp_max:
            temp_max = r.T
        if temp_max > out.T[0]+400 and temp_grad[n] < 0.1* temp_grad_max:
            #print 'Loop stops at', t, 's'
            break
           #return out
        ## end of Rudy's mods.

    del sim
    del r
    
    #FIND INDUCTION TIME - MAXIMUM TEMPERATURE GRADIENT
    k = 0; MAX = max(temp_grad); d = temp_grad[0]; HMPWt = zeros(2,float)
    if d == MAX:
        print 'Initial Gradient is Maximum - post shock temperature may be too low'
        return gas
    while d < MAX:
        k = k + 1; d = temp_grad[k];
    out.ind_time = out.time[k]; k1 = k; k = 0;
    MAX10 = 0.1*MAX; d = temp_grad[0];
    while(d < MAX10 and k < b-1):
        k = k + 1; d = temp_grad[k];
    if(k == b):
        print 'MAX10 may be incorrect - reached end of array'
    out.ind_time_10 = out.time[k]; k = 0;
    MAX90 = 0.9*MAX; d = temp_grad[0];
    while(d < MAX90 and k < b-1):
        k = k + 1; d = temp_grad[k];
    if(k == b-1):
        print 'MAX90 may be incorrect - reached end of array'
    out.ind_time_90 = out.time[k];
    #print 'MAX_temp_grad at time=', out.ind_time

    #MAX_TEMP_GRAD_TIME = out.ind_time
    # find exothermic time
    half_T_flag1 = 0;
    half_T_flag2 = 0;
    #Go into a loop to find two times when Temperature is half its maximum
    for j in range(b) : 
        if (half_T_flag1 == 0 ):
            if (temp_grad[j] >= (0.5* MAX)):
                half_T_flag1 = 1;
                tstep1 = j;
        else:
            if (half_T_flag2 == 0) :
                if (temp_grad[j] <= (0.5* MAX) ):
                    half_T_flag2 = 1;
                    tstep2 = j;
                else:
                    tstep2 = 0;
    ##Exothermic time for constant volume explosion
    ##old method:
    out.exo_time = out.time[tstep2] - out.time[tstep1];
    ##new method by WRudy (linear approximation)
    a1 = (temp_grad[tstep1]-temp_grad[tstep1-1])/(out.time[tstep1]-out.time[tstep1-1])
    b1 = temp_grad[tstep1]-a1*out.time[tstep1]
    # calculate x1 (time), the hipothetical intersection point between 0.5*max(Thermicity) and line with equation y=a1*x+b1
    x1 = (0.5*MAX-(temp_grad[tstep1]-a1*out.time[tstep1]))/a1
    # similar calculations as above but for decreasing thermicity function
    a2 = (temp_grad[tstep2]-temp_grad[tstep2-1])/(out.time[tstep2]-out.time[tstep2-1])
    b2 = temp_grad[tstep2]-a2*out.time[tstep2]
    x2 = (0.5*MAX-(temp_grad[tstep2]-a2*out.time[tstep2]))/a2
    #and finally:
    out.exo_time2 = x2-x1

    #print 'outExo_time =', out.exo_time
    #OUTEXO_TIME = out.exo_time
    
    if fname==0:
        return out

    else:
        k = 0; MAX = max(out.T); d = out.T[0];
        while d < MAX:
            k = k + 1; d = out.T[k];
        if out.time[k] == 0:
            maxt = out.ind_time*5;
        elif out.time[k] >= out.ind_time*50:
            maxt = out.ind_time*5;
        else:
            maxt = out.time[k] + 0.1*out.time[k];
        mint = 0;
        maxT = max(out.T)+0.1*min(out.T); minT = min(out.T)-0.1*min(out.T); 
        maxP = max(out.P)+0.1*min(out.P); minP = min(out.P)-0.1*min(out.P); 
        maxpw = HMPWt[1] + 0.1*HMPWt[1]; minpw = HMPWt[0] - 0.1*HMPWt[0]; 
        maxTG = max(temp_grad) + 0.1*abs(max(temp_grad));
        minTG = min(temp_grad)-0.1*abs(min(temp_grad));
        d = datetime.date.today(); P = out.P[0]/OneAtm;

        fn = fname + '_CVprofile.plt';
        outfile = file(fn, 'w');
        outfile.write('# CONSTANT VOLUME PROFILES\n');
        outfile.write('# CALCULATION RUN ON %s\n\n' % d);
        outfile.write('# Maximum time calculated = %.4e\n' % max(out.time))
        outfile.write('# t_min = %.4f, t_max = %.4e\n' % (mint, maxt))
        outfile.write('# T_min = %.2f, T_max = %.2f\n' % (minT, maxT))
        outfile.write('# P_min = %.2f, P_max = %.2f\n' % (minP, maxP))
        outfile.write('# TG_min = %.2f, TG_max = %.2f\n' % (minTG, maxTG))
        outfile.write('# TG_min = %.2f, TG_max = %.2f\n' % (minTG, maxTG))
        outfile.write('# maxTgrad_time = %.2f\n'  % (out.ind_time))
        outfile.write('# maxExo_time = %.2f\n'  % (out.exo_time))
        outfile.write('# THE OUTPUT DATA COLUMNS ARE:\n');
        outfile.write('Variables = "Time", "Temperature", "Pressure", "temp_grad"\n');
        for i in range(b):
            outfile.write('%.4E \t %.4E \t %.4E \t %.4E\n'% (out.time[i],out.T[i],out.P[i],temp_grad[i]));

        return out

def ea_r(Pinit, Tinit, q, mech):
    gas1 = Solution(mech);
    gas2 = Solution(mech);
    gas3 = Solution(mech);
    
    # chapman - jouguet detonation
    [cj_speed,_] = CJspeed(Pinit, Tinit, q, mech, 0);  
    
    gas1 = PostShock_fr(cj_speed, Pinit, Tinit, q, mech)
    T1 = gas1.T#its Von Neumann temp
    
    # 1.6 D
    gas2 = PostShock_fr(1.6*cj_speed, Pinit, Tinit, q, mech)
    T2 = gas2.T
    
    # 1.3 = (1.6+1.0)/2.0 D
    gas3 = PostShock_fr(1.3*cj_speed, Pinit, Tinit, q, mech)
    T3 = gas3.T
    
    # SOLVE CONSTANT VOLUME EXPLOSION ODES
    # b must be large
    b = 1000000; j = gas1.n_species;
    out1 = cvoutput(b,j)
    out1 = explosion(gas1,0,out1,b,j);
    ind_time1 = out1.ind_time
    exo_time1 = out1.exo_time2

    j = gas2.n_species;
    out2 = cvoutput(b,j)
    out2 = explosion(gas2,0,out2,b,j);
    ind_time2 = out2.ind_time
    exo_time2 = out2.exo_time2
    
    #Effective Activation Energy * Tvn
    Ea_R = log(ind_time1/ind_time2) / ((1.0/T1) - (1.0/T2))
    #Post shock temp
    Tps = T3
    #Von neumann temp
    Tvn = T1
    
    return [Ea_R,Tps,Tvn]
