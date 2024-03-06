#!/usr/bin/env python
# coding: utf-8

# In[12]:


import numpy as np
import copy 
G = 6.67480e-18 # cm3/gs
pi = 3.14159 
c= 2.99792e10 # cm/s
AU = 1.49598e13 # cm
M = 1.98854e33 # g

# Let G as 1 then find constant for translating degree to km/s , s
def to_phy_vel (v_sim):
    return v_sim/100000.0 *np.sqrt(G*M/AU)
def to_sim_vel(v_phy):
    return v_phy*100000.0 / np.sqrt(G*M/AU)
def to_phy_time (t_sim):
    return t_sim/ np.sqrt(G*M/AU**3)
def to_sim_time (t_phy):
    return t_phy * np.sqrt(G*M/AU**3)
#ts =  G * M *(c**-3)
#print(tv, "\n",ts)

def init(m1 , m2 , a , e ):  
    tot_m = m1 + m2
    mu = (m1 *m2) / tot_m
    k = m1 *m2
    r = a*(1+e)
    l = np.sqrt(mu*k*a*(1. -e * e))
    theta_dot = l/mu/r/r
    
    x1 = -r * m2 / tot_m
    x2 = r + x1
    # x2 = r*m1 /tot_m
    v1 =  theta_dot * x1
    v2 = theta_dot * x2
    
    period = np.sqrt(4. * pi**2 *mu /k * a *a *a)
    
    return x1,x2,v1,v2,period
def new_body (mass,pos,vel):
    return [mass,pos,vel]  # store mass to index 0, pos to ind 1, vel to ind2
def vec_sub (v1,v2):
    return np.array([v1[0]-v2[0], v1[1]-v2[1]],dtype = np.float64) 
def dot_p (v1,v2): 
    term = v1[0]*v2[0] + v1[1]*v2[1]           
    return term     
def size (v): return np.sqrt(v[0]**2 + v[1]**2) 
def acce (body, system):
    res = np.array([0.,0.])
    for i in range(len(system)):
            if i != body:    
                pos_rel = vec_sub(system[int(body)][1],system[i][1])
                r = size(pos_rel)
                res[0] -= system[i][0] * pos_rel[0] /r**3
                res[1] -= system[i][0] * pos_rel[1] /r**3
     
    return res        

def acce_pn (body, system):
    res = np.array([0.,0.])
    for i in range(len(system)):
        if i != body:
            pos_rel = vec_sub(system[body][1],system[i][1])
            vel_rel = vec_sub(system[body][2],system[i][2])
            r = size(pos_rel)
            v = size(vel_rel)
            tot_m = system[body][0] + system[i][0]
            eta = system[body][0] * system[i][0] / tot_m**2
            r_dot = dot_p(pos_rel,vel_rel) /r
            pn_a = 8./5. * eta * tot_m / r * r_dot * (17./3. * tot_m / r + 3. * v**2) 
            pn_b = - 8./5. * eta * tot_m / r * (3. * tot_m/ r + v * v )
            
            res += system[i][0] /r**2 * (pn_a * pos_rel /r + pn_b * vel_rel )
    res = res / to_sim_vel(c)**5        
    return res      

def acce_dot (body, system):
    res = np.array([0.,0.])
    for i in range(len(system)):
        if i != body:
            pos_rel = vec_sub(system[body][1],system[i][1])
            r = size(pos_rel)
            vel_rel = vec_sub(system[body][2],system[i][2])
            r_dot = dot_p(pos_rel, vel_rel) /r 
            res -= system[i][0] * (vel_rel / r**3 - 3. * r_dot* pos_rel/r**4 )
    return res          

def acce_pn_dot (body, system):
    res = np.array([0.,0.])
    for i in range(len(system)):
        if i != body:
            pos_rel = vec_sub(system[body][1],system[i][1])
            vel_rel = vec_sub(system[body][2],system[i][2])
            acce_rel = acce(body, system)
            r = size(pos_rel)
            v = size(vel_rel)
            tot_m = system[body][0] + system[i][0]
            eta = system[body][0] * system[i][0] / tot_m**2
            r_dot = dot_p(pos_rel,vel_rel) /r 
            v_dot = dot_p(vel_rel, acce_rel) / v
            pn_a = 8./5. * eta * tot_m / r * r_dot * (17./3. * tot_m/r + 3.*(v**2) )
            pn_b = -8./5. * eta * tot_m /r * (3*tot_m / r + v**2)
            r_ddot = (v**2 + dot_p(pos_rel, acce_rel)- r_dot**2) / r
            pn_a_dot = 8. * tot_m * eta / 5. /r * (17. * tot_m / 3. / r + 3. * v**2)* r_ddot 
            - 8.*tot_m*eta /5. /r/r *(17.* tot_m/ 3./r + 3.* v**2)* r_dot**2
            + 8.*tot_m*eta / 5./ r * (-17. *tot_m* r_dot / 3./ r /r + 6. * v* v_dot)*r_dot
            pn_b_dot = 8./5. * eta * tot_m /r /r * r_dot * (3. * tot_m/ r + v**2 )
            -8./5. * eta *tot_m /r * (-3.* tot_m /r/ r * r_dot + 2. *dot_p(vel_rel, acce_rel))
            
            res += -2. * system[i][0] /r**3 * (pn_a* pos_rel /r + pn_b *vel_rel)*r_dot
            + system[i][0]/ r**2 * (-pn_a*pos_rel/r**2 * r_dot + pn_a *vel_rel/r + pn_b*acce_rel + pos_rel*pn_a_dot/r + vel_rel*pn_b_dot)
    
    res = res / to_sim_vel(c)**5        
    return res    

def integrater( system, dt, pn):
        a0 = []
        a0_dot = []
        for i in range(len(system)):
                a0.append( acce(i ,system ) )
                a0_dot.append( acce_dot(i,system ) )
        a0 = np.array(a0)
        a0_dot = np.array(a0_dot)        
            
        if pn is True:
            for i in range(len(system)):
                a0[i] += acce_pn(i,system )
                a0_dot[i] += acce_pn_dot(i, system )    
        #print(a0,'\n', a0_dot)
        system_p = copy.deepcopy(system)
        
        for i in range (len(system)):
            system_p[i][1] += system[i][2] * dt + a0[i] * dt**2 /2. + a0_dot[i] * dt**3 / 6.
            system_p[i][2] += a0[i] * dt + a0_dot[i] * dt**2 /2.
        
        ap = []
        ap_dot =  []
        for i in range(len(system)):
                ap.append( acce(i ,system_p ) )
                ap_dot.append( acce_dot(i ,system_p ) )
        ap= np.array(ap)
        ap_dot = np.array(ap_dot)         
        if pn is True:
            for i in range(len(system)):
                ap[i] += acce_pn(i, system_p )
                ap_dot[i] += acce_pn_dot(i, system_p )
                    
        
        a0_ddot = [] 
        a0_dddot = [] 
        
        for i in range(len(system)):
            a0_ddot.append( -6.* (a0[i]- ap[i])/dt**2 -2.*(2*a0_dot[i]+ ap_dot[i])/dt) 
            
            a0_dddot.append( 12.*(a0[i]-ap[i] )/dt**3 + 6.*(a0_dot[i] + ap_dot[i])/dt**2)
            
        a0_ddot = np.array(a0_ddot)
        a0_dddot = np.array(a0_dddot)    
        for i in range(len(system)):
            system_p[i][1] += a0_ddot[i]* dt**4 / 24. + a0_dddot[i] * dt**5 / 120.
            system_p[i][2] += a0_ddot[i]* dt**3  / 6. + a0_dddot[i]* dt**4 / 24.
            
        return system_p
     
def solver(m1,m2,a,e, n_period, n_step, pn):
    x1,x2,v1,v2,period  = init(m1,m2,a,e )
    dt = 0.1
    tot_t = period * n_period
    acc_t = 0.
    
    x1_his =  []
    x2_his = [] 
    y1_his = [] 
    y2_his = [] 
    system = [new_body(m1, np.array([x1,0.]),np.array([0.,v1])),
              new_body(m2, np.array([x2,0.]),np.array([0.,v2]) )]           
    
    sch_r = 2.*(m1+m2) / to_sim_vel(c)**2
    print('initial value:')
    print('body1 mass:', system[0][0], '\n', 'initial_pos: ', system[0][1],'init vel:', system[0][2]) 
    print('body2 mass:', system[1][0], '\n', 'initial_pos: ', system[1][1],'init vel:', system[1][2])
    print('seperation:', system[0][1]- system[1][1])
    print('target time: ', tot_t)
    print('R_sch of total mass: ', sch_r)
    
    dt_alpha = 30. * (1. + e)*(1. + int(pn))
    if pn is False:
        pos = size(system[0][1]- system[1][1]) 
        vel = size(system[0][2] - system[1][2])
        est_dt1 = pos / vel / dt_alpha
        est_dt2 = (2.*a  - pos )**2 / vel / pos /dt_alpha
        est_n1 = tot_t / (est_dt1 + est_dt2) * 2.
        est_n2 = tot_t / (est_dt1 * est_dt2 )**0.5
        est_n  = (est_n1 + est_n2) * 0.5
        print("Est N :", est_n)
        if abs(est_n-float(n_step))/float(n_step) * 100. > 30.:
            dt_alpha *= float(n_step)/ est_n / 1.2 
    
    
    print("Start Solve")
    early_break = False
    
    for it in range(1, n_step+1):
        pos = size(system[0][1]- system[1][1]) 
        if pos < sch_r:
            early_break = True
            break
        
        vel = size(system[0][2]- system[1][2])   
        dt = pos / vel / dt_alpha
        if acc_t + dt > tot_t:
            dt = acc_t - tot_t
            
        new_system = integrater(system,dt, pn)
        system , new_system = new_system, system
        
        if it % 1 ==0:
            x1_his.append( system[0][1][0])
            y1_his.append( system[0][1][1])
            x2_his.append( system[1][1][0])
            y2_his.append( system[1][1][1])
             
        
        acc_t +=dt 
        if acc_t >= tot_t:
            print('Total t achived at iter', it)
            early_break = True
            break
        
    if early_break is False:
        print('complete all iter', n_step)
        
    pos = size(system[0][1]-system[1][1] )
    print('Final value')
    print('body1 mass:', system[0][0], '\n', 'fianl_pos: ', system[0][1],'final vel:', system[0][2]) 
    print('body2 mass:', system[1][0], '\n', 'final_pos: ', system[1][1],'final vel:', system[1][2])
    print('seperation:', pos)
    print('total time: ', acc_t, to_phy_time(acc_t))
    print('body1 vel', to_phy_vel(size(system[0][2]))) 
    print('body2 vel', to_phy_vel(size(system[1][2])))
    
    f = open(  "result.txt", "w+")
    
    #for i in range(len(system)):
        #f.write("body%d mass: %f \n body%d pos: %f \n body%d vel: %f" %i %system[i][0] %i %system[i][1] %i %system[i][2] )
    
    f.write("body1 mass: %f\n body1 pos: %s\n body1 vel: %s\n"  %(system[0][0],system[0][1],system[0][1]))
    f.write("body2 mass: %s\n body2 pos: %s\n body2 vel: %s"  %(system[1][0],system[1][1],system[1][1]))
    
    f.close    
     
    x1_his= np.array(x1_his)
    y1_his = np.array(y1_his)
    x2_his = np.array(x2_his)
    y2_his = np.array(y2_his)
    
    return[ x1_his,y1_his, x2_his, y2_his]

    
    
k = 2.22*(10**17) # if mass is given with g translate to solar mass M
m_h = k/M
sun_x_his, sun_y_his, h_x_his, h_y_his = solver(1.,m_h, 17.8,0.967,50,40000, False)


# In[ ]:




