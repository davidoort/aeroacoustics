from math import *
import numpy as np
from scipy.integrate import trapz
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

maxcount = 40

verbose = True

## 3D velocity induced by a vortex filament
def velocity_3D_from_vortex_filament(gamma,X1,X2,XP,core):
  """
    params:
      gamma - vorticity strength along filament
      X1,X2 - 3D vectors of ends of filament
      XP - 3D vector of control point
  """
  # function to calculate the velocity induced by a straight 3D vortex filament
  # with circulation gamma at a point XP. The geometry of the vortex filament
  # is defined by its edges: the filament starts at X1 and ends at X2.
  # the input core defines a vortex core radius, inside which the velocity
  # is defined as a solid body rotation.
  # The function is adapted from the algorithm presented in:
  #                Katz, Joseph, and Allen Plotkin. Low-speed aerodynamics.
  #                Vol. 13. Cambridge university press, 2001.


  # calculate geometric relations for integral of the velocity induced by filament
  R1 = sqrt((XP[0]-X1[0])**2 + (XP[1]-X1[1])**2 + (XP[2]-X1[2])**2)
  R2 = sqrt((XP[0]-X2[0])**2 + (XP[1]-X2[1])**2 + (XP[2]-X2[2])**2)
  R1XR2_X = (XP[1]-X1[1])*(XP[2]-X2[2]) - (XP[2]-X1[2])*(XP[1]-X2[1])
  R1XR2_Y = -(XP[0]-X1[0])*(XP[2]-X2[2]) + (XP[2]-X1[2])*(XP[0]-X2[0])
  R1XR2_Z = (XP[0]-X1[0])*(XP[1]-X2[1]) - (XP[1]-X1[1])*(XP[0]-X2[0])
  R1XR_SQR = R1XR2_X**2 + R1XR2_Y**2 + R1XR2_Z**2
  R0R1 = (X2[0]-X1[0])*(XP[0]-X1[0]) + (X2[1]-X1[1])*(XP[1]-X1[1]) + (X2[2]-X1[2])*(XP[2]-X1[2])
  R0R2 = (X2[0]-X1[0])*(XP[0]-X2[0]) + (X2[1]-X1[1])*(XP[1]-X2[1]) + (X2[2]-X1[2])*(XP[2]-X2[2])
  #// check if target point is in the vortex filament core,
  #// and modify to solid body rotation
  if (R1XR_SQR < core**2):
    R1XR_SQR = core**2
    #// gamma = 0
  if (R1 < core):
    R1 = core
    #// gamma = 0
  if (R2 < core):
    R2 = core
    #// gamma = 0
  #// determine scalar
  K = gamma/(4.0*pi*R1XR_SQR)*(R0R1/R1 - R0R2/R2)

  # output results, vector with the three velocity components
  return [K*R1XR2_X, K*R1XR2_Y, K*R1XR2_Z]
# end velocity_3D_from_vortex_filament

def solve_lifting_line_system_matrix_approach(wake_system,wind,RPM,rotorradius,rho=1.225):
    # this codes solves a lifting line model of a horizontal axis rotor
    """
      params:
        rotor_wake_system: data structure that contains the geometry of the horseshoe vortex rings,
                            and the control points at the blade
        wind: unperturbed windstream velocity, also known as U_infinity
        Omega: rotational velocity of the rotor rad/s
        rotorradius: the radius of the rotor
    """
    ## get controlpoints data structure
    controlpoints = wake_system.controlPoints
    ## get horseshoe vortex rings data structure
    rings = wake_system.rings
    Omega = RPM*pi/30.0
    fwdVel = np.linalg.norm(wind)
    angleOfAttack = np.arctan2(np.linalg.norm(wind[1:]),wind[0])
    ## initialize variables that we will use during the calculation
    velocity_induced =[] # velocity induced by a horseshoe vortex ring at a control point
    #up = []
    #vp = []
    #wp = [] # components of the velocity induced by one horseshoe vortex ring
    radialposition, azimdir = 0, 0 # radial position of the control point
    alpha = 0 # angle of attack
    GammaNew = [0]*len(controlpoints) # new estimate of bound circulation, initialized as zeroes
    Gamma = [0]*len(controlpoints) # current solution of bound circulation
    vel1, vmag, vaxial, vazim, temploads = 0,0,0,0,0 # velocities and loads at controlpoint
    MatrixU = np.zeros((len(controlpoints),len(rings))) # matrix of induction, for velocity component in x-direction
    MatrixV = np.zeros((len(controlpoints),len(rings))) # matrix of induction, for velocity component in y-direction
    MatrixW = np.zeros((len(controlpoints),len(rings))) # matrix of induction, for velocity component in z-direction

    # the variables below are to setup the maximum number of iterations and convergence criteria
    Niterations = 500
    errorlimit = 0.01
    error = 1.0
    refererror = 1.0
    ConvWeight = 0.2
    thrust = 0
    prevThrust = 0

    # initalize and calculate matrices for velocity induced by horseshoe vortex rings
    # nested for loops, each row varying with controlpoint "icp", each column varying with
    # horseshoe vortex ring "jring"
    for icp in range(0,len(controlpoints)):
      #MatrixU[icp] = new Array(); // new line of matrix
      #MatrixV[icp] = new Array(); // new line of matrix
      #MatrixW[icp] = new Array(); // new line of matrix
      for jring in range(0,len(rings)):
        # set ring strenth to unity, to calculate velocity induced by horseshoe vortex ring "jring"
        # at controlpoint "icp"
        rings[jring] = update_Gamma_single_ring(rings[jring],1.0,1.0)
        velocity_induced = velocity_induced_single_ring(rings[jring],controlpoints[icp])
        # add compnent of velocity per unit strength of circulation to induction matrix
        MatrixU[icp][jring] = velocity_induced[0]
        MatrixV[icp][jring] = velocity_induced[1]
        MatrixW[icp][jring] = velocity_induced[2]
      # end for
    # end for

    # calculate solution through an iterative process
    for kiter in range(1,Niterations+1):
      # copy the array, since arrays are mutable
      #for ig in range(0,len(GammaNew)):
      #  Gamma[ig] = GammaNew[ig] # update current bound circulation with new estimate
      # end for
      Gamma = np.copy(GammaNew)

      ## also need to update filament locations based on new induced downwash...
      Nblades = wake_system.rotor.blades
      ## re-initialize output variables
      a_temp = [[] for b in range(Nblades)] # output vector for axial induction
      aline_temp = [[] for b in range(Nblades)]  # output vector for azimuthal induction
      r_temp = [[] for b in range(Nblades)]  # output vector for radial position
      Fnorm_temp = [[] for b in range(Nblades)]  # output vector for axial force
      Ftan_temp = [[] for b in range(Nblades)]  # output vector for tangential force
      Gamma_temp = [[] for b in range(Nblades)]  # output vector for circulation
      local_temp = [[] for b in range(Nblades)]  # output vector for circulation
      torque_temp = [[] for b in range(Nblades)]
      Cl_temp = [[] for b in range(Nblades)]
      Cd_temp = [[] for b in range(Nblades)]

      # calculate velocity, circulation and loads at the controlpoints
      for icp in range(0,len(controlpoints)):
        # determine radial position of the controlpoint;
        # radialposition = np.sqrt(dot(controlpoints[icp].coordinates, controlpoints[icp].coordinates))
        radialposition = controlpoints[icp].radius
        u, v, w = 0, 0, 0 # initialize velocity
        # multiply icp line of Matrix with vector of circulation Gamma to calculate velocity at controlpoint
        for jring in range(0,len(rings)):
          u += MatrixU[icp][jring]*Gamma[jring] # axial component of velocity
          v += MatrixV[icp][jring]*Gamma[jring] # y-component of velocity
          w += MatrixW[icp][jring]*Gamma[jring] # z-component of velocity
        # end for rings
        # calculate total perceived velocity FIX MEEEEE doesn't work for reverse rotation!!!!
        vrot = cross3([-Omega*wake_system.direction, 0, 0], controlpoints[icp].coordinates ) # rotational velocity
        vel1 = [wind[0] + u + vrot[0], wind[1] + v + vrot[1], wind[2] + w + vrot[2]] # total perceived velocity at section
        # calculate azimuthal and axial velocity
        azimdir = cross3([-1.0*wake_system.direction/radialposition, 0, 0], controlpoints[icp].coordinates ) # rotational direction
        vazim = dot(azimdir, vel1) # azimuthal direction
        vaxial =  dot([1.0, 0, 0], vel1) # axial velocity
        # calculate loads using blade element theory
        temploads = loadBladeElement(vaxial,vazim,radialposition,wake_system.rotor,wake_system.pitch,rho)
        # new point of new estimate of circulation for the blade section
        GammaNew[icp] = temploads[2]*wake_system.direction
        # update output vector
        b = controlpoints[icp].bladeID
        a_temp[b].append(u + vrot[0]) #-(u + vrot[0])/wind[0]
        aline_temp[b].append(vazim/(radialposition*Omega)-1)
        r_temp[b].append(radialposition)
        Fnorm_temp[b].append(temploads[0])
        Ftan_temp[b].append(temploads[1])
        Gamma_temp[b].append(temploads[2]) ## ?why? already have GammaNew
        local_temp[b].append(temploads[3])
        torque_temp[b].append(temploads[1]*radialposition)
        Cl_temp[b].append(temploads[4])
        Cd_temp[b].append(temploads[5])
      # end loop control points

      # check convergence of solution
      refererror = np.max(np.absolute(GammaNew)) #gets the max Gamma value
      refererror = max(refererror,0.001) # define scale of bound circulation, 0.001 to make sure it isn't divide by 0
      diff = np.average(np.absolute(np.subtract(GammaNew,Gamma))) # difference betweeen iterations
      error = diff/refererror # relative error

      # use thrust to calculate new average induced velocity, and update wake system
      prevThrust = thrust
      thrust = 0
      for i in range(len(r_temp)):
        thrust += trapz(Fnorm_temp[i],r_temp[i])

      tError = np.absolute(thrust-prevThrust)/thrust
      ## for checking convergence
      #if kiter < 20:
        #print(kiter,error,tError,(thrust-prevThrust)/thrust)
        #plt.plot(r_temp[0],Cd_temp[0])

      # check if GAMMA is sufficiently converged
      if error < errorlimit:
        # check to see if inflow is close enough
        # get the induced velocity/inflow ratio
        Vi = wake_system.rotor.Vi_hover(thrust,rho)
        inflowRatio = wake_system.rotor.inflowRatio(fwdVel,angleOfAttack,Vi,RPM)
        if np.absolute(inflowRatio-wake_system.inflowRatio)/inflowRatio < 0.01 and (tError < errorlimit):
          # if error smaller than limit, stop iteration cycle
          # kiter < Niterations
          print(kiter,inflowRatio,wake_system.inflowRatio,(inflowRatio-wake_system.inflowRatio)/inflowRatio)
          break
        else:
          # update the wake model
          print("update wake ",kiter,inflowRatio,wake_system.inflowRatio,error,tError,thrust,prevThrust)
          # requires some weight for good convergence!!!
          wake_system.update_Inflow(0.5*inflowRatio+0.5*wake_system.inflowRatio)

      # this seems to work better for all cases... correct to catch overshoots
      maxCW = 0.05
      minCW = 0.03
      ConvWeight = min(maxCW,error)
      ConvWeight = max(ConvWeight,minCW)
      # set new estimate of bound circulation
      for ig in range(0,len(GammaNew)):
        GammaNew[ig] = (1 - ConvWeight)*Gamma[ig] + ConvWeight*GammaNew[ig]
      # end for
    # end iteration loop
    if kiter == Niterations: 
      print("max iterations!",kiter,error,tError,thrust,wake_system.inflowRatio)
    
    # get some total values, lift, torque, power...
    totalLift = 0
    totalTorque = 0
    totalPower = 0
    for i in range(len(r_temp)):
      totalLift += trapz(Fnorm_temp[i],r_temp[i])
      totalTorque += trapz(torque_temp[i],r_temp[i])
      totalPower += trapz(torque_temp[i],r_temp[i])*Omega

    #plt.show()
    #// output results of converged solution
    return({"Vi":a_temp,"Vi_azim":aline_temp,"r":r_temp,"Fnorm":Fnorm_temp,"Ftan":Ftan_temp,"Gamma":Gamma_temp,"local":local_temp,"Power":totalPower,"Lift":totalLift,"Torque":totalTorque});


class LLCoaxOutput:

  def __init__(self,Nblades):
    self.ax = [[] for b in range(Nblades)] # output vector for axial induction
    self.az = [[] for b in range(Nblades)]  # output vector for azimuthal induction
    self.r = [[] for b in range(Nblades)]  # output vector for radial position
    self.Fnorm = [[] for b in range(Nblades)]  # output vector for axial force
    self.Ftan = [[] for b in range(Nblades)]  # output vector for tangential force
    self.Gamma = [[] for b in range(Nblades)]  # output vector for circulation
    self.local = [[] for b in range(Nblades)]  # output vector for local angle of attack
    self.torque = [[] for b in range(Nblades)]  # output vector for torque (per element)
    self.Cl = [[] for b in range(Nblades)]  # output vector for Cl per element
    self.Cd = [[] for b in range(Nblades)]  # output vector for Cd per element

def init_Vi_Matrices(controlpoints,rings,MatrixU,MatrixV,MatrixW):
  ## Matrices mst be mutable!!!

  # initalize and calculate matrices for velocity induced by horseshoe vortex rings
  # nested for loops, each row varying with controlpoint "icp", each column varying with
  # horseshoe vortex ring "jring"
  for icp in range(0,len(controlpoints)):
    #MatrixU[icp] = new Array(); // new line of matrix
    #MatrixV[icp] = new Array(); // new line of matrix
    #MatrixW[icp] = new Array(); // new line of matrix
    for jring in range(0,len(rings)):
      # set ring strenth to unity, to calculate velocity induced by horseshoe vortex ring "jring"
      # at controlpoint "icp"
      rings[jring] = update_Gamma_single_ring(rings[jring],1.0,1.0)
      velocity_induced = velocity_induced_single_ring(rings[jring],controlpoints[icp])
      # add compnent of velocity per unit strength of circulation to induction matrix
      MatrixU[icp][jring] = velocity_induced[0]
      MatrixV[icp][jring] = velocity_induced[1]
      MatrixW[icp][jring] = velocity_induced[2]
    # end for
  # end for
  return

def calc_ControlPoints(wakeSystem,controlpoints,rings1,Gamma1,rings2,Gamma2,MatrixU,MatrixV,MatrixW,Omega,wind,LLoutput,rho=1.225):
  # calculate velocity, circulation and loads at the controlpoints
  GammaNew = [0]*len(controlpoints)
  for icp in range(0,len(controlpoints)):
    # determine radial position of the controlpoint;
    # radialposition = np.sqrt(dot(controlpoints[icp].coordinates, controlpoints[icp].coordinates))
    radialposition = controlpoints[icp].radius
    u, v, w = 0, 0, 0 # initialize velocity
    # multiply icp line of Matrix with vector of circulation Gamma to calculate velocity at controlpoint
    for jring in range(0,len(rings1)):
      u += MatrixU[icp][jring]*Gamma1[jring] # axial component of velocity
      v += MatrixV[icp][jring]*Gamma1[jring] # y-component of velocity
      w += MatrixW[icp][jring]*Gamma1[jring] # z-component of velocity
    for jring in range(0,len(rings2)):
      u += MatrixU[icp][jring]*Gamma2[jring] # axial component of velocity
      v += MatrixV[icp][jring]*Gamma2[jring] # y-component of velocity
      w += MatrixW[icp][jring]*Gamma2[jring] # z-component of velocity
    # end for rings

    # calculate total perceived velocity
    vrot = cross3([-Omega*wakeSystem.direction, 0, 0], controlpoints[icp].coordinates) # rotational velocity
    vel1 = [wind[0] + u + vrot[0], wind[1] + v + vrot[1], wind[2] + w + vrot[2]] # total perceived velocity at section
    # calculate azimuthal and axial velocity
    azimdir = cross3([-1.0/radialposition, 0 , 0], controlpoints[icp].coordinates ) # rotational direction
    vazim = dot(azimdir , vel1) # azimuthal direction
    vaxial =  dot([1.0, 0, 0] , vel1) # axial velocity
    # calculate loads using blade element theory
    temploads = loadBladeElement(vaxial,vazim,radialposition,wakeSystem.rotor,wakeSystem.pitch,rho)
    # new point of new estimate of circulation for the blade section
    GammaNew[icp] = temploads[2]
    # update output vector
    b = controlpoints[icp].bladeID
    LLoutput.ax[b].append(u + vrot[0]) #-(u + vrot[0])/wind[0]
    LLoutput.az[b].append(vazim/(radialposition*Omega)-1)
    LLoutput.r[b].append(radialposition)
    LLoutput.Fnorm[b].append(temploads[0])
    LLoutput.Ftan[b].append(temploads[1])
    LLoutput.Gamma[b].append(temploads[2]) ## ?why? already have GammaNew
    LLoutput.local[b].append(temploads[3])
    LLoutput.torque[b].append(temploads[1]*radialposition)
    LLoutput.Cl[b].append(temploads[4])
    LLoutput.Cd[b].append(temploads[5])
  # end loop control points
  return GammaNew

def solveLiftingLineCoaxial(wakeSystemTop,wakeSystemBot,wind,RPM,rho=1.225):
    # this codes solves a lifting line model of a horizontal axis coaxial rotor
    """
      params:
        wakeSystemTop,wakeSystemBot: data structures that contain the geometry of the horseshoe vortex rings,
                                    and the control points along the blades, each blade must have unique bladeID
        wind: unperturbed windstream velocity, also known as U_infinity (3D vector)
        RPM: rotational velocity of the rotor rot/min
        rho: optional air density. assumed sea-level if not provided
    """
    ## get controlpoints data structure
    controlpointsTop = wakeSystemTop.controlPoints
    controlpointsBot = wakeSystemBot.controlPoints
    ## get horseshoe vortex rings data structure
    ringsTop = wakeSystemTop.rings
    ringsBot = wakeSystemBot.rings
    Omega = RPM*pi/30.0
    fwdVel = np.linalg.norm(wind)
    angleOfAttack = np.arctan2(np.linalg.norm(wind[1:]),wind[0])
    ## initialize variables that we will use during the calculation
    velocity_induced = [] # velocity induced by a horseshoe vortex ring at a control point (3 components u,v,w)
    radialposition, azimdir = 0, 0 # radial position of the control point
    alpha = 0 # angle of attack

    GammaNewTop = [0]*len(controlpointsTop) # new estimate of bound circulation, initialized as zeroes
    GammaTop = [0]*len(controlpointsTop) # current solution of bound circulation
    GammaNewBot = [0]*len(controlpointsBot) # new estimate of bound circulation, initialized as zeroes
    GammaBot = [0]*len(controlpointsBot) # current solution of bound circulation

    #vel1, vmag, vaxial, vazim, temploads = 0,0,0,0,0 # velocities and loads at controlpoint... no need to declare here
    MatrixUT = np.zeros((len(controlpointsTop),len(ringsTop))) # matrix of induction, for velocity component in x-direction
    MatrixVT = np.zeros((len(controlpointsTop),len(ringsTop))) # matrix of induction, for velocity component in y-direction
    MatrixWT = np.zeros((len(controlpointsTop),len(ringsTop))) # matrix of induction, for velocity component in z-direction
    MatrixUB = np.zeros((len(controlpointsBot),len(ringsBot))) # matrix of induction, for velocity component in x-direction
    MatrixVB = np.zeros((len(controlpointsBot),len(ringsBot))) # matrix of induction, for velocity component in y-direction
    MatrixWB = np.zeros((len(controlpointsBot),len(ringsBot))) # matrix of induction, for velocity component in z-direction

    # the variables below are to setup the maximum number of iterations and convergence criteria
    Niterations = 500
    errorlimit = 0.01
    error = 1.0
    refererror = 1.0
    ConvWeight = 0.2
    thrust = 0
    prevThrust = 0

    # initalize and calculate matrices for velocity induced by horseshoe vortex rings
    init_Vi_Matrices(controlpointsTop,ringsTop,MatrixUT,MatrixVT,MatrixWT)
    init_Vi_Matrices(controlpointsBot,ringsBot,MatrixUB,MatrixVB,MatrixWB)

    # calculate solution through an iterative process
    for kiter in range(1,Niterations+1):
      # copy the arrays, since arrays are mutable
      GammaTop = np.copy(GammaNewTop)
      GammaBot = np.copy(GammaNewBot)

      ## re-initialize output variables
      Nblades = wakeSystemTop.rotor.blades + wakeSystemBot.rotor.blades # wakeSystems' rotors should have unique bladeID, otherwise rotors will overwrite in the output!! FIX
      tempOut = LLCoaxOutput(Nblades)

      # calculate velocity, circulation and loads at the controlpoints, and update Gamma
      # also updates tempOut with loads/data for each blade
      GammaNewTop = calc_ControlPoints(wakeSystemTop,controlpointsTop,ringsTop,GammaTop,ringsBot,GammaBot,MatrixUT,MatrixVT,MatrixWT,Omega,wind,tempOut,rho=rho)
      GammaNewBot = calc_ControlPoints(wakeSystemBot,controlpointsBot,ringsBot,GammaBot,ringsTop,GammaTop,MatrixUB,MatrixVB,MatrixWB,Omega,wind,tempOut,rho=rho)

      #################################
      # check convergence of solution #
      refererrorT = np.max(np.absolute(GammaNewTop)) #gets the max Gamma value
      refererrorB = np.max(np.absolute(GammaNewBot))
      refererrorT = max(refererrorT,0.001) # define scale of bound circulation, 0.001 to make sure it isn't divide by 0
      refererrorB = max(refererrorB,0.001)
      diffT = np.average(np.absolute(np.subtract(GammaNewTop,GammaTop))) # difference betweeen iterations
      diffB = np.average(np.absolute(np.subtract(GammaNewBot,GammaBot)))
      error = max(diffT/refererrorT,diffB/refererrorB) # max relative error

      # use thrust to calculate new average induced velocity, and update wake system(s)
      prevThrust = thrust
      thrust = 0
      for i in range(len(tempOut.r)):
        thrust += trapz(tempOut.Fnorm[i],tempOut.r[i])

      tError = np.absolute(thrust-prevThrust)/thrust
      ## for checking convergence during code testing
      if kiter < 20:
        
        thrustT = trapz(tempOut.Fnorm[0],tempOut.r[0]) + trapz(tempOut.Fnorm[1],tempOut.r[1])
        thrustB = trapz(tempOut.Fnorm[2],tempOut.r[2]) + trapz(tempOut.Fnorm[3],tempOut.r[3])
        print(kiter,error,tError,(thrust-prevThrust)/thrust,thrustT,thrustB)
        #plt.plot(r_temp[0],Cd_temp[0])

      # check if GAMMA is sufficiently converged
      if error < errorlimit:
        # check to see if inflow is close enough
        # get the induced velocity/inflow ratio
        # this is the SINGLE rotor method... will need more complex method for coaxial!!! FIX THIS!
        Vi = wakeSystemTop.rotor.Vi_hover(thrust,rho)
        inflowRatio = wakeSystemTop.rotor.inflowRatio(fwdVel,angleOfAttack,Vi,RPM)
        if np.absolute(inflowRatio-wakeSystemTop.inflowRatio)/inflowRatio < 0.01 and (tError < errorlimit):
          # if error smaller than limit, stop iteration cycle
          # kiter < Niterations
          print(kiter,inflowRatio,wakeSystemTop.inflowRatio,(inflowRatio-wakeSystemTop.inflowRatio)/inflowRatio)
          break
        else:
          # update the wake model
          print(inflowRatio,wakeSystemTop.inflowRatio,kiter,error,tError,thrust,prevThrust)
          wakeSystemTop.update_Inflow(inflowRatio) # keeping track with the TOP rotor... bottom ignored until method is fixed

      # this seems to work better for all cases... correct to catch overshoots
      maxCW = 0.05
      minCW = 0.03
      ConvWeight = min(maxCW,error)
      ConvWeight = max(ConvWeight,minCW)
      # set new estimate of bound circulation
      for ig in range(0,len(GammaNewTop)):
        GammaNewTop[ig] = (1 - ConvWeight)*GammaTop[ig] + ConvWeight*GammaNewTop[ig]
      # end for
      for ig in range(0,len(GammaNewBot)):
        GammaNewBot[ig] = (1 - ConvWeight)*GammaBot[ig] + ConvWeight*GammaNewBot[ig]
      # end for
    # end iteration loop
    if kiter == Niterations: 
      print("max iterations!",kiter,error,tError,thrust,wakeSystemTop.inflowRatio)
    
    # get some total values, lift, torque, power...
    totalLift = 0
    totalTorque = 0
    totalPower = 0
    for i in range(len(tempOut.r)):
      totalLift += trapz(tempOut.Fnorm[i],tempOut.r[i])
      totalTorque += trapz(tempOut.torque[i],tempOut.r[i])
      totalPower += trapz(tempOut.torque[i],tempOut.r[i])*Omega

    #plt.show()
    #// output results of converged solution
    #return({"Vi":a_temp,"Vi_azim":aline_temp,"r":r_temp,"Fnorm":Fnorm_temp,"Ftan":Ftan_temp,"Gamma":Gamma_temp,"local":local_temp,"Power":totalPower,"Lift":totalLift,"Torque":totalTorque});
    return tempOut

def update_Gamma_single_ring(ring,GammaNew,WeightNew):
  for i in range(0,len(ring.filaments)):
    ring.filaments[i].gamma = ring.filaments[i].gamma*(1-WeightNew) + WeightNew * GammaNew
  # end for
  return ring 
# end update_Gamma_sinle_ring

def velocity_induced_single_ring(ring,controlpoint):
  """
    params:
      ring - list of filament location and vorticity strength
      controlpoint - 3D vector of control point
  """
  #tempvel1 = [0,0,0]
  vel_ind = [0,0,0]
  CORE = 0.00001
  for i in range(0,len(ring.filaments)):

    GAMMA = ring.filaments[i].gamma
    XV1 = ring.filaments[i].X1
    XV2 = ring.filaments[i].X2
    XVP = controlpoint.coordinates
    tempvel1 = velocity_3D_from_vortex_filament(GAMMA,XV1,XV2,XVP,CORE)
    # console.log(i);
    # console.log(tempvel1)
    # console.log("XV1 " + XV1 + " XV2 " + XV2 + " XVP " + XVP );
    # console.log(velind)
    # vel_ind = math.add(velind,tempvel1);

    vel_ind[0] += tempvel1[0]
    vel_ind[1] += tempvel1[1]
    vel_ind[2] += tempvel1[2]
    # print vel_ind
  # end for

  return vel_ind
# end velocity_induced_single_ring

def loadBladeElement(Vnorm,Vtan,r,rotor,pitch,rho=1.225):
  Vmag2 = (Vnorm**2 + Vtan**2)
  inflowAngle = atan(Vnorm/Vtan) # atan2 ???
  # console.log('inflow angle ' + InflowAngle)

  ## get the twist and chord of the blade at the desired radial location
  #temp = geoBlade(r_R) # returns [chord, twist+pitch] #### old code...
  chord = rotor.chord(r)#temp[0]
  twist = rotor.twist(r)#temp[1]

  ## calculate the local angle of attack of the blade, then lift and drag coefficient of airfoil
  # console.log('twist ' + twist)
  alpha = twist - inflowAngle + pitch # it was + ??? turbines ???

  # get Reynolds and Mach number for blade element
  Velem = np.sqrt(Vmag2)
  mach = Velem/343.2 # assume 20 C temp for now
  reynolds = Velem*chord*rho*10**6/18.27

  #Cl,Cd,Cm = rotor.airfoil.getCoef(alpha,reynolds,mach)

  # now we check for airfoil distribution or single airfoil
  if rotor.airfoil.type == "airfoil":
    # single airfoil is easy!
    Cl,Cd,Cm = rotor.airfoil.getCoef(alpha,reynolds,mach)
  elif rotor.airfoil.type == "airfoils":
    # find section of distribution
    loI = 0
    hiI = 0
    for xr in range(0,len(rotor.airfoil.distribution)-1):
      if r >= rotor.airfoil.distribution[xr] and r <= rotor.airfoil.distribution[xr+1]:
        loI = xr
        hiI = xr+1
        break
    # determine if the two airfoils are the same
    if loI == hiI or rotor.airfoil.aflist[loI] == rotor.airfoil.distribution[hiI]:
      Cl,Cd,Cm = rotor.airfoil.aflist[loI].getCoef(alpha,reynolds,mach)
    else: # if not, get coefficients from both and interpolate
      Cl1,Cd1,Cm1 = rotor.airfoil.aflist[loI].getCoef(alpha,reynolds,mach)
      Cl2,Cd2,Cm2 = rotor.airfoil.aflist[hiI].getCoef(alpha,reynolds,mach)
      hw = (r - rotor.airfoil.distribution[loI])/(rotor.airfoil.distribution[hiI]-rotor.airfoil.distribution[loI])
      lw = 1 - hw
      Cl = Cl1*lw + Cl2*hw
      Cd = Cd1*lw + Cd2*hw
      Cm = Cm1*lw + Cm2*hw

  ## get lift drag moment coefficients

  Lift = 0.5*Vmag2*Cl*chord*rho
  Drag = 0.5*Vmag2*Cd*chord*rho
  Fnorm = Lift*np.cos(inflowAngle)-Drag*np.sin(inflowAngle) ### CHECK SIGNS!!!
  Ftan = Lift*np.sin(inflowAngle)+Drag*np.cos(inflowAngle) ### CHECK SIGNS!!!
  Gamma = 0.5*Velem*Cl*chord*rho
  result = [Fnorm,Ftan,Gamma,alpha*180.0/pi,Cl,Cd,Cm]
  # console.log('Result ' + result)

  return result
# end loadBLadeElement

def cross3(vec1, vec2): #this is faster than numpy
    """ Calculate the cross product of two 3d vectors. """
    result = [0]*3
    a1, a2, a3 = float(vec1[0]), float(vec1[1]), float(vec1[2])
    b1, b2, b3 = float(vec2[0]), float(vec2[1]), float(vec2[2])
    result[0] = a2 * b3 - a3 * b2
    result[1] = a3 * b1 - a1 * b3
    result[2] = a1 * b2 - a2 * b1
    return result

def dot(vec1, vec2): #this is faster than numpy
    """ Calculate the cross product of two 3d vectors. """
    if (len(vec1) != len(vec2)): ## not matching vectors.
      return 0
    result = 0
    for i,n in enumerate(vec1):
      result += vec1[i]*vec2[i]
    return result

def getControlPoints(Npoints,hubWidth,cosDist=True):
  # Npoints number of control points along blade
  # hubWidth as a percentage of blade length
  # set number of control points, then spread in cosine, more lines near the tip/root
  controlPoints =  np.arange(0,Npoints+1,1)/float(Npoints)
  # cosine uniform distribution is default, else uniform distribution
  if cosDist: controlPoints = (-np.cos(np.pi*(controlPoints))+1)/2
  #adjust the control points to account for the hub
  controlPoints = controlPoints*(1-hubWidth)+hubWidth

  return controlPoints

def rotateYZ(y,z,angle):
  newY = y*np.cos(angle) - z*np.sin(angle)
  newZ = y*np.sin(angle) + z*np.cos(angle)
  return newY,newZ

class Filament:
  def __init__(self,x1,y1,z1,x2,y2,z2,theta,gamma=1):
    #declare a filament
    self.X1 = [x1,y1,z1]
    self.X2 = [x2,y2,z2]
    self.gamma = gamma
    self.theta = theta
  def update_Gamma(self,gammaNew=1,weightNew=1):
    self.gamma = self.gamma*(1 - weightNew) + weightNew * gammaNew
  def update_Inflow(self,inflowOld,inflowNew,rotorOffsetX):
    self.X1[0] = (self.X1[0] - rotorOffsetX)*inflowNew/inflowOld + rotorOffsetX
    self.X2[0] = (self.X2[0] - rotorOffsetX)*inflowNew/inflowOld + rotorOffsetX

class FilamentRing:
  # for now assume the flow comes from negative x/y direction
  def __init__(self,r1,r2,radius,inflowRatio,advanceRatio,PSIangle,blade,rotorOffsetX=0,rotorDirection=1,rotations=5,filPerRot=36):
    self.filaments = []
    self.blade = blade
    self.rotorX = rotorOffsetX
    self.rotorDirection = rotorDirection
    self.inflowRatio = inflowRatio
    self.advanceRatio = advanceRatio
    #self.vel = vel_vec
    # setup values for rotation matrix
    cosrot = cos(PSIangle)
    sinrot = sin(PSIangle)
    # determines number/resolution of filaments
    theta_list = np.arange(0,filPerRot*rotations+1,1)*(2*pi/filPerRot)
    # create bound filament
    bound = Filament(-self.rotorX,r1*cos(PSIangle),r1*sin(PSIangle),-self.rotorX,r2*cos(PSIangle),r2*sin(PSIangle),0) # gamma ommitted for now
    self.filaments.append(bound)
    x = -theta_list*inflowRatio*radius - self.rotorX
    y1 = np.cos(self.rotorDirection*theta_list)*r1
    z1 = np.sin(self.rotorDirection*theta_list)*r1
    y2 = np.cos(self.rotorDirection*theta_list)*r2
    z2 = np.sin(self.rotorDirection*theta_list)*r2

    Y1,Z1 = rotateYZ(y1,z1,PSIangle)
    Y2,Z2 = rotateYZ(y2,z2,PSIangle)

    Z1 = theta_list*advanceRatio*radius + Z1 
    Z2 = theta_list*advanceRatio*radius + Z2

    for i,t in enumerate(theta_list):
      if (i == 0): continue
      fil1 = Filament(x[i],Y1[i],Z1[i],x[i-1],Y1[i-1],Z1[i-1],t) # gamma ommitted for now
      fil2 = Filament(x[i-1],Y2[i-1],Z2[i-1],x[i],Y2[i],Z2[i],t)
      self.filaments.insert(0,fil1)
      self.filaments.append(fil2)
  def update_Gamma(self,gammaNew=1,weightNew=1):
    for fil in self.filaments:
      fil.update_Gamma(gammaNew,weightNew)
    return gammaNew*weightNew
  def update_Inflow(self,inflowNew):
    for fil in self.filaments:
      fil.update_Inflow(self.inflowRatio,inflowNew,self.rotorX)
    self.inflowRatio = inflowNew
    return self.inflowRatio

  def update_filaments(self,vel_vecNew):
    #assume the filaments are points on a circle, skewed by vel_vec
    #change the skew ratio by the ratio of vel_vec:vel_vecNew
    return
  def get_filament_points(self):
    points = [[],[],[]]
    for f in self.filaments:
      points[0].append(f.X1[0])
      points[1].append(f.X1[1])
      points[2].append(f.X1[2])
    #get last point
    points[0].append(f.X2[0])
    points[1].append(f.X2[1])
    points[2].append(f.X2[2])
    return points

class WakeSystem:
  def __init__(self,rotor,controlPointDist,inflowRatio,advanceRatio,pitch,PSIangle=0,rotorOffsetX=0,rotorDirection=1,startID=0):
    # Nblades - int - number of blades
    # radius - float - length of blade
    # controlPoints - list/array of points along blade, between 0 and 
    # chordDist, twistDist and pitch can be replaced by Rotor class eventually... FIX THIS
    self.rings = []
    self.controlPoints = []
    self.inflowRatio = inflowRatio
    self.advanceRatio = advanceRatio
    self.rotor = rotor
    self.pitch = pitch
    self.direction = rotorDirection
    for N in range(0,rotor.blades):
      cosrot = cos(N*2.0*pi/rotor.blades+PSIangle)
      sinrot = sin(N*2.0*pi/rotor.blades+PSIangle)
      for i,cp in enumerate(controlPointDist):
        if (i == 0): continue
        self.rings.append(FilamentRing(controlPointDist[i-1]*rotor.radius,controlPointDist[i]*rotor.radius,rotor.radius,inflowRatio,advanceRatio,N*2*pi/rotor.blades+PSIangle,N,rotorOffsetX,rotorDirection))
        r = rotor.radius*(controlPointDist[i-1]+controlPointDist[i])/2.0
        c = rotor.chord(r)
        t = rotor.twist(r)
        self.controlPoints.append(ControlPoint([-rotorOffsetX,r*cosrot,r*sinrot],r,c,[cos(t+pitch),-sin(t+pitch)*rotorDirection*sinrot,sin(t+pitch)*rotorDirection*cosrot],[-sin(t+pitch),-cos(t+pitch)*rotorDirection*sinrot,cos(t+pitch)*rotorDirection*cosrot],bladeID=N+startID))
  def update_Inflow(self,inflowNew):
    for ring in self.rings:
      ring.update_Inflow(inflowNew)
    self.inflowRatio = inflowNew
    return self.inflowRatio

def plotWakeSystems(wakeSystems,lims=(-2,2)):

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  colors = ["red","blue","green"]
  for WS in wakeSystems:

    for ring in WS.rings:
      p = ring.get_filament_points()
      ax.plot(p[0],p[1],p[2],color=colors[ring.blade])

    for cp in WS.controlPoints:
      x,y,z = cp.get_N_T_Vectors(0.5)
      ax.plot(x,y,z,color="black", alpha=0.5)

  ax.set_xlim3d(lims[0],lims[1])
  ax.set_ylim3d(lims[0],lims[1])
  ax.set_zlim3d(lims[0],lims[1])

class ControlPoint:
  #{coordinates: [ 0 , (span_array[i]+span_array[i+1])/2   , 0 ] , chord: 1, normal: [ 0,  0, 1] , tangential: [1, 0, 0]   } 
  def __init__(self,coord,radius,chord,normal=[0,0,1.0],tangential=[1.0,0,0],bladeID=0):
    self.coordinates = coord
    self.radius = radius # distance from control point to rotor center
    self.chord = chord # blade chord at location
    self.normal = normal # normal direction (lift)
    self.tangential = tangential # tangential direction (drag)
    self.bladeID = bladeID
  def get_N_T_Vectors(self,m=1.0):
    #return [self.coordinates[0]+self.normal[0]*m,self.coordinates[0]],[self.coordinates[1]+self.normal[1]*m,self.coordinates[1]],[self.coordinates[2]+self.normal[2]*m,self.coordinates[2]]
    return [self.coordinates[0]+self.normal[0]*m,self.coordinates[0],self.coordinates[0]+self.tangential[0]*m],[self.coordinates[1]+self.normal[1]*m,self.coordinates[1],self.coordinates[1]+self.tangential[1]*m],[self.coordinates[2]+self.normal[2]*m,self.coordinates[2],self.coordinates[2]+self.tangential[2]*m]


"""
def create_rotor_geometry(span_array, radius, tipspeedratio, Uinf, theta_array, nblades) {
  // create rotor eometry and rotor-wake circulation system
  var filaments = [];
  var ring = [];
  var controlpoints = [];
  var bladepanels = [];
  var geodef; var geodef2;
  var r; var angle;  var angle2 ; var dx; var dz; var dy; var dtheta; var xt; var yt; var zt;
  var temp1; var temp2;

  for (var krot = 0; krot < nblades; krot++) {

    var angle_rotation = 2*Math.PI/nblades*krot;
    var cosrot = Math.cos(angle_rotation);
    var sinrot = Math.sin(angle_rotation);

    for (var i = 0; i < span_array.length-1; i++) {
      r = (span_array[i]+span_array[i+1])/2;
      geodef = geoBlade(r/radius);
      angle = geodef[1]*Math.PI/180;
      // define controlpoints
      temp1 = {coordinates: [ 0 ,  r  , 0 ] , chord: geodef[0], normal: [ Math.cos(angle),  0, -1*Math.sin(angle)] , tangential: [-1*Math.sin(angle), 0, -1*Math.cos(angle)]   };
      // rotate blade to position
      temp1.coordinates = [ 0 ,  temp1.coordinates[1]*cosrot -  temp1.coordinates[2]*sinrot , temp1.coordinates[1]*sinrot +  temp1.coordinates[2]*cosrot ];
      temp1.normal = [ temp1.normal[0],  temp1.normal[1]*cosrot -  temp1.normal[2]*sinrot , temp1.normal[1]*sinrot +  temp1.normal[2]*cosrot ];
      temp1.tangential = [ temp1.tangential[0] ,  temp1.tangential[1]*cosrot -  temp1.tangential[2]*sinrot , temp1.tangential[1]*sinrot +  temp1.tangential[2]*cosrot ];

      controlpoints.push( temp1 );
      // define bound vortex filament
      temp1= {x1:0 , y1:span_array[i], z1:0, x2:0 , y2:span_array[i+1], z2:0, Gamma: 0  }   ;
      // rotate filament to position

      filaments.push(temp1);
      // create trailing filaments, at x1 of bound filament
      geodef = geoBlade(span_array[i]/radius);
      angle = geodef[1]*Math.PI/180;
      temp1= {x1: geodef[0]*Math.sin(-1*angle), y1:span_array[i], z1: -1*geodef[0]*Math.cos(angle), x2:0 , y2:span_array[i], z2:0, Gamma: 0  }   ;
      filaments.push(temp1);
      for (var j = 0; j < theta_array.length-1; j++) {
        xt = filaments[filaments.length-1].x1;
        yt = filaments[filaments.length-1].y1;
        zt = filaments[filaments.length-1].z1;
        dy = (Math.cos(-theta_array[j+1])-Math.cos(-theta_array[j])) * span_array[i];
        dz = (Math.sin(-theta_array[j+1])-Math.sin(-theta_array[j])) * span_array[i];
        dx = (theta_array[j+1]-theta_array[j])/tipspeedratio*radius;

        temp1= {x1: xt+dx, y1: yt+dy, z1: zt+dz, x2: xt , y2:yt, z2:zt, Gamma: 0  }   ;
        // rotate filament to position

        filaments.push(temp1);
      };

      // create trailing filaments, at x2 of bound filament
      geodef = geoBlade(span_array[i+1]/radius);
      angle = geodef[1]*Math.PI/180;
      temp1= {x1:0 , y1:span_array[i+1], z1: 0, x2:geodef[0]*Math.sin(-1*angle) , y2:span_array[i+1], z2:-1*geodef[0]*Math.cos(angle), Gamma: 0  }   ;
      filaments.push(temp1);
      for (var j = 0; j < theta_array.length-1; j++) {
        xt = filaments[filaments.length-1].x2;
        yt = filaments[filaments.length-1].y2;
        zt = filaments[filaments.length-1].z2;
        dy = (Math.cos(-theta_array[j+1])-Math.cos(-theta_array[j])) * span_array[i+1];
        dz = (Math.sin(-theta_array[j+1])-Math.sin(-theta_array[j])) * span_array[i+1];
        dx = (theta_array[j+1]-theta_array[j])/tipspeedratio*radius;

        temp1= {x1: xt, y1: yt, z1: zt, x2: xt+dx , y2:yt+dy, z2:zt+dz, Gamma: 0  }   ;
        // rotate filament to position

        filaments.push(temp1);
      };

      for (var ifil = 0; ifil < filaments.length; ifil++) {
        temp1=filaments[ifil];
        temp2 = [ temp1.y1*cosrot -  temp1.z1*sinrot , temp1.y1*sinrot +  temp1.z1*cosrot , temp1.y2*cosrot -  temp1.z2*sinrot , temp1.y2*sinrot +  temp1.z2*cosrot ];
        temp1.y1 = temp2[0];
        temp1.z1 = temp2[1];
        temp1.y2 = temp2[2];
        temp1.z2 = temp2[3];
        filaments[ifil] = temp1;
      }

      ring.push({filaments: filaments});
      filaments = [];

      // panel of the blade section
      geodef = geoBlade(span_array[i]/radius);
      angle = geodef[1]*Math.PI/180;
      geodef2 = geoBlade(span_array[i+1]/radius);
      angle2 = geodef2[1]*Math.PI/180;

      temp1= {
        p1: [-0.25*geodef[0]*Math.sin(-1*angle) , span_array[i], 0.25*geodef[0]*Math.cos(angle)],
        p2: [-0.25*geodef2[0]*Math.sin(-1*angle2) , span_array[i+1], 0.25*geodef2[0]*Math.cos(angle2)],
        p3: [0.75*geodef2[0]*Math.sin(-1*angle2) , span_array[i+1], -0.75*geodef2[0]*Math.cos(angle2)],
        p4: [0.75*geodef[0]*Math.sin(-1*angle) , span_array[i], -0.75*geodef[0]*Math.cos(angle)]
      };
      temp1.p1 = [ 0 ,  temp1.p1[1]*cosrot -  temp1.p1[2]*sinrot , temp1.p1[1]*sinrot +  temp1.p1[2]*cosrot ];
      temp1.p2 = [ 0 ,  temp1.p2[1]*cosrot -  temp1.p2[2]*sinrot , temp1.p2[1]*sinrot +  temp1.p2[2]*cosrot ];
      temp1.p3 = [ 0 ,  temp1.p3[1]*cosrot -  temp1.p3[2]*sinrot , temp1.p3[1]*sinrot +  temp1.p3[2]*cosrot ];
      temp1.p4 = [ 0 ,  temp1.p4[1]*cosrot -  temp1.p4[2]*sinrot , temp1.p4[1]*sinrot +  temp1.p4[2]*cosrot ];

      bladepanels.push(temp1);
    };

  };
  return({controlpoints: controlpoints ,  rings: ring, bladepanels:bladepanels});
};

#"""
"""
function create_straight_wing_geometry(span_array, Alpha, WakeLength){
  var temp1
  var temp2
  var filaments = [];
  var ring = [];
  var infinity =WakeLength;
  var controlpoints = [];


  for (var i = 0; i < span_array.length-1; i++) {
    controlpoints.push( {coordinates: [ 0 , (span_array[i]+span_array[i+1])/2   , 0 ] , chord: 1, normal: [ 0,  0, 1] , tangential: [1, 0, 0]   } );
    temp1= {x1: infinity , y1:span_array[i], z1: infinity*Math.sin(Alpha*Math.PI/180), x2: 1.25 , y2:span_array[i], z2:0, Gamma: 0  }   ;
    filaments.push(temp1);
    temp1= {x1: 1.25 , y1:span_array[i], z1:0, x2:0 , y2:span_array[i], z2:0, Gamma: 0  }   ;
    filaments.push(temp1);
    temp1= {x1:0 , y1:span_array[i], z1:0, x2:0 , y2:span_array[i+1], z2:0, Gamma: 0  }   ;
    filaments.push(temp1);
    temp1= {x1:0 , y1:span_array[i+1], z1:0, x2:1.25 , y2:span_array[i+1], z2:0, Gamma: 0  }   ;
    filaments.push(temp1);
    temp1= {x1:1.25 , y1:span_array[i+1], z1:0, x2:infinity , y2:span_array[i+1], z2:infinity*Math.sin(Alpha*Math.PI/180), Gamma: 0  }   ;
    filaments.push(temp1);
    ring.push({filaments: filaments});
    filaments = [];
  };
  return({controlpoints: controlpoints ,  rings: ring});
};



#"""


