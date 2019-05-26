function lambda_bot = get_lambda_bot(F,r,pitch,coaxial,flowfield,lambda_u)
%{
FORWARD Flight
This function is a simple code implementation of the equation found in
literature (see below) for the inflow distribution of the bottom rotor in two 
distinctive regions (inner - ingests top rotor downwash and outer - sees free-stream).

NOTE! Make sure that lambda bottom inner uses lambda_u from root to tip instead
of from root to rd - "Note that the contracting streamtubes must be mapped
from originating points on the upper rotor to receiving points on the
lower rotor within the inner area that is affected by the wake, i.e., a
streamtube emanating from radial point r on the upper rotor will map to
new radial point ar on the lower rotor" - Leishman

Inputs:
    F - (array) Prandtl tip loss function (not necessarily converged) 
    on the bottom rotor. 

    r - (array) of radial positions from dr to 1-dr

    pitch - (array) of pitch of each of the blade elements with respect to
    the rotor plane. This will be influenced by geometrical twist of the
    blade as well as by the collective setting.

    rotor - (struct) containing geometrical properties for both rotors such 
    as pitch distribution, radius, rpm, blade number, etc 

    flowfield - (struct) containing lambda_inf (array) for both rotors 
    (normalization is different for each rotor since tip speed may be different)

Outputs:
    lambda_bot - (array) of inflow distribution (not necessarily converged)
    along the span for the bottom rotor.

Other m-files required: none

MAT-files required: none

Literature referenced: 
    ! An optimum Coaxial Rotor System for Axial Flight. Leishman, 2008. See
    equation (13) & (14)

    Unmanned coaxial rotor helicopter dynamics and system parameter
    estimation. Rashid et al. Springer, 2014.
    
    Modelling and robust control of an unmanned coaxial rotor helicopter
    with unstructured uncertainties. Dong et al. Advances in Mechanical
    Engineering, 2017, Vol. 9(I) 1-14

Assumptions: 
    ! Viscosity/dissipation procedure: used to get similar looking plots to
    Leishman for bottom rotor. What I am doing is removing a part of the
    high induced velocities near the tip of the top rotor which simulates a
    dissipation effect between the two rotors.

    ! Wake contraction has to be assumed, contained in rd. Currently, the
    value used by Leishman of 82% is being used.

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
May 2019; Last revision: 26-May-2019
%}

%------------- BEGIN CODE --------------

rotor = coaxial.rotor;

lambda_inf = flowfield(2).lambda_inf;

rd =coaxial.params.rd;

rotor = rotor(2); %in case the generic rotor struct is given
sigma = rotor.solidity;
cl_a = rotor.aero.cl_alpha;
inner = r<rd; %inside the downwash circle condition (to select elements in array later)
outer = r>=rd; %outside the downwash circle condition (to select elements in array later)
%idx_inner = find(r<rd);


%Procedure has been verified in the command window. lambda_u now has the
%same shape but with a coarser resolution to match the dimensions of
%r(inner).

%% sensitive area - instead of defining ri from 0 to 1 translate exactly the lambda up plot (starts at r(1) and ends at r(end))
interp_length = length(r(inner)); 
visc = 1.2; %equivalent to dissipation in my eyes
ri = linspace(r(1),r(end)/visc,interp_length);
lambda_u = interp1(r,lambda_u,ri);

%test - weird plots happening 
%lambda_u = lambda_u(inner);

%%
%since lambda inf is uniform, it can just be truncated so lambda_inf(inner)

lambda_bot_inner = sqrt((sigma*cl_a*1./(16*F(inner))-lambda_u/(2*rd^2)-lambda_inf(inner)/2).^2+...
    sigma*cl_a*pitch(inner).*r(inner)*1./(8*F(inner)))-sigma*cl_a*1./(16*F(inner))+lambda_u/(2*rd^2)+lambda_inf(inner)/2;

lambda_bot_outer = sqrt((sigma*cl_a*1./(16*F(outer))-lambda_inf(outer)/2).^2+sigma*cl_a*pitch(outer).*r(outer)*1./(8*F(outer)))-sigma*cl_a*1./(16*F(outer))+lambda_inf(outer)/2;

lambda_bot = [lambda_bot_inner lambda_bot_outer];

end