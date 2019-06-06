function [Fcf_u, lambda_u, Fcf_l, lambda_l] = convergeflowfield(flowfield, r, epsilon, coaxial)
%{
AXIAL Flight
convergeflowfield - This function solves for the Prandtl tip loss function
along the span of the blade on the upper and lower rotor as well as the
inflow ratio lambda on the upper and lower rotor. It is called
"converge" because there is an iterative component of getting to the
results since lambda and F mutually depend on each other.

NOTE: Usually converges within 5 iterations.

Inputs:
    flowfield - (struct) containing lambda_inf (array) for both rotors 
    (normalization is different for each rotor since tip speed may be different)

    r - (array) of radial positions from dr to 1-dr

    epsilon - accuracy of convergence (scalar)

    coaxial - (struct object) containing geometrical properties for both rotors such 
    as pitch distribution, radius, rpm, blade number, etc 

Outputs:
    Fcf_u - (array) containing the converged (for a given collective pitch
    setting) Prandtl tip loss function at every blade span for upper rotor

    lambda_u - (array) with converged inflow distribution as a function of
    blade span for the upper rotor

    Fcf_l - (array) containing the converged (for a given collective pitch
    setting) Prandtl tip loss function at every blade span for upper rotor

    lambda_l - (array) with converged inflow distribution as a function of
    blade span for the upper rotor

Other m-files required (in path): 
    get_lambda_up
    Prandtl_tip_loss
    get_lambda_bot

MAT-files required: none

Literature referenced: 
    ! An optimum Coaxial Rotor System for Axial Flight. Leishman, 2008.

    Unmanned coaxial rotor helicopter dynamics and system parameter
    estimation. Rashid et al. Springer, 2014.
    
    Modelling and robust control of an unmanned coaxial rotor helicopter
    with unstructured uncertainties. Dong et al. Advances in Mechanical
    Engineering, 2017, Vol. 9(I) 1-14

Assumptions:
    ! Top rotor not influenced by bottom rotor: its inflow can be converged
    independently.

    Start at F=1 provides fast convergence.

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
March 2019; Last revision: 21-April-2019
%}

%------------- BEGIN CODE --------------


rotor = coaxial.rotor;

pitch_u = rotor(1).pitch;
pitch_l = rotor(2).pitch;


Fcf0_u = zeros(1,length(r)); %dummy, to start iteration
Fcf0_l = zeros(1,length(r)); %dummy, to start iteration
Fcf_u = ones(1,length(r)); %"Fcf started with an initial value of 1" - in standard procedure also by Leishman
Fcf_l = ones(1,length(r)); %"Fcf started with an initial value of 1"

lambda0_u = zeros(1,length(r)); %dummy, to start iteration

lambda_tot_u = get_lambda_up(Fcf_u,r,pitch_u,rotor,flowfield); %remember that this is lambda_u_induced+lambda_inf


%% converge upper rotor - can be done without converging bottom rotor (assumption of no influence of bottom on top)
i = 0;
while norm(Fcf_u-Fcf0_u)>epsilon || norm(lambda_tot_u-lambda0_u)>epsilon
    Fcf0_u = Fcf_u;
    lambda0_u = lambda_tot_u;
    
    Fcf_u = Prandtl_tip_loss(r,lambda_tot_u,rotor(1));
    lambda_tot_u = get_lambda_up(Fcf_u,r,pitch_u,rotor,flowfield);
    i = i+1;
    if i > 1e5
        error("Couldn't converge flowfield upper rotor")
    end
end

%% Intermezzo
lambda_u = lambda_tot_u-flowfield(1).lambda_inf;

lambda0_l = zeros(1,length(r)); %dummy, to start iteration
lambda_tot_l = get_lambda_bot(Fcf_l,r,pitch_l,coaxial,flowfield,lambda_u);
%% converge bottom rotor 
i = 0;
while norm(Fcf_l-Fcf0_l)>epsilon || norm(lambda_tot_l-lambda0_l)>epsilon
    Fcf0_l = Fcf_l;
    lambda0_l = lambda_tot_l;
    
    Fcf_l = Prandtl_tip_loss(r,lambda_tot_l,rotor(2));
    lambda_tot_l = get_lambda_bot(Fcf_l,r,pitch_l,coaxial,flowfield,lambda_u);
    i = i+1;
    if i > 1e3
        error("Couldn't converge flowfield lower rotor")
    end
end
    
lambda_l = lambda_tot_l-flowfield(2).lambda_inf; %From Leishman paper

%------------- END OF CODE --------------

end