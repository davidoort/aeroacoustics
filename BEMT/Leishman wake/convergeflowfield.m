function [Fcf_u, lambda_u, Fcf_l, lambda_l] = convergeflowfield(flowfield, r, epsilon, rotor)

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
end

%% Intermezzo
lambda_u = lambda_tot_u-flowfield(1).lambda_inf;

lambda0_l = zeros(1,length(r)); %dummy, to start iteration
lambda_tot_l = get_lambda_bot(Fcf_l,r,pitch_l,rotor,flowfield,lambda_u);
%% converge bottom rotor 
i = 0;
while norm(Fcf_l-Fcf0_l)>epsilon || norm(lambda_tot_l-lambda0_l)>epsilon
    Fcf0_l = Fcf_l;
    lambda0_l = lambda_tot_l;
    
    Fcf_l = Prandtl_tip_loss(r,lambda_tot_l,rotor(2));
    lambda_tot_l = get_lambda_bot(Fcf_l,r,pitch_l,rotor,flowfield,lambda_u);
    i = i+1;
end
    
lambda_l = lambda_tot_l-flowfield(2).lambda_inf; %From Leishman paper
end