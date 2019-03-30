function [] = convergeflowfield()

Fcf0_u = zeros(1,length(r)); %dummy, to start iteration
Fcf0_l = zeros(1,length(r)); %dummy, to start iteration
Fcf_u = ones(1,length(r)); %"Fcf started with an initial value of 1"
Fcf_l = ones(1,length(r)); %"Fcf started with an initial value of 1"

lambda0_u = zeros(1,length(r)); %dummy, to start iteration
lambda0_l = zeros(1,length(r)); %dummy, to start iteration
lambda_u = get_lambda_up(Fcf_u,r,pitch_u,rotor);
lambda_l = get_lambda_bot(Fcf_l,r,pitch_l,rotor);


%converge upper rotor - can be done without converging bottom rotor (assumption of no influence of bottom on top)
i = 0;
while norm(Fcf_u-Fcf0_u)>epsilon || norm(lambda_u-lambda0_u)>epsilon
    Fcf0_u = Fcf_u;
    lambda0_u = lambda_u;
    
    Fcf_u = Prandtl_tip_loss(r,lambda_u,rotor(1));
    lambda_u = get_lambda_up(Fcf_u,r,pitch_u,rotor);
    i = i+1;
end

%converge bottom rotor 
i = 0;
while norm(Fcf_l-Fcf0_l)>epsilon || norm(lambda_l-lambda0_l)>epsilon
    Fcf0_l = Fcf_l;
    lambda0_l = lambda_l;
    
    Fcf_l = Prandtl_tip_loss(r,lambda_l,rotor(2));
    lambda_l = get_lambda_bot(Fcf_l,r,pitch_l,rotor);
    i = i+1;
end
    

end