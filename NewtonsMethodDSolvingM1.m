function [M1] = NewtonsMethodDSolvingM1(M2, M1_init, gamma, A1, A2)
%Used for finding the previous mach number use the mach number function, D.

[D2, ~, ~, ~, ~, ~] = MachNumberFunctions(M2, gamma);
error = 1; %Initialize error value
i = 1; %Iteration Number

while abs(error) >= 1e-10
    if i == 1
        M1_k(i) = M1_init;
    else
        M1_k(i) = M1_kPlus1(i-1);
    end 
    [D_k(i), dDOverdM_k(i), ~, ~, ~, ~] = MachNumberFunction(M1_k(i), gamma);
    
    M1_kPlus1(i) = M1_k(i) + ((A2.*D2)./A1 - D_k(i))./(dDOverdM_k(i)); %Calculate the next one
    
    error(i) = M1_kPlus1(i) - M1_k(i); %Calculate the error between the guess and calculated
    i = i+1;

end
M1 = M1_kPlus1(end);