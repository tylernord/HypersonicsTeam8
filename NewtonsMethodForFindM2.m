function [M2] = NewtonsMethodForFindM2(M1, M2_init, gamma)

%Find N1
[~, ~, ~, ~, N1, ~] = MachNumberFunctions(M1, gamma);
error = 1; %Initialize error value
i = 1; %Iteration Number

while abs(error) >= 1e-10
    if i == 1
        M2_k(i) = M2_init;
    else
        M2_k(i) = M2_kPlus1(i-1);
    end 
    %Find N2 and its derivative
    [~, ~, ~, ~, N_k(i), dNOverdM_k(i)] = MachNumberFunctions(M2_k(i), gamma);
    %M2_k(i) = M2_k;
    M2_kPlus1(i) = M2_k(i) + (N1 - N_k(i))./(dNOverdM_k(i)); %Calculate the next one
    
    error(i) = M2_kPlus1(i) - M2_k(i); %Calculate the error between the guess and calculated
    i = i+1;

end
M2 = M2_kPlus1(end);