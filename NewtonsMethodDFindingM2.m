function [M2] = NewtonsMethodDFindingM2(M1, M2_init, gamma, A1, A2)

[D1, ~, ~, ~, ~, ~] = MachNumberFunctions(M1, gamma);
error = 1; %Initialize error value
i = 1; %Iteration Number

while abs(error) >= 1e-10
    if i == 1
        M2_k(i) = M2_init;
    else
        M2_k(i) = M2_kPlus1(i-1);
    end 
    [D_k(i), dDOverdM_k(i), ~, ~, ~, ~] = MachNumberFunctions(M2_k(i), gamma);
    %M2_k(i) = M2_k;
    M2_kPlus1(i) = M2_k(i) + ((A1.*D1)./A2 - D_k(i))./(dDOverdM_k(i)); %Calculate the next one
    
    error(i) = M2_kPlus1(i) - M2_k(i); %Calculate the error between the guess and calculated

    i = i + 1;
end
M2 = M2_kPlus1(end);
end