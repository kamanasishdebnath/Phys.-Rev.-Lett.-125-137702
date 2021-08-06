function dydt = equations(t, y, g, N, wc, delta, kappa, Gamma, num_of_ensemble, pulse1, pulse2, T1, T2, T3, T4)


number_of_equations= 1 + 2*(num_of_ensemble);
dydt= zeros(number_of_equations, 1);


if t>=T1 && t<=T2
    dydt(1)= -1j*wc*y(1) - 0.5*kappa*y(1) - 1j*pulse1;
    
elseif t>=T3 && t<=T4
    dydt(1)= -1j*wc*y(1) - 0.5*kappa*y(1) - 1j*pulse2;
    
else
    dydt(1)= -1j*wc*y(1) - 0.5*kappa*y(1);
    
end

for i= 0:(num_of_ensemble-1)
    index= 2*i + 2;
    wa= delta(i+1);
    coup= g(i+1);
   
    dydt(1)= dydt(1) - 1j*coup*N*y(index);
    dydt(index)= -2j*wa*y(index) + 1j*coup*y(1)*y(index+1) - (Gamma)*y(index);
    dydt(index+1)= -2j*coup*y(1)*conj(y(index)) + 2j*coup*y(index)*conj(y(1));
    
end

end
