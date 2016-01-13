function omega_list = evaluate_omega(t_list, params)

cntr = 1;

for t = t_list
    tspan = [params.start_t t]
    y0 = [params.S_init params.I_init (params.N - params.S_init - params.I_init)];
    [t_vals, y_vals] = ode45(@(t,y)simulateEpiModel(t,y,params), tspan, y0);

    %plot(t_vals, y_vals); hold on;
    
    t_end = t_vals(end)
    %S_end = y_vals(end, 1);
    
    I_end = y_vals(end, 2);
    
    omega_t = params.k*I_end
   
    omega_list(cntr) = omega_t;
    %figure; plot(t_end, omega_t, '.'); hold on;
    cntr = cntr + 1;
    
    t_vals = [];
    y_vals = [];
   
end

%figure; 
plot(t_list, omega_list, '.'); hold on;

t_list
omega_list
%pause