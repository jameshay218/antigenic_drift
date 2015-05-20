x_eq1 = x_eq;
x_eq = [];
x_eq(1:100) = x_eq1(1:100);
x_eq(101:200) = x_eq1(301:400);
x_eq(201:300) = x_eq1(601:700);
save('dat_x_eq_ini.mat', 'x_eq');