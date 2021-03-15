function dae = linearTangentDae(sol);

a = 100;

t = sol.time;
x = sol.state;
u = sol.control;
x1dot = x(:,3);
x2dot = x(:,4);
x3dot = a.*u(:,1);
x4dot = a.*u(:,2);
path = u(:,1).^2+u(:,2).^2;

dae = [x1dot x2dot x3dot x4dot path];

%----------------------------------%
% END: function linearTangentDae.m %
%----------------------------------%
