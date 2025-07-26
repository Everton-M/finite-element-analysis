function uNext = TimeMarch(u, N, dx, dt, vel, visc, method, eq_type, diff_scheme)
% Avança a solução no tempo usando o esquema especificado

switch method
    case 'euler' % Euler explícito (1 primeira ordem)
        uNext = u + dt * RHS(u, N, dx, vel, visc, eq_type, diff_scheme);

    case 'rk2' % Runge-Kutta 2 segunda ordem
        k1 = RHS(u, N, dx, vel, visc, eq_type, diff_scheme);
        u_temp = u + dt * k1;
        k2 = RHS(u_temp, N, dx, vel, visc, eq_type, diff_scheme);
        uNext = u + 0.5*dt*(k1 + k2);

    case 'rk4' % Runge-Kutta clássico (4 quarta ordem)
        k1 = RHS(u, N, dx, vel, visc, eq_type, diff_scheme);
        k2 = RHS(u + 0.5*dt*k1, N, dx, vel, visc, eq_type, diff_scheme);
        k3 = RHS(u + 0.5*dt*k2, N, dx, vel, visc, eq_type, diff_scheme);
        k4 = RHS(u + dt*k3, N, dx, vel, visc, eq_type, diff_scheme);
        uNext = u + dt*(k1 + 2*k2 + 2*k3 + k4)/6;

    case 'maccormack' % Preditor-corretor
        % Preditor (backward)
        u_pred = u + dt * RHS(u, N, dx, vel, visc, eq_type, 'backward');
        % Corretor (forward)
        uNext = 0.5*(u + u_pred + dt * RHS(u_pred, N, dx, vel, visc, eq_type, 'forward'));

    otherwise
        error('Esquema temporal não implementado!');
end
end
