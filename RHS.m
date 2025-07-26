function y = RHS(u, N, dx, vel, visc, eq_type, diff_scheme)

% Calcula as equações diferenciais

switch eq_type
    case 'wave'
        y = -vel * Diff(u, N, dx, diff_scheme);

    case 'diffusion'
        y = visc * Diff2(u, N, dx, 'centered4');  % Difusão sempre com centered4

    case 'burgers_nao_viscoso'
        flux = 0.5*u.^2;
        y = -Diff(flux, N, dx, diff_scheme);

    case 'burgers_viscoso'
        flux = 0.5*u.^2;
        conv = -Diff(flux, N, dx, diff_scheme);
        visc_term = visc * Diff2(u, N, dx, 'centered4');
        y = conv + visc_term;

    case 'kuramoto_sivashinsky'
        flux = 0.5*u.^2;
        transport = -Diff(flux, N, dx, diff_scheme);
        anti_diff = -Diff2(u, N, dx, 'centered4');
        hyper_diff = -visc * Diff4(u, N, dx, 'centered4');
        y = transport + anti_diff + hyper_diff;

    otherwise
        error('Equação não implementada!');
end
end
