function d4u = Diff4(u, SizeU, dx, scheme)
% Calcula a 4a derivada com esquemas centrados (2a ou 4a ordem)

d4u = zeros(1, SizeU);

switch scheme
    case 'centered' % Esquema 2a ordem
        for i = 3:SizeU-2
            d4u(i) = (u(i+2) - 4*u(i+1) + 6*u(i) - 4*u(i-1) + u(i-2))/dx^4;
        end
        % Bordas periódicas
        d4u(1) = (u(3) - 4*u(2) + 6*u(1) - 4*u(end-1) + u(end-2))/dx^4;
        d4u(2) = (u(4) - 4*u(3) + 6*u(2) - 4*u(1) + u(end-1))/dx^4;
        d4u(end-1) = (u(2) - 4*u(end) + 6*u(end-1) - 4*u(end-2) + u(end-3))/dx^4;
        d4u(end) = d4u(1);

    case 'centered4' % Esquema 4a ordem
        for i = 4:SizeU-3
            d4u(i) = (-u(i+3) + 12*u(i+2) - 39*u(i+1) + 56*u(i) ...
                     -39*u(i-1) + 12*u(i-2) - u(i-3))/(6*dx^4);
        end
        % Bordas periódicas
        d4u(1) = (-u(4) + 12*u(3) - 39*u(2) + 56*u(1) ...
                 -39*u(end-1) + 12*u(end-2) - u(end-3))/(6*dx^4);
        d4u(2) = (-u(5) + 12*u(4) - 39*u(3) + 56*u(2) ...
                 -39*u(1) + 12*u(end-1) - u(end-2))/(6*dx^4);
        d4u(3) = (-u(6) + 12*u(5) - 39*u(4) + 56*u(3) ...
                 -39*u(2) + 12*u(1) - u(end-1))/(6*dx^4);
        d4u(end-2) = (-u(2) + 12*u(end) - 39*u(end-1) + 56*u(end-2) ...
                     -39*u(end-3) + 12*u(end-4) - u(end-5))/(6*dx^4);
        d4u(end-1) = (-u(3) + 12*u(2) - 39*u(end) + 56*u(end-1) ...
                     -39*u(end-2) + 12*u(end-3) - u(end-4))/(6*dx^4);
        d4u(end) = d4u(1);

    otherwise
        error('Differentiation scheme not implemented! Use "centered" or "centered4"');
end
end
