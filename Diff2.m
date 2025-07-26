function ddu = Diff2(u, SizeU, dx, DifferentiationSchemeIdentifier)

% Inicializa du
ddu = zeros(1, SizeU);

switch DifferentiationSchemeIdentifier

    case 'centered'

        % Calcula derivadas na qual nõ precisam de condição de contorno
        for i = 2:SizeU-1
          ddu(i) = (u(i+1)- 2*u(i)+u(i-1))/(dx)^2;
        end

        % Condições de contorno: ddu = 0 no contorno
         ddu(1) = (u(2) - 2*u(1) + u(SizeU-1)) / dx^2;
         ddu(SizeU) = ddu(1);  % Periodicidade: u(SizeU+1) ? u(1)

    case 'centered4'
            % Ordem de precisão: O(?x^4)
            for i = 3:SizeU-2
                ddu(i) = (-u(i+2) + 16*u(i+1) - 30*u(i) + 16*u(i-1) - u(i-2)) / (12*dx^2);
            end
            % Condições de contorno periodicas
            ddu(1) = (-u(3) + 16*u(2) - 30*u(1) + 16*u(SizeU-1) - u(SizeU-2)) / (12*dx^2);
            ddu(2) = (-u(4) + 16*u(3) - 30*u(2) + 16*u(1) - u(SizeU-1)) / (12*dx^2);
            ddu(SizeU-1) = (-u(2) + 16*u(SizeU) - 30*u(SizeU-1) + 16*u(SizeU-2) - u(SizeU-3)) / (12*dx^2);
            ddu(SizeU)   = ddu(1);

    % Caso o usuario não escolha um esquema de marcha no tempo implementado
    otherwise
        error(['Differentiation scheme identified as ', ...
            DifferentiationSchemeIdentifier, ' not implemented! (Second Derivative)']);



end
