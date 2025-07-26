% Retorna a derivada de u calculada usando um esquema numérico específicado

function du = Diff(u, SizeU, dx, DifferentiationSchemeIdentifier)

% Inicializa du
du = zeros(1, SizeU);

switch DifferentiationSchemeIdentifier
    case 'backward'

         % Esquema retrógrado (backward) de 1ª ordem
        for i = 2:SizeU
            du(i) = (u(i) - u(i-1)) / dx;
        end
        % Condição periódica no primeiro ponto (i=1): u(0) ≡ u(SizeU)
        du(1) = du(SizeU);

    case 'centered'

        % Calculate derivatives which don't need boundary conditions
        for i = 2:SizeU-1
           du(i) = ( u(i+1) - u(i-1) ) / ( 2*dx ) ;
        end

        % Calcula a derivada no contorno assumindo contorno periódico
        % Condições.  u0 = uN - 1, e u1 = uN
        du(1) = ( u(2) - u(SizeU - 1) ) / ( 2*dx );

        du(SizeU) = du(1); % Lembrar que: u1 = uN


    case 'centered4'
        %Esquema centrado de quarta ordem
        for i = 3:SizeU-2
            du(i) = ( -u(i+2) + 8*u(i+1) - 8*u(i-1) + u(i-2) ) / (12*dx);
        end

        % Condições de contornos periodicas.
        du(1) = ( -u(3) + 8*u(2) - 8*u(SizeU-1) + u(SizeU-2) ) / (12*dx);
        du(2) = ( -u(4) + 8*u(3) - 8*u(1) + u(SizeU-1) ) / (12*dx);
        du(SizeU-1) = ( -u(2) + 8*u(SizeU) - 8*u(SizeU-2) + u(SizeU-3) ) / (12*dx); %Não sei se é 1 ou 2
        du(SizeU)   = du(1);  % Periodicidade: u(1) ≡ u(SizeU+1)



    case 'forward'

        % Esquema avançado (forward) de 1ª ordem
        for i = 1:SizeU-1
            du(i) = (u(i+1) - u(i)) / dx;
        end
        % Condição periódica no último ponto (i=SizeU): u(SizeU+1) ≡ u(1)
        du(1) = du(SizeU);

    case 'backward2'

        % Esquema backward de 2ª ordem
        for i = 3:SizeU
            du(i) = (3*u(i) - 4*u(i-1) + u(i-2)) / (2*dx);
        end

        % Condições periódicas para os dois primeiros pontos
        du(2) = (3*u(2) - 4*u(1) + u(SizeU-1)) / (2*dx);
        du(1) = (3*u(1) - 4*u(SizeU-1) + u(SizeU-2)) / (2*dx);

    % Caso o usuário não escolha um esquema de marcha no tempo implementado
    otherwise
        error(['Differentiation scheme identified as ', ...
            TimeMarchingSchemeIdentifier,  ' not implemented!']);



end
