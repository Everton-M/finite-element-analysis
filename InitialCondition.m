
function y = InitialCondition(Identifier, x, k)


switch Identifier
    case 'sine'
        y = sin(k*x);

    case 'gaussian'
        y = 1./(10*k*x.^2+1); % Curva aproximada de uma gaussiana

    case 'step'
        y = tanh(x*k)+1; % Curva aproximada de uma função degrau

    case 'doubleexp'
        y = exp(-((x+50)/10).^2) - exp(-((x-50)/10).^2); %Curva aproximada de uma dupla exponencial

    otherwise
        error(['Initial condition ', Identifier, ' not implemented'])
end
