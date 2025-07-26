% Configurar variáveis básicas necessárias antes da simulação, que são
% x (grade), u (matriz que armazenará a saída do solucionador) e parâmetros de simulação
% derivados das entradas do usuário (delta x e delta t)

function [x, u, dx, dt] = Setup(Inputs)

% Cria grade
x = linspace(Inputs.x0,Inputs.xn,Inputs.SizeX);

% Tamanho de delta x
dx = ( Inputs.xn - Inputs.x0 ) / ( Inputs.SizeX - 1 );

% Tamanho de delta t
if (strcmp(Inputs.RHSIdentifier,'diffusion'))
   dt = Inputs.CFLViscous * dx^2 / abs(Inputs.Viscosity); % Viscous CFL
else
    dt = Inputs.CFL * dx / abs(Inputs.Velocity); % Convective CFL
end

% Inicializa o array u que conterá a saída do solucionador (matriz nx por nt).
% Ou seja, um snapshot corresponde a uma coluna

u = zeros(Inputs.SizeX,Inputs.SizeT);

end

