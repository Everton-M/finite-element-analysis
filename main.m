clc; clear all; close all;

%% Configurações da Simulação
% --- Controle principal da simulação ---

% Domínio espacial
Inputs.x0 = -100;
Inputs.xn =  100;
Inputs.SizeX = 201; % Número de pontos na malha

% Tempo de simulação
Inputs.SizeT = 200; % Número de iterações temporais

% Parâmetros Físicos
Inputs.Velocity = 1; % Velocidade de propagação
Inputs.Viscosity = 1; % Coeficiente de viscosidade

% Controle do passo de tempo
Inputs.CFL = 0.1; % Número CFL para termos convectivos
Inputs.CFLViscous = 0.1; % CFL para termos viscosos
Inputs.CFL_hyperviscous = 0.1; % CFL para hiperviscosidade

% Configurações da simulação
Inputs.k = 1; % Parâmetro de onda inicial
Inputs.InitialConditionIdentifier = 'sine'; % Opções 'gaussian', 'sine', 'step', 'doubleexp'
Inputs.RHSIdentifier = 'burgers_viscoso'; % 'diffusion,'burgers_nao_viscoso', 'burgers_viscoso', 'kuramoto_sivashinsky' % Equação a ser resolvida
Inputs.TimeMarchingSchemeIdentifier = 'rk4'; % 'rk2', 'rk4', 'euler', 'maccormack'
Inputs.DifferentiationSchemeIdentifier = 'centered4'; % 'centered', 'centered4', 'backward', 'backward2', 'forward', etc

%% Preparação da Simulação
[x, u, dx] = Setup(Inputs);
u(:,1) = InitialCondition(Inputs.InitialConditionIdentifier, x, Inputs.k);

% Vetor para armazenar o tempo
time_vector = zeros(1, Inputs.SizeT);

%% Execução da Simulação
for n = 1:Inputs.SizeT-1
    % Cálculo dinâmico do passo de tempo para equações não-lineares
    isNonLinear = strcmp(Inputs.RHSIdentifier, 'burgers_nao_viscoso') || ...
                  strcmp(Inputs.RHSIdentifier, 'burgers_viscoso') || ...
                  strcmp(Inputs.RHSIdentifier, 'kuramoto_sivashinsky');

    if isNonLinear
        u_max = max(abs(u(:,n)));
        if u_max < 1e-9, u_max = 1e-9; end

        dt_conv = Inputs.CFL * dx / u_max;
        dt_visc = Inputs.CFLViscous * dx^2 / Inputs.Viscosity;
        dt_hyper = Inputs.CFL_hyperviscous * dx^4 / Inputs.Viscosity;

        dt_final_para_passo = min([dt_conv, dt_visc, dt_hyper]);
    else
        dt_final_para_passo = dt;
    end

    time_vector(n+1) = time_vector(n) + dt_final_para_passo;

    % Avanço no tempo
    u(:,n+1) = TimeMarch(u(:,n)', Inputs.SizeX, dx, dt_final_para_passo, Inputs.Velocity, ...
                     Inputs.Viscosity, Inputs.TimeMarchingSchemeIdentifier, ...
                     Inputs.RHSIdentifier, Inputs.DifferentiationSchemeIdentifier);

    % Visualização durante a simulação
    if mod(n,20) == 0
        %{figure(1)
        %plot(x, u(:,n), 'b-', 'LineWidth', 1.0)
        %ylim([min(u(:,1)), max(u(:,1))])
        %title_text = sprintf('Iteração: %d', n);
        %title(title_text);
        %pause(0.01)

        figure(1)
        plot(x, u(:,n), 'r-', 'LineWidth', 0.5)
        ylim([min(u(:,1)), max(u(:,1))])

        % Construção do título
        title_text = sprintf('Iteração: %d', n);
        title(title_text);

        % Informações detalhadas na legenda
        legend_text = {
            ['Equação: ', Inputs.RHSIdentifier ], ...
            ['Cond. Inicial: ', 'Testeeee'], ...
            ['Marcha no tempo: ', '.'], ...
            ['Derivada: ', '.']
        };
        legend(legend_text, 'Location', 'northeast', 'FontSize', 12);

        % Texto com as informações
        info_text2 = sprintf(['Equação: %s\nCond. Inicial: %s\nMarcha no tempo: %s\nEsq. de diferenciação: %s'], ...
            Inputs.RHSIdentifier, Inputs.InitialConditionIdentifier, ...
            Inputs.TimeMarchingSchemeIdentifier, Inputs.DifferentiationSchemeIdentifier);

        % Adiciona a caixa de anotação no canto ...
        annotation('textbox', [0.02, 0.02, 0.02, 0.02], ...  % [x y width height] (valores entre 0 e 1, frações da figura)
            'String', info_text2, ...
            'FitBoxToText', 'on', ...
            'BackgroundColor', 'white', ...
            'EdgeColor', 'black', ...
            'FontSize', 13, ...
            'HorizontalAlignment', 'center', ...
            'LineWidth', 0.5, ...
            'Units', 'normalized', ...
            'Tag', 'InfoBox'
            );

        pause(0.01)


    end
end

%% Plot Final
% Figura 1 - Comparacao entre condicao inicial e final

f = figure(1);
set(gcf, 'Color', 'w', 'Position', [100, 100, 1200, 500]);

% Curvas: condicao inicial (preta tracejada) e final (azul solida)
hold on
plot(x, u(:,1), 'k--', 'LineWidth', 1, 'DisplayName', ['Condicao Inicial: ', Inputs.InitialConditionIdentifier ]);
plot(x, u(:,end), 'b-', 'LineWidth', 0.6, 'DisplayName', 'Estado Final');
hold off

% Titulo com nome da equacao (sem "EQ")
title_text = sprintf('%s | %s | %s', ...
    replace(Inputs.RHSIdentifier, '_', ' '), ...
    upper(Inputs.TimeMarchingSchemeIdentifier), ...
    upper(Inputs.DifferentiationSchemeIdentifier));
title(title_text, 'FontSize', 14);

xlabel('x');
ylabel('u');
xlim([Inputs.x0, Inputs.xn]);
grid on;
legend('show', 'Location', 'best');


%% Figura 2 - Mapa Espaco-Tempo

f2 = figure(2);
set(gcf, 'Color', 'w', 'Position', [150, 150, 1000, 600]);

% Plot: espaco no eixo X, tempo no eixo Y
pcolor(x, time_vector, u'); % transposta para alinhar dimensoes
shading interp;
colorbar;

n = 160; % pontos por segmento (para total ~160)

custom_colormap = [
    linspace(0,0.4,n)', linspace(0,0.7,n)', linspace(0.5,1,n)';   % azul escuro -> azul claro
    linspace(0.4,1,n)', linspace(0.7,1,n)', linspace(1,0,n)';     % azul claro -> amarelo
    linspace(1,1,n)', linspace(1,0,n)', linspace(0,0,n)'          % amarelo -> vermelho
];
colormap(custom_colormap);

% Eixos e titulo
xlabel('Espaco (x)', 'FontSize', 12);
ylabel('Tempo (t)', 'FontSize', 12);
title(['Mapa no Espaco-Tempo | ', ...
       replace(Inputs.RHSIdentifier, '_', ' '), ' | ', ...
       upper(Inputs.TimeMarchingSchemeIdentifier), ' | ', ...
       upper(Inputs.DifferentiationSchemeIdentifier)], 'FontSize', 14);

% Escala de cor centralizada em zero para melhor contraste
caxis([-1 1]*max(abs(u(:))));
