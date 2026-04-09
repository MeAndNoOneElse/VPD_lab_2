clear; clc; close all;

%% ========================
% 1. ПАРАМЕТРЫ И ДАННЫЕ U-I
% ========================
csv_file = '../U_I.csv';  % путь к вашему CSV-файлу с U и I

% Геометрические параметры ротора (измеренные)
m = 15.7e-3;        % кг
r_rotor = 11.5e-3;  % м
gear_ratio = 48;    % передаточное отношение редуктора EV3
L = 0.0047;         % Гн (индуктивность, дана в методичке)
Umax = 9.0;         % В (максимальное напряжение питания)

% Чтение вольтамперной характеристики
T = readtable(csv_file, 'Delimiter', ';');
U = str2double(strrep(string(T{:,2}), ',', '.'));
I = str2double(strrep(string(T{:,3}), ',', '.'));

% Определение сопротивления методом наименьших квадратов (формула 41)
R = sum(U .* I) / sum(I.^2);

% Расчёт момента инерции
J_ed = (m * r_rotor^2) / 2;      % момент инерции якоря
J = gear_ratio^2 * J_ed;          % приведённый момент инерции

% Сохраняем в workspace для Simulink
assignin('base', 'J', J);
assignin('base', 'L', L);
assignin('base', 'R', R);

fprintf('===== ПАРАМЕТРЫ =====\n');
fprintf('R = %.4f Ohm\n', R);
fprintf('J_ed = %.6f kg*m^2\n', J_ed);
fprintf('J = %.6f kg*m^2\n', J);
fprintf('L = %.4f H\n', L);

%% ========================
% 2. ОБРАБОТКА ЭКСПЕРИМЕНТАЛЬНЫХ ДАННЫХ (переходные процессы)
% ========================
data_folder = '../data';
graphs_folder = 'graphs';
if ~exist(graphs_folder, 'dir'), mkdir(graphs_folder); end

% Создаём подпапки для графиков
subfolders = {'angles', 'velocities', 'comparison', 'parameters', 'static'};
for i = 1:length(subfolders)
    folder_path = fullfile(graphs_folder, subfolders{i});
    if ~exist(folder_path, 'dir'), mkdir(folder_path); end
end

% Уровни управляющего сигнала (ШИМ, %)
voltages = [-50, -45, -40, -35, -30, -25, -20, -15, -10, ...
            10, 15, 20, 25, 30, 35, 40, 45, 50];
n_files = length(voltages);

% Массивы для хранения параметров
k_all = NaN(n_files, 1);
Tm_all = NaN(n_files, 1);
omega_steady_deg_all = NaN(n_files, 1);   % установившаяся скорость (град/с)
U_ctrl_all = NaN(n_files, 1);              % напряжение для каждого ШИМ

%% ========================
% 3. ГРАФИКИ: ВСЕ УГЛЫ И ВСЕ СКОРОСТИ НА ОДНОМ РИСУНКЕ
% ========================
figure(1); clf; hold on;
colors = jet(n_files);
figure(2); clf; hold on;

leg_angle = cell(n_files, 1);
leg_vel = cell(n_files, 1);

% Сопоставление ШИМ и напряжения (из таблицы ВАХ)
% Индексы: voltages(1) = -50% -> U(1), voltages(2) = -45% -> U(2), ...
% Для положительных напряжений — соответствующие положительные значения из таблицы
U_from_table = [-2.75; -2.53; -2.27; -1.98; -1.68; -1.39; -1.19; -0.94; -0.60;
                 0.61;  0.95;  1.19;  1.36;  1.68;  1.92;  1.99;  2.53;  2.75];

for i = 1:n_files
    U_pr = voltages(i);
    U_ctrl_all(i) = U_from_table(i);
    
    filename = fullfile(data_folder, sprintf('data%d.txt', U_pr));
    if ~isfile(filename)
        filename = fullfile(data_folder, sprintf('data%d', U_pr));
        if ~isfile(filename)
            warning('Файл %d не найден', U_pr);
            continue;
        end
    end

    data = readmatrix(filename);
    if isempty(data), continue; end

    time = data(:, 1);
    angle_deg = data(:, 2);
    omega_deg = data(:, 3);
    angle_rad = angle_deg * pi / 180;

    % --- Определение установившейся скорости (последние 10% времени) ---
    t_end = max(time);
    idx_steady = time > 0.9 * t_end;
    omega_steady_deg_all(i) = mean(omega_deg(idx_steady));
    
    % --- Графики всех углов и скоростей ---
    figure(1);
    plot(time, angle_deg, 'Color', colors(i,:), 'LineWidth', 1.2);
    leg_angle{i} = sprintf('U = %d%%', U_pr);

    figure(2);
    plot(time, omega_deg, 'Color', colors(i,:), 'LineWidth', 1.2);
    leg_vel{i} = sprintf('U = %d%%', U_pr);

    % --- Аппроксимация переходного процесса угла ---
    % Модель: θ(t) = k * U * (t - Tm * (1 - exp(-t/Tm)))
    fun = @(par, t) U_pr * par(1) * (t - par(2) * (1 - exp(-t/par(2))));
    opts = optimoptions('lsqcurvefit', 'Display', 'off');

    try
        % Начальные приближения: k ~ 0.18, Tm ~ 0.08
        p = lsqcurvefit(fun, [0.18, 0.08], time, angle_rad, [], [], opts);
        k_all(i) = p(1);
        Tm_all(i) = p(2);
    catch
        warning('Аппроксимация для %d%% не удалась', U_pr);
    end
end

% --- Оформление графика углов ---
figure(1);
xlabel('Time, s'); ylabel('Angle, deg');
title('Transition processes of rotation angle');
grid on;
legend(leg_angle, 'Location', 'eastoutside', 'FontSize', 7);
saveas(gcf, fullfile(graphs_folder, 'all_angles.png'));

% --- Оформление графика скоростей ---
figure(2);
xlabel('Time, s'); ylabel('Angular velocity, deg/s');
title('Transition processes of angular velocity');
grid on;
legend(leg_vel, 'Location', 'eastoutside', 'FontSize', 7);
saveas(gcf, fullfile(graphs_folder, 'all_velocities.png'));

%% ========================
% 4. ИНДИВИДУАЛЬНЫЕ ГРАФИКИ УГЛА И СКОРОСТИ (с аппроксимацией)
% ========================
for i = 1:n_files
    if isnan(k_all(i)), continue; end

    U_pr = voltages(i);
    filename = fullfile(data_folder, sprintf('data%d.txt', U_pr));
    if ~isfile(filename), filename = fullfile(data_folder, sprintf('data%d', U_pr)); end

    data = readmatrix(filename);
    time = data(:, 1);
    angle_deg = data(:, 2);
    omega_deg = data(:, 3);
    t_end = max(time);
    
    % Аппроксимированные кривые
    t_apr = linspace(0, t_end, 200);
    theta_apr_deg = U_pr * k_all(i) * (t_apr - Tm_all(i) * (1 - exp(-t_apr/Tm_all(i)))) * 180 / pi;
    omega_apr_deg = k_all(i) * U_pr * (1 - exp(-t_apr/Tm_all(i))) * 180 / pi;

    % График угла
    figure(10 + i);
    plot(time, angle_deg, 'b.-', 'LineWidth', 1.2, 'MarkerSize', 6); hold on;
    plot(t_apr, theta_apr_deg, 'r-', 'LineWidth', 1.5);
    title(sprintf('Angle U = %d%%, k = %.4f, Tm = %.4f s', U_pr, k_all(i), Tm_all(i)));
    xlabel('Time, s'); ylabel('Angle, deg');
    legend('Experiment', 'Approximation', 'Location', 'best');
    grid on;
    saveas(gcf, fullfile(graphs_folder, 'angles', sprintf('angle_%d.png', U_pr)));
    close(gcf);

    % График скорости
    figure(30 + i);
    plot(time, omega_deg, 'b.-', 'LineWidth', 1.2, 'MarkerSize', 6); hold on;
    plot(t_apr, omega_apr_deg, 'r-', 'LineWidth', 1.5);
    title(sprintf('Speed U = %d%%, k = %.4f, Tm = %.4f s', U_pr, k_all(i), Tm_all(i)));
    xlabel('Time, s'); ylabel('Angular velocity, deg/s');
    legend('Experiment', 'Approximation', 'Location', 'best');
    grid on;
    saveas(gcf, fullfile(graphs_folder, 'velocities', sprintf('velocity_%d.png', U_pr)));
    close(gcf);
end

%% ========================
% 5. ВОЛЬТ-АМПЕРНАЯ ХАРАКТЕРИСТИКА
% ========================
I_line = linspace(min(I), max(I), 200);
U_line = R * I_line;

figure('Name','U(I)');
plot(I, U, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 12); hold on;
plot(I_line, U_line, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Current I, A'); ylabel('Voltage U, V');
title('Volt-ampere characteristic');
legend('Experimental', sprintf('U = %.4f \\cdot I', R), 'Location', 'best');
saveas(gcf, fullfile(graphs_folder, 'U_I_dependence.png'));
close(gcf);

%% ========================
% 6. СТАТИЧЕСКАЯ ХАРАКТЕРИСТИКА U(ω_уст) — ОПРЕДЕЛЕНИЕ k_e (п. 2.6 ЗАДАНИЯ)
% ========================
% Убираем NaN
valid_steady = ~isnan(omega_steady_deg_all) & ~isnan(U_ctrl_all);
omega_steady_deg_valid = omega_steady_deg_all(valid_steady);
U_ctrl_valid = U_ctrl_all(valid_steady);

% Переводим скорость в рад/с
omega_steady_rad = omega_steady_deg_valid * pi / 180;

% Аппроксимация прямой через 0: U = k_e * ω
% Метод наименьших квадратов для модели y = a*x
k_e_static = sum(U_ctrl_valid .* omega_steady_rad) / sum(omega_steady_rad.^2);

fprintf('\n===== ОПРЕДЕЛЕНИЕ k_e ИЗ СТАТИЧЕСКОЙ ХАРАКТЕРИСТИКИ =====\n');
fprintf('k_e = %.4f V·s/rad\n', k_e_static);

% График U(ω_уст)
figure('Name','U(omega_steady)');
plot(omega_steady_rad, U_ctrl_valid, 'bo', 'MarkerSize', 10, 'LineWidth', 1.5); hold on;
omega_fit = linspace(0, max(omega_steady_rad), 100);
U_fit = k_e_static * omega_fit;
plot(omega_fit, U_fit, 'r-', 'LineWidth', 1.5);
xlabel('Angular velocity \omega_{steady}, rad/s');
ylabel('Voltage U, V');
title(sprintf('Static characteristic U(\\omega_{steady}), k_e = %.4f V·s/rad', k_e_static));
legend('Experiment', sprintf('U = %.4f \\cdot \\omega', k_e_static), 'Location', 'best');
grid on;
saveas(gcf, fullfile(graphs_folder, 'static', 'U_vs_omega_steady.png'));
close(gcf);

%% ========================
% 7. ЗАВИСИМОСТИ ПАРАМЕТРОВ ОТ НАПРЯЖЕНИЯ
% ========================
valid = ~isnan(k_all) & ~isnan(Tm_all);
U_valid = U_ctrl_all(valid);
k_valid = k_all(valid);
Tm_valid = Tm_all(valid);
omega_steady_rad_all = omega_steady_deg_all * pi / 180;
omega_steady_valid = omega_steady_rad_all(valid);

% Зависимость k(U)
figure('Name','k(U)');
plot(U_valid, k_valid, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Voltage U, V'); ylabel('Coefficient k, rad/(s·%)');
title('Dependence of transfer coefficient k on voltage');
grid on;
saveas(gcf, fullfile(graphs_folder, 'parameters', 'k.png'));
close(gcf);

% Зависимость Tm(U)
figure('Name','Tm(U)');
plot(U_valid, Tm_valid, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Voltage U, V'); ylabel('Time constant T_m, s');
title('Dependence of time constant T_m on voltage');
grid on;
saveas(gcf, fullfile(graphs_folder, 'parameters', 'Tm.png'));
close(gcf);

% Зависимость ω_уст(U)
figure('Name','omega_steady(U)');
plot(U_valid, omega_steady_valid * 180/pi, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Voltage U, V'); ylabel('Steady-state angular velocity, deg/s');
title('Dependence of steady-state speed on voltage');
grid on;
saveas(gcf, fullfile(graphs_folder, 'parameters', 'omega_steady.png'));
close(gcf);

% Средние значения
k_mean = mean(k_valid);
Tm_mean = mean(Tm_valid);

fprintf('\n===== СРЕДНИЕ ПАРАМЕТРЫ =====\n');
fprintf('k_mean  = %.4f rad/(s·%%)\n', k_mean);
fprintf('Tm_mean = %.4f s\n', Tm_mean);

%% ========================
% 8. ПЕРЕХОДНЫЙ ПРОЦЕСС ТОКА (ODE45 для полной модели)
% ========================
% Используем k_e, полученный из статической характеристики
ke = k_e_static;
km = ke;  % для ДПТ с постоянными магнитами km = ke

assignin('base', 'ke', ke);
assignin('base', 'km', km);

% Выбираем одно напряжение для демонстрации (например, +50% -> 2.75 В)
U_demo = 2.75;
tspan = [0 1];

% Система дифференциальных уравнений полной модели (вход-состояние)
% dx(1) = omega, dx(2) = I
f = @(t, x) [ (km/J) * x(2);
              (1/L)*U_demo - (ke/L)*x(1) - (R/L)*x(2) ];

[t_ode, x_ode] = ode45(f, tspan, [0; 0]);

figure('Name','Current_transient');
plot(t_ode, x_ode(:,2), 'b-', 'LineWidth', 1.5);
grid on; 
xlabel('t, s'); 
ylabel('I, A');
title(sprintf('Current transient I(t) for U = %.2f V (full model)', U_demo));
saveas(gcf, fullfile(graphs_folder, 'ode_current.png'));
close(gcf);

%% ========================
% 9. SIMULINK-МОДЕЛЬ (ПОЛНАЯ МОДЕЛЬ «ВХОД-СОСТОЯНИЕ»)
% ========================
if license('test', 'Simulink')
    modelName = 'EV3_motor_full_model';
    
    % Закрываем модель, если открыта
    if bdIsLoaded(modelName), close_system(modelName, 0); end
    if isfile([modelName '.slx']), delete([modelName '.slx']); end
    
    % Создаём полную модель
    create_full_motor_model(modelName);
    
    % Для нескольких напряжений (выборочно)
    test_voltages = [-2.75, -1.68, 1.68, 2.75];
    
    for i = 1:length(test_voltages)
        U_test = test_voltages(i);
        
        % Явно передаём переменные в workspace перед симуляцией
        assignin('base', 'U_sim', U_test);
        assignin('base', 'ke', k_e_static);
        assignin('base', 'km', k_e_static);
        assignin('base', 'R', R);
        assignin('base', 'L', L);
        assignin('base', 'J', J);
        
        % Находим соответствующий ШИМ для подписи
        idx = find(abs(U_ctrl_all - U_test) < 0.01, 1);
        if ~isempty(idx)
            pwm_label = voltages(idx);
        else
            pwm_label = round(U_test / 9 * 100);
        end
        
        % Запуск симуляции
        t_end_sim = 1.0;
        set_param(modelName, 'StopTime', num2str(t_end_sim));
        
        try
            % Очищаем старые переменные
            clear omega theta current;
            
            % Запускаем симуляцию
            simOut = sim(modelName, 'SaveOutput', 'on', 'OutputSaveName', 'yout');
            
            % Извлекаем данные
            if exist('omega', 'var')
                t_sim = omega.Time;
                omega_sim = omega.Data;
                theta_sim = theta.Data;
                current_sim = current.Data;
            else
                % Альтернативный способ: из simOut
                t_sim = simOut.tout;
                % Нужно понять, где какие сигналы — зависит от порядка блоков
                if size(simOut.yout, 2) >= 3
                    omega_sim = simOut.yout(:,1);
                    theta_sim = simOut.yout(:,2);
                    current_sim = simOut.yout(:,3);
                else
                    warning('Не удалось извлечь данные для U = %.2f', U_test);
                    continue;
                end
            end
            
            % Загружаем экспериментальные данные для сравнения
            pwm_idx = find(voltages == pwm_label, 1);
            if ~isempty(pwm_idx) && ~isnan(k_all(pwm_idx))
                filename = fullfile(data_folder, sprintf('data%d.txt', voltages(pwm_idx)));
                if isfile(filename)
                    exp_data = readmatrix(filename);
                    time_exp = exp_data(:, 1);
                    angle_exp = exp_data(:, 2);
                    omega_exp = exp_data(:, 3);
                    
                    % Аппроксимация для этого уровня
                    k_curr = k_all(pwm_idx);
                    Tm_curr = Tm_all(pwm_idx);
                    t_apr = linspace(0, min(t_end_sim, max(time_exp)), 200);
                    theta_apr_deg = voltages(pwm_idx) * k_curr * (t_apr - Tm_curr * (1 - exp(-t_apr/Tm_curr))) * 180 / pi;
                    omega_apr_deg = k_curr * voltages(pwm_idx) * (1 - exp(-t_apr/Tm_curr)) * 180 / pi;
                    
                    % График сравнения
                    figure(100 + i);
                    
                    subplot(2,1,1);
                    plot(time_exp, angle_exp, 'b-', 'LineWidth', 1.2); hold on;
                    plot(t_sim, theta_sim * 180/pi, 'r-', 'LineWidth', 1.5);
                    plot(t_apr, theta_apr_deg, 'g--', 'LineWidth', 1.5);
                    xlabel('Time, s'); ylabel('Angle, deg');
                    title(sprintf('Angle comparison, U = %.2f V (%d%%)', U_test, pwm_label));
                    legend('Experiment', 'Simulink (full model)', 'Approximation (1st order)', 'Location', 'best');
                    grid on;
                    
                    subplot(2,1,2);
                    plot(time_exp, omega_exp, 'b-', 'LineWidth', 1.2); hold on;
                    plot(t_sim, omega_sim * 180/pi, 'r-', 'LineWidth', 1.5);
                    plot(t_apr, omega_apr_deg, 'g--', 'LineWidth', 1.5);
                    xlabel('Time, s'); ylabel('Angular velocity, deg/s');
                    title(sprintf('Speed comparison, U = %.2f V (%d%%)', U_test, pwm_label));
                    legend('Experiment', 'Simulink (full model)', 'Approximation (1st order)', 'Location', 'best');
                    grid on;
                    
                    saveas(gcf, fullfile(graphs_folder, 'comparison', sprintf('comparison_%d.png', pwm_label)));
                    close(gcf);
                end
            end
        catch ME
            warning('Ошибка при симуляции для U = %.2f: %s', U_test, ME.message);
        end
    end
    
    % Закрываем модель
    if bdIsLoaded(modelName), close_system(modelName, 0); end
else
    fprintf('Simulink не найден, верификация пропущена.\n');
end

%% ========================
% 10. ВЫВОД ИТОГОВЫХ ЗНАЧЕНИЙ ДЛЯ ОТЧЁТА
% ========================
fprintf('\n===== ИТОГОВЫЕ ПАРАМЕТРЫ ДЛЯ ОТЧЁТА =====\n');
fprintf('R   = %.4f Ohm\n', R);
fprintf('J   = %.6f kg*m^2\n', J);
fprintf('L   = %.4f H (из методички)\n', L);
fprintf('k   = %.4f rad/(s·%%)\n', k_mean);
fprintf('Tm  = %.4f s\n', Tm_mean);
fprintf('ke  = %.4f V·s/rad (из статической характеристики U(ω_уст))\n', k_e_static);
fprintf('km  = %.4f N·m/A (принято равным ke)\n', k_e_static);

disp(' ');
disp('Обработка завершена. Все графики сохранены в папку graphs/');

%% ========================
% 11. ФУНКЦИЯ СОЗДАНИЯ ПОЛНОЙ МОДЕЛИ «ВХОД-СОСТОЯНИЕ»
% ========================
function create_full_motor_model(modelName)
    % Создаёт Simulink-модель, соответствующую схеме на рис. 6 методички
    % Полная модель: dω/dt = (km/J)*I, dI/dt = (1/L)*U - (ke/L)*ω - (R/L)*I
    
    % Создаём новую систему
    new_system(modelName);
    open_system(modelName);
    set_param(modelName, 'StopTime', '1.0', 'Solver', 'ode45');
    
    % 1. Блок Constant для входного напряжения
    add_block('simulink/Sources/Constant', [modelName '/U_ctrl']);
    set_param([modelName '/U_ctrl'], 'Value', 'U_sim');
    
    % 2. Блок Gain для 1/L
    add_block('simulink/Math Operations/Gain', [modelName '/Gain_1_L']);
    set_param([modelName '/Gain_1_L'], 'Gain', '1/L');
    
    % 3. Блок Gain для ke/L
    add_block('simulink/Math Operations/Gain', [modelName '/Gain_ke_L']);
    set_param([modelName '/Gain_ke_L'], 'Gain', 'ke/L');
    
    % 4. Блок Gain для R/L
    add_block('simulink/Math Operations/Gain', [modelName '/Gain_R_L']);
    set_param([modelName '/Gain_R_L'], 'Gain', 'R/L');
    
    % 5. Блок Gain для km/J
    add_block('simulink/Math Operations/Gain', [modelName '/Gain_km_J']);
    set_param([modelName '/Gain_km_J'], 'Gain', 'km/J');
    
    % 6. Сумматор для тока (3 входа)
    add_block('simulink/Math Operations/Sum', [modelName '/Sum_I']);
    set_param([modelName '/Sum_I'], 'Inputs', '++-', 'IconShape', 'round');
    
    % 7. Сумматор для скорости (2 входа)
    add_block('simulink/Math Operations/Sum', [modelName '/Sum_omega']);
    set_param([modelName '/Sum_omega'], 'Inputs', '+-', 'IconShape', 'round');
    
    % 8. Интегратор тока
    add_block('simulink/Continuous/Integrator', [modelName '/Integ_I']);
    set_param([modelName '/Integ_I'], 'InitialCondition', '0');
    
    % 9. Интегратор скорости
    add_block('simulink/Continuous/Integrator', [modelName '/Integ_omega']);
    set_param([modelName '/Integ_omega'], 'InitialCondition', '0');
    
    % 10. Интегратор угла
    add_block('simulink/Continuous/Integrator', [modelName '/Integ_theta']);
    set_param([modelName '/Integ_theta'], 'InitialCondition', '0');
    
    % 11. To Workspace для скорости
    add_block('simulink/Sinks/To Workspace', [modelName '/ToWorkspace_omega']);
    set_param([modelName '/ToWorkspace_omega'], 'VariableName', 'omega', 'SaveFormat', 'Structure');
    
    % 12. To Workspace для угла
    add_block('simulink/Sinks/To Workspace', [modelName '/ToWorkspace_theta']);
    set_param([modelName '/ToWorkspace_theta'], 'VariableName', 'theta', 'SaveFormat', 'Structure');
    
    % 13. To Workspace для тока
    add_block('simulink/Sinks/To Workspace', [modelName '/ToWorkspace_I']);
    set_param([modelName '/ToWorkspace_I'], 'VariableName', 'current', 'SaveFormat', 'Structure');
    
    % ================ СОЕДИНЕНИЯ ================
    
    % Подключаем U_ctrl -> Gain_1_L
    add_line(modelName, 'U_ctrl/1', 'Gain_1_L/1');
    
    % Подключаем Gain_1_L -> Sum_I (первый вход, +)
    add_line(modelName, 'Gain_1_L/1', 'Sum_I/1');
    
    % Подключаем Sum_I -> Integ_I
    add_line(modelName, 'Sum_I/1', 'Integ_I/1');
    
    % Подключаем Integ_I -> Gain_km_J
    add_line(modelName, 'Integ_I/1', 'Gain_km_J/1');
    
    % Подключаем Gain_km_J -> Sum_omega (первый вход, +)
    add_line(modelName, 'Gain_km_J/1', 'Sum_omega/1');
    
    % Подключаем Sum_omega -> Integ_omega
    add_line(modelName, 'Sum_omega/1', 'Integ_omega/1');
    
    % Обратная связь по току: Integ_I -> Gain_R_L -> Sum_I (третий вход, -)
    add_line(modelName, 'Integ_I/1', 'Gain_R_L/1');
    add_line(modelName, 'Gain_R_L/1', 'Sum_I/3');
    
    % Обратная связь по скорости: Integ_omega -> Gain_ke_L -> Sum_I (второй вход, -)
    add_line(modelName, 'Integ_omega/1', 'Gain_ke_L/1');
    add_line(modelName, 'Gain_ke_L/1', 'Sum_I/2');
    
    % Нулевой сигнал на второй вход Sum_omega (можно оставить неподключенным)
    % или добавить Constant 0. Оставляем как есть — Simulink сам подставит 0.
    
    % Интегратор угла
    add_line(modelName, 'Integ_omega/1', 'Integ_theta/1');
    
    % Подключение To Workspace
    add_line(modelName, 'Integ_omega/1', 'ToWorkspace_omega/1');
    add_line(modelName, 'Integ_theta/1', 'ToWorkspace_theta/1');
    add_line(modelName, 'Integ_I/1', 'ToWorkspace_I/1');
    
    % Сохраняем модель
    save_system(modelName);
    
    fprintf('Simulink-модель "%s" успешно создана.\n', modelName);
end

