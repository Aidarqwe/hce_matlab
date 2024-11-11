clear;
clc;

% Заданные параметры
N = 50; % разрешение сетки
L = 0.014; % размер пластины в метрах (14 мм)
h = L / N; % размер ячейки сетки (м)

x = linspace(0, L, N);
y = linspace(0, L, N);

[X, Y] = meshgrid(x, y);

% Начальные условия
T0 = 300; % начальная температура (К)
T = T0 * ones(size(X)); % инициализация всей пластины

% Параметры временного шага
t = 0; % начальное время (с)
tmax = 60; % конечное время моделирования (с)
Nt = 5000; % количество временных шагов
tau = (tmax - t) / Nt; % временной шаг (с)

% Граничные условия
Tb = 300; % температура на нижней границе (К)
T(1, :) = Tb; % задание температуры на нижней границе

% Свойства материала
rho = 7800; % плотность (кг/м^3)
cp = 460; % удельная теплоемкость (Дж/кг·К)
k = 46; % теплопроводность (Вт/м·К)

% Параметры конвекции и излучения
h_c = 5; % коэффициент теплопередачи (Вт/м^2·К)
T_ref = 300; % температура окружающей среды (К)
epsilon = 0.8; % коэффициент излучательной способности
sigma_sb = 5.67e-8; % постоянная Стефана-Больцмана (Вт/м^2·К^4)

% Параметры лазера
P = 1000; % мощность лазера (Вт)
a = 1;
b = 3;
r = 0.01; % радиус области воздействия лазера (м)
S = pi * r^2; % площадь пятна лазера
x_c = L / 2; % центр по x
y_c = L / 2; % центр по y

% Цикл по времени
alpha = (k * tau) / (rho * cp * h^2);
while t < tmax
    Tprev = T;  % сохранение предыдущего распределения температур
    
    % Метод Гаусса-Зейделя для решения СЛАУ
    for GS = 1:100
        for ix = 2:N-1
            for iy = 2:N-1
                sum = T(ix+1, iy) + T(ix-1, iy) + T(ix, iy+1) + T(ix, iy-1);
                T(ix, iy) = (alpha * sum + Tprev(ix, iy)) / (1 + 4 * alpha);
            end
        end
    end

    % Рассчитываем тепловой поток от конвекции и излучения на верхней границе
    q_conv = -h_c * (Tprev(end, :) - T_ref);
    q_rad = -epsilon * sigma_sb * (Tprev(end, :).^4 - T_ref^4);

    % Воздействие лазера как функция Гаусса на верхней границе
    q_laser = (P / S) * a * exp(-b * ((x - x_c).^2) / r^2); % формула для мощности лазера
    
    % Уравнение для верхней границы с учетом конвекции, излучения и лазера
    T(end, :) = Tprev(end-1, :) + h * (q_laser + q_conv + q_rad) / h_c;

    % Визуализация
    contourf(X, Y, T, 50, 'EdgeColor', 'none');
    xlabel('x, м');
    ylabel('y, м');
    title(sprintf('Температура пластины, К, t = %.2f с', t));
    colorbar;
    drawnow;

    % Обновление времени
    t = t + tau;
end