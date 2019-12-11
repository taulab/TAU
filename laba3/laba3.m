% y'' - 2y' - 2y = u
% u = -psi*y
% psi = 4, y*s > 0; psi = -6, y*s <= 0
%s = y'+3y

% y'' - 2y' - 2y = -psi*y

%если s*y>0 <=>  y*y' + 3y^2 > 0
% Сделаем замену
% x1 = y; x2 = x1'
%x1(x2 + 3x1) > 0
%решение неравенства : ( x1 > 0 && x2 > -3x1) || (x1 < 0 && x2 < -3x1)
%в этой области уравнение принимает вид y'' - 2y' - 2y = -4*y <=> 
%y'' - 2y' + 2y = 0
%После замены  x1 = y; x2 = x1' получается система
%x2' = 2x2 - 2x1
%x1' = x2
% Разделим первое уравнение на второе
% dx2 / dx1 = 2 - 2*x1/x2
% Решение:
% arctg(-x2/x1 + 1) - (1/2)*ln(x2^2/x1^2 - 2*x2/x1 + 2) = ln(x1) + C



%если s*y <= 0 <=>  y*y' + 3y^2 <= 0
% Сделаем замену
% x1 = y; x2 = x1'
%x1(x2 + 3x1) <= 0
%решение неравенства : ( x1 <= 0 && x2 >= -3x1) || (x1 >= 0 && x2 <= -3x1)
%в этой области уравнение принимает вид y'' - 2y' - 2y = 6*y <=> 
%y'' - 2y' - 8y = 0
%После замены  x1 = y; x2 = x1' получается система
%x2' = 2x2 + 8x1
%x1' = x2
% Разделим первое уравнение на второе
% dx2 / dx1 = 2 + 8*x1/x2
% Решение:
%(-2/3)ln(-x2/x1 + 4)- (1/3)ln(x2/x1 + 2) = ln(x) + C


hold on;

axis([-1.5 1.5 -3 3]); 

%запишу функции f и g - решения ду для psi=4 и psi=-6 соответственно с
%учетом вариации константы С

f1 = @(x1, x2) log(abs(x1)) + log(abs(2 - 2.*x2./x1 + x2.^2./(x1.^2))) ./ 2 - atan(1 - x2./x1);
f2 = @(x1, x2) -0.25 + log(abs(x1)) + log(abs(2 - 2.*x2./x1 + x2.^2./(x1.^2))) ./ 2 - atan(1 - x2./x1);
f3 = @(x1, x2) -0.5 + log(abs(x1)) + log(abs(2 - 2.*x2./x1 + x2.^2./(x1.^2))) ./ 2 - atan(1 - x2./x1);
f4 = @(x1, x2) -0.75 + log(abs(x1)) + log(abs(2 - 2.*x2./x1 + x2.^2./(x1.^2))) ./ 2 - atan(1 - x2./x1);
f5 = @(x1, x2) -1 + log(abs(x1)) + log(abs(2 - 2.*x2./x1 + x2.^2./(x1.^2))) ./ 2 - atan(1 - x2./x1);


g1 = @(x1, x2) 2.25 + log(abs(x1)) + (2./3).*log(abs(-(x2./x1) + 4)) + (1/3)*log(abs((x2./x1)+2));
g2 = @(x1, x2) 2 + log(abs(x1)) + (2./3).*log(abs(-(x2./x1) + 4)) + (1/3)*log(abs((x2./x1)+2));
g3 = @(x1, x2) 1.75 + log(abs(x1)) + (2./3).*log(abs(-(x2./x1) + 4)) + (1/3)*log(abs((x2./x1)+2));
g4 = @(x1, x2) 1.5 + log(abs(x1)) + (2./3).*log(abs(-(x2./x1) + 4)) + (1/3)*log(abs((x2./x1)+2));
g5 = @(x1, x2) 1.25 + log(abs(x1)) + (2./3).*log(abs(-(x2./x1) + 4)) + (1/3)*log(abs((x2./x1)+2));


syms x1 x2



% Найдём точки пересечения x2 + 3*x1 = 0 и f1 = 0 для вычисления области
% построения
eqns = [x2 + 3*x1 == 0, f1 + x1*0 == 0];
S = solve(eqns,[x1, x2]);
f1_inter_1 = [-1.5 0 -2 double(abs(S.x2))];
f1_inter_2 = [0 1.5 -double(abs(S.x2)) 2];

% Найдём точки пересечения x2 + 3*x1 = 0 и f2 = 0 для вычисления области
% построения
eqns = [x2 + 3*x1 == 0, f2 + x1*0 == 0];
S = solve(eqns,[x1, x2]);
f2_inter_1 = [-1.5, 0, -2, double(abs(S.x2))];
f2_inter_2 = [0, 1.5, -double(abs(S.x2)), 2];

% Найдём точки пересечения x2 + 3*x1 = 0 и f3 = 0 для вычисления области
% построения
eqns = [x2 + 3*x1 == 0, f3 + x1*0 == 0];
S = solve(eqns,[x1, x2]);
f3_inter_1 = [-1.5, 0, -2, double(abs(S.x2))];
f3_inter_2 = [0, 1.5, -double(abs(S.x2)), 2];

% Найдём точки пересечения x2 + 3*x1 = 0 и f4 = 0 для вычисления области
% построения
eqns = [x2 + 3*x1 == 0, f4 + x1*0 == 0];
S = solve(eqns,[x1, x2]);
f4_inter_1 = [-1.5, 0, -2, double(abs(S.x2))];
f4_inter_2 = [0, 1.5, -double(abs(S.x2)), 2];

% Найдём точки пересечения x2 + 3*x1 = 0 и f5 = 0 для вычисления области
% построения
eqns = [x2 + 3*x1 == 0, f5 + x1*0 == 0];
S = solve(eqns,[x1, x2]);
f5_inter_1 = [-1.5, 0, -2, double(abs(S.x2))];
f5_inter_2 = [0, 1.5, -double(abs(S.x2)), 2];

% Найдём точки пересечения x2 + 3*x1 = 0 и g1 = 0 для вычисления области
% построения
eqns = [x2 + 3*x1 == 0, g1 + x1*0 == 0];
S = solve(eqns,[x1, x2]);
g1_inter_1 = [-double(abs(S.x1(1))), 0, 0 ,2];
g1_inter_2 = [0, double(abs(S.x1(1))), -2, 0];

% Найдём точки пересечения x2 + 3*x1 = 0 и g2 = 0 для вычисления области
% построения
eqns = [x2 + 3*x1 == 0, g2 + x1*0 == 0];
S = solve(eqns,[x1, x2]);
g2_inter_1 = [-double(abs(S.x1(1))), 0, 0 , 2];
g2_inter_2 = [0, double(abs(S.x1(1))), -2, 0];

% Найдём точки пересечения x2 + 3*x1 = 0 и g3 = 0 для вычисления области
% построения
eqns = [x2 + 3*x1 == 0, g3 + x1*0 == 0];
S = solve(eqns,[x1, x2]);
g3_inter_1 = [-double(abs(S.x1(1))), 0, 0, 2];
g3_inter_2 = [0, double(abs(S.x1(1))), -2, 0];

% Найдём точки пересечения x2 + 3*x1 = 0 и g4 = 0 для вычисления области
% построения
eqns = [x2 + 3*x1 == 0, g4 + x1*0 == 0];
S = solve(eqns,[x1, x2]);
g4_inter_1 = [-double(abs(S.x1(1))), 0,0, 2];
g4_inter_2 = [0, double(abs(S.x1(1))), -2, 0];

% Найдём точки пересечения x2 + 3*x1 = 0 и g5 = 0 для вычисления области
% построения
eqns = [x2 + 3*x1 == 0, g5 + x1*0 == 0];
S = solve(eqns,[x1, x2]);
g5_inter_1 = [-double(abs(S.x1(1))), 0, 0, 2];
g5_inter_2 = [0, double(abs(S.x1(1))), -2, 0];

hold on;
% Построим и подпишем оси координат
plot([0 0] , [-3 3], 'k');
plot([-1.5 1.5], [0 0], 'k'); 
txt1 = {'x1'};
text(1.4, -0.1, txt1)
txt2 = {'x2'};
text(0.05, 2.9, txt2)

% Выведем сепаратрису; второй сепаратрисой является x1=0 (см. рис.3)
plot([-1 1], [3 -3], 'k'); 
% Изобразим фазовые траектории системы

fimplicit(f1, f1_inter_1);
fimplicit(f1, f1_inter_2);
fimplicit(f2, f2_inter_1);
fimplicit(f2, f2_inter_2);
fimplicit(f3, f3_inter_1);
fimplicit(f3, f3_inter_2);
fimplicit(f4, f4_inter_1);
fimplicit(f4, f4_inter_2);
fimplicit(f5, f5_inter_1);
fimplicit(f5, f5_inter_2);

fimplicit(g1, g1_inter_1);
fimplicit(g1, g1_inter_2);
fimplicit(g2, g2_inter_1);
fimplicit(g2, g2_inter_2);
fimplicit(g3, g3_inter_1);
fimplicit(g3, g3_inter_2);
fimplicit(g4, g4_inter_1);
fimplicit(g4, g4_inter_2);
fimplicit(g5, g5_inter_1);
fimplicit(g5, g5_inter_2);

% Укажем направления траекторий; x1'=x2, x2'=2*x2 - 4*x1 - psi*x1; psi
% определяем с помощью функции my_f
[x1,x2] = meshgrid(-1.5:0.2:1.5,-3:0.2:3);
u = x2;
v = 2*x2 - 4*x1 - my_f(x1, x2).*x1;
quiver(x1,x2,u,v)






% Из фазового портрета системы видно, что она неустойчива по Ляпунову
disp(0);

% Функция для определения psi
function answ=my_f(x, y)
answ = zeros(size(x, 1), size(x, 2));
for i=1:size(x, 1)
    for j=1:size(x, 2)
        if ((x(i, j)>0) && (y(i, j)>-3*x(i, j)) || (x(i, j)<0) && (y(i, j)<-3*x(i, j)))
            answ(i, j) = 4;
        else
            answ(i, j) = -6;
        end  
    end
end
end
