beta=[1 1.48];
alpha1=[1 -0.57];
alpha0=[1 -1];
alpha=conv(alpha1, alpha0);
T = 0.1;

W = tf(beta, alpha, T);
beta1=[1];
beta0=[1 1.48];
r=2; %порядок астатизма
k=0; %кратность корня z=1
B = transpose(roots(beta));
AA = transpose(roots(alpha));

for i=1:length(AA)
    if (A(i)==1)
        k = k+1;
    end
end

for i=1:length(B)
    if (B(i)==1)
        k = k+1;
    end
end
    
display(k);

l=r-k;
o =((length(alpha0)-1) + l - 1 - (length(beta1)-1) + (length(alpha)-1));
nv = 3;
nm = 1;
nn = 1;
%их я ручками посчитала, потому что не поняла, как тут неравенства решать
% z^3 = (m1*z + m0)*(z+1.48) + (n1*z + n0)*(z-1)*(z-1) - решаем данное уравнение

A = [[1 0 0 0]; [-2 1 1 0]; [1 -2 1.48 1]; [0 1 0 1.48]];
B = [1 0 0 0];
X = inv(A) * transpose(B);  
n1 = X(1);
n0 = X(2);
m1 = X(3);
m0 = X(4);
M = [m1 m0];
N = [n1 n0];
R = tf(conv(alpha1, M), conv(conv(beta1, N), [1 -1]), T);

t = 0 : 0.1 : 1;
WW = feedback(series(R, W),[1], -1);

g = [1 2 1];
[n, d] = tfdata(feedback([1],series(R, W), -1));
n = n{1};
d = d{1};

c0 =(polyval(n, 1) / polyval(d, 1));

c1 = T * polyval(deconv( n, [1 -1]),1) /polyval(d, 1);

c2 = T * T * polyval(deconv(n, [1 -2 1]), 1) / polyval(d, 1);

%g'' = 2 и с0 и с1 равны 0 => ошибка равна 2*с2 

%извините, я не хочу выносить это в отдельную m-функцию, потому что не вижу в этом никакого смысла
disp(2*c2);

subplot(2,1, 1), step(WW, t);
subplot(2, 1, 2), lsim(feedback([1],series(R, W), -1), t.^2 + 2*t +1, t);

% e(inf) = 0.0391
