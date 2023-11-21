# TransferChan
> 喵喵喵

### > Save_me_plz.mlx
Github 歧视 matlab 是叭（恼

```maltlab
% 解方程 %
clear;clc;

syms x a0 a1 a2 a3 a4 a5

A = [1 1 1 1 1 1 ;
1 2 4 8 16 32;
0 1 2 3 4 5;
0 1 4 12 32 80;
0 0 2 6 12 20;
0 0 2 12 48 160];

b = [0.5; 1; 0.75; 0; 0; 0];

a = [a0; a1; a2; a3; a4; a5];

solutions = solve(A * a == b, a)


% 拉普拉斯 %
clear;clc;

syms t s a

y = t^2
disp('拉普拉斯: ');
Y(s) = laplace(y, t, s)  % 拉普拉斯


syms s t a

Y = (s + 1 + s^(-1))/(s^2 + 2*s + 5)  % 1/(s^2 + 2*s + 5)
disp('反拉普拉斯: ');
y(t) = ilaplace(Y, s, t)  % 反拉普拉斯


% 微积分 %
clear;clc;

A=[1 2; 3 4]

syms t
% y = exp{At}:
d_y = diff(exp(A*t))  % 一阶导
d4_y = diff(exp(A*t), t, 4)  % 四阶导

syms x;

% 使用 solve 求根
y = (x+3)*(x-1)*(x-5)*(x-8);  % 直接写表达式即可，这里只是方便看根
x_root = solve(y, x)

% 使用 roots 求根
% 假设 y=x^3+3x^2+3x
% 将系数 由高到低 填入
p = [1 3 3 0];
roots = roots(p)


% 矩阵性质 %
clear;clc;

A = [0 1 0 0;
     0 0 1 0;
     0 0 0 1;
    -680 -176 -86 -6];

[V, D] = eig(A);
disp('特征值：');
disp(D);
disp('特征向量：');
disp(V);

rank = rank(A);
disp(['矩阵的秩：', num2str(rank)]);

det_A = det(A)  % 行列式
conj_A = conj(A.')  % 共轭转置, 其中A.'为非共轭转置

invert_A = A^(-1)  % 矩阵的逆

[U, S, V] = svd(A)  % 奇异值分解

NullSpace = null(A)  % 零空间，默认右零空间
LeftNullSpace = null(conj(A.'))  % 左零空间

[V, D] = eig(A);  % 求特征值和特征向量

% 解状态空间方程 %
clear;clc;

A = [-2 0;
     0 -3];
B = [1; 1];
C = [-3 4];
D = 0;

% 阶跃信号（单位阶跃函数）：heaviside(t)
% 正弦信号：sin(t)
% 方波信号：square(t)
% 脉冲信号：dirac(t)
% 指数信号：exp(t)

syms s t

% u(t) = heaviside(t);
u(t) = 2 * exp(t);
% u(t) = 0;  % zero input response

% x0 = [2/3; 1/2];
x0 = [0; 0];  % zero initial state response

U_s = laplace(u(t), t, s);  % 对输入信号u(t)进行拉普拉斯变换

X_s = inv(s * eye(size(A)) - A) * x0 + inv(s * eye(size(A)) - A) * B * U_s  % 计算X(s)的频域表达式
Y_s = (C * inv(s*eye(size(A)) - A) * B + D) * U_s + C * inv(s*eye(size(A)) - A) * x0  % 计算Y(s)的频域表达式

X_t = ilaplace(X_s, s, t)  % 对X(s)进行反拉普拉斯变换得到时域解析解
Y_t = ilaplace(Y_s, s, t)  % 对Y(s)进行反拉普拉斯变换得到时域解析解


% 状态方程 <-> 传递函数 %
clear;clc;

A = [0 1 0 0;
     0 0 1 0;
     0 0 0 1;
    -680 -176 -86 -6];
B = [0; 0; 0; 1];
C = [100 20 10 0];
D = 0;

sys = ss(A, B, C, D);
H = tf(sys)  % 求传递函数 G(s)

G_num = [0 0 10 20 100];
G_den = [1 6 86 176 680];
[gA, gB, gC, gD] = tf2ss(G_num,G_den)  % 求出状态矩阵 // 状态变量的选取原则不一定和实验者一致。

% Diagonal Canonical Form ABCD 对角标准型 %

[V, D] = eig(A);  % 求特征值和特征向量

A_dcf = V^(-1) * A * V
B_dcf = V^(-1) * B
C_dcf = C * V
D_dcf = D


% 能控性和能控标准型 %
clear;clc;

A = [0, 1;
     4, 3];  % 系统矩阵

B = [0; 1];  % 输入矩阵

C = [100 20];

D = 0;

% 计算能控性矩阵
M_c = ctrb(A, B)

% 检查能控性
if rank(M_c) == size(A, 1)
    disp('系统是能控的');

    % 计算能控标准型 (SI)
    T = M_c;
    % T = fliplr(M_c);
    A_tilde = T \ A * T
    B_tilde = T \ B
    C_tilde = C * T
    D_tilde = D
    
else
    disp('系统不是能控的');
end


% 计算伴随矩阵 %
clear;clc;

A = [4, 3;  
     6, 3];

adjoint_A = adjointMatrix(A)


function adjA = adjointMatrix(A)
    % 计算矩阵A的伴随矩阵
    n = size(A, 1);
    adjA = zeros(n); % 初始化伴随矩阵

    for i = 1:n
        for j = 1:n
            %  创建余子矩阵
            subA = A([1:i-1, i+1:end], [1:j-1, j+1:end]);
            %  计算代数余子式
            adjA(j, i) = (-1)^(i+j) * det(subA);
        end
    end
end

% end %

