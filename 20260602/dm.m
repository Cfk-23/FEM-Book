% ============================================================
% 算例二：二维两杆桁架
% 节点: 1(1,0), 2(0,0), 3(1,1)
% 单元1:1-3, 单元2:2-3, E=1, A=1
% 边界条件: 节点1和节点2固定
% 载荷: 节点3 Fx=10, Fy=0
% ============================================================
clear; clc;

%% --------------------- 1. 前处理 ---------------------------
nsd = 2;            % 空间维数
ndof = 2;           % 每节点自由度数
nnp = 3;            % 节点数
nel = 2;            % 单元数
nen = 2;            % 每单元节点数

% 节点坐标 (x, y)
x = [0; 1; 2];
y = [0; 0; 0];

% 单元连接 [节点1, 节点2]
IEN = [1, 2;
       2, 3];

% 材料与截面
E = [100; 200];         % 弹性模量
A = [1; 1];         % 截面积

% 位移边界条件: 节点1和节点2完全固定
fixed_dof = [1, 2, 4, 6];   % u1, v1, u2, v2
fixed_val = [0, 0, 0, 0];

% 载荷: 节点3 水平力 Fx=10
force_dof = [5, 6];         % u3, v3
force_val = [10, 0];

n_global_dof = nnp * ndof;

fprintf('========== 算例二：二维两杆桁架 ==========\n');
fprintf('节点数: %d, 单元数: %d, 总自由度数: %d\n', nnp, nel, n_global_dof);
fprintf('固定自由度: %s\n', mat2str(fixed_dof));
fprintf('载荷自由度: %s\n', mat2str(force_dof));
fprintf('==========================================\n\n');

%% --------------------- 2. 生成对号矩阵 LM ------------------
LM = zeros(ndof*nen, nel);
for e = 1:nel
    nodes = IEN(e, :);
    for a = 1:nen
        for i = 1:ndof
            row = (a-1)*ndof + i;
            global_dof = (nodes(a)-1)*ndof + i;
            LM(row, e) = global_dof;
        end
    end
end
fprintf('对号矩阵 LM (每列对应一个单元的全局自由度编号):\n');
disp(LM);

%% --------------------- 3. 组装总体刚度矩阵 K ------------------
K = zeros(n_global_dof, n_global_dof);

for e = 1:nel
    n1 = IEN(e,1); n2 = IEN(e,2);
    x1 = x(n1); y1 = y(n1);
    x2 = x(n2); y2 = y(n2);
    
    dx = x2 - x1;
    dy = y2 - y1;
    L = sqrt(dx^2 + dy^2);
    c = dx / L;
    s = dy / L;
    
    k_local = E(e) * A(e) / L;   % 修改点：使用 E*A/L
    k_e = k_local * [ c^2,  c*s, -c^2, -c*s;
                      c*s,  s^2, -c*s, -s^2;
                     -c^2, -c*s,  c^2,  c*s;
                     -c*s, -s^2,  c*s,  s^2 ];
    
    for a = 1:4
        for b = 1:4
            global_row = LM(a, e);
            global_col = LM(b, e);
            K(global_row, global_col) = K(global_row, global_col) + k_e(a, b);
        end
    end
end

% 输出总体刚度矩阵
fprintf('\n总体刚度矩阵 K (6x6):\n');
disp(K);
fprintf('对称性检查 (最大不对称量): %e\n', max(max(abs(K-K'))));
if rank(K) < n_global_dof
    fprintf('施加边界条件前总体刚度矩阵奇异（秩亏），符合预期。\n\n');
end

%% --------------------- 4. 方程求解 (缩减法) ------------------
all_dof = 1:n_global_dof;
free_dof = setdiff(all_dof, fixed_dof);
n_free = length(free_dof);

K_FF = K(free_dof, free_dof);
K_FE = K(free_dof, fixed_dof);
K_EF = K(fixed_dof, free_dof);
K_EE = K(fixed_dof, fixed_dof);

F = zeros(n_global_dof, 1);
F(force_dof) = force_val;
F_F = F(free_dof);
F_E = F(fixed_dof);

d_E = fixed_val(:);
d_F = K_FF \ (F_F - K_FE * d_E);

d = zeros(n_global_dof, 1);
d(free_dof) = d_F;
d(fixed_dof) = d_E;

r_E = K_EF * d_F + K_EE * d_E - F_E;

fprintf('========== 求解结果 ==========\n');
fprintf('节点位移 (全局自由度顺序):\n');
for i = 1:n_global_dof
    fprintf('d_%d = %12.6f\n', i, d(i));
end
fprintf('\n约束反力 (固定自由度上的力):\n');
for i = 1:length(fixed_dof)
    fprintf('r_%d = %12.6f\n', fixed_dof(i), r_E(i));
end

if rank(K_FF) == n_free
    fprintf('\n缩减后刚度矩阵 K_FF 非奇异，求解成功。\n');
else
    fprintf('\n警告：缩减后刚度矩阵仍然奇异，请检查边界条件！\n');
end

%% --------------------- 5. 后处理 ------------------
fprintf('\n========== 单元结果 ==========\n');
for e = 1:nel
    n1 = IEN(e,1); n2 = IEN(e,2);
    dx = x(n2) - x(n1);
    dy = y(n2) - y(n1);
    L = sqrt(dx^2 + dy^2);
    c = dx / L;
    s = dy / L;
    
    node1_dof = (n1-1)*ndof+1 : (n1-1)*ndof+ndof;
    node2_dof = (n2-1)*ndof+1 : (n2-1)*ndof+ndof;
    de = [d(node1_dof); d(node2_dof)];
    
    strain = ((de(3)-de(1))*c + (de(4)-de(2))*s) / L;
    stress = E(e) * strain;
    axial_force = stress * A(e);
    
    fprintf('单元 %d: 长度=%.4f, c=%.4f, s=%.4f\n', e, L, c, s);
    fprintf('        应变=%.6f, 应力=%.6f, 轴力=%.6f\n', strain, stress, axial_force);
end

fprintf('\n========== 程序运行结束 ==========\n');