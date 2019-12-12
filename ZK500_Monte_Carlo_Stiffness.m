clear;clc
NUM = 5000;
ratio_num = 50;
joint_angle = joint_angle_generator(NUM);
% endpoint = zeros(NUM,3);
ratio = linspace(0.8,1.3,ratio_num);
stiffness = zeros(1,ratio_num);
index = zeros(1, ratio_num);
parfor r = 1:1:ratio_num
%     temp = 0;
%     for i = 1:1:NUM
%         temp = temp + Jacobian_calc(joint_angle(i, :), ratio(r));
%     end
    s_idx = avg_calc(joint_angle, ratio(r));
    stiffness(r) = s_idx(1);
    index(r) = s_idx(2);
end
% scatter3(endpoint(:, 1), endpoint(:, 2), endpoint(:, 3));
figure(1)
plot(ratio, stiffness)
figure(2)
plot(ratio, index)

% th = [0,0,0,0,0];
% Jacobian_calc(th,1.0)

function avg = avg_calc(joint_angle,r)
    temp1 = 0;
    idx = 0;
    num = size(joint_angle);
    for i = 1:1:num(1)
        % T_matrix = T_matrix_calc(joint_angle(i, :));
        % endpoint(i, :) = T_matrix(1:3, 4);
        temp = Jacobian_calc(joint_angle(i, :), r);
        temp1 = temp1 + temp(1);
        idx = idx + temp(2);
    end
    avg = temp1/num(1);
    avg = [avg, idx/num(1)];
end

%% func for jth generating
    % args: NUM: number of end points, seed: the random seed
    % return: a matrix with shape of NUMx5 
function joint_angle = joint_angle_generator(NUM)
%     range = [pi*2, (145/180*pi), (230/180*pi), pi, (240/180*pi)];
%     bias = [-pi, -(90/180*pi), -(50/180*pi), -pi/2, -deg2rad(120)];
    
    range = [pi*2, (145/180*pi), (230/180*pi), pi*2, (120/180*pi)];
    bias = [-pi, -(90/180*pi), -(70/180*pi), -pi, 0];
    
%     rng(72, 'twister');
    joint_angle = rand(NUM,5);
    % using for testing the range and bias calc
    % joint_angle = ones(NUM,5);
    % range + bias
    
    range = diag(range);
    joint_angle = joint_angle * range;
    joint_angle = joint_angle +  bias;
end
%%

%% func for end point calc
    % args: th: joint angle
function T_matrix = T_matrix_calc(th, ratio)
    th1=th(1);th2=th(2);th3=th(3);th4=th(4);th5=th(5);

    l1=1.0;l2=0.4;l6=0.5;e=0.37;
    l3 = 1.3*ratio;
    %l3 = 1.3;
    l4 = 1.3/2;
    l5 = l4;

    T1 = [cos(th1) -sin(th1) 0 0;
        sin(th1) cos(th1) 0 0;
        0 0 1 l1;0 0 0 1];
    T2 = [-sin(th2) -cos(th2) 0 l2;
        0 0 -1 0;
        cos(th2) -sin(th2) 0 0;
        0 0 0 1];
    T3 = [cos(th3) -sin(th3) 0 l3;
        sin(th3) cos(th3) 0 0;
        0 0 1 0;0 0 0 1];
    T4 = [cos(th4) -sin(th4) 0 e;
        0 0 -1 -l4-l5;
        sin(th4) cos(th4) 0 0;
        0 0 0 1;];
    T5 = [-cos(th5) sin(th5) 0 0;
        0 0 1 0;
        sin(th5) cos(th5) 0 0;
        0 0 0 1;];
    Tt = [-1 0 0 0;
        0 0 1 l6;
        0 1 0 0;
        0 0 0 1;];
    T_matrix = T1*T2*T3*T4*T5*Tt;
end
%%

%% func for orientation stiffness calc
% this function is only calc the stiffness about force, not including
% torque
    % args: E: eig values(3d array)
    %       R: eig vector(3x3 matrix)
    %       Vector: orientation of external force, a 3x1 array
function ori_stiffness = orientation_stiffness(E, R, Vector)
    n_vec = R\Vector;
%     E = sqrt(E);
    R_matrix = R_calc(n_vec');
    theta = 0.6435;
    temp = 0;
    for i=1:1:201
        t = 2*pi/200*(i-1);
        m_vec = [cos(t)*sin(theta), sin(t)*sin(theta), cos(theta)];
        m_vec = R_matrix * m_vec';
        
        apha0 = (m_vec(1)/E(1))^2+(m_vec(2)/E(2))^2+(m_vec(3)/E(3))^2;
        k0 = m_vec(3)^2+m_vec(2)^2+m_vec(1)^2;
        temp = temp + sqrt(k0/apha0);
        % ori_stiffness = sqrt(ori_stiffness);
    end
    ori_stiffness = temp/201;
end


% function ori_stiffness = orientation_stiffness(E, R, Vector)
%     n_vec = R\Vector;
%     E = sqrt(E);
%     apha0 = (n_vec(1)/E(1))^2+(n_vec(2)/E(2))^2+(n_vec(3)/E(3))^2;
%     k0 = n_vec(3)^2+n_vec(2)^2+n_vec(1)^2;
%     ori_stiffness = sqrt(k0/apha0);
%     % ori_stiffness = sqrt(ori_stiffness);
% end
%%

%%
function deflect = Jacobian_calc(th, ratio)
    th1=th(1);th2=th(2);th3=th(3);th4=th(4);th5=th(5);

    l1=1.05;l2=0.5;l6=0.39;e=0.15;
    l3 = 1.2*ratio;
    l4 = 1.2/2;
    l5 = l4;

    T1 = [cos(th1) -sin(th1) 0 0;
        sin(th1) cos(th1) 0 0;
        0 0 1 l1;0 0 0 1];
    T2 = [-sin(th2) -cos(th2) 0 l2;
        0 0 -1 0;
        cos(th2) -sin(th2) 0 0;
        0 0 0 1];
    T3 = [cos(th3) -sin(th3) 0 l3;
        sin(th3) cos(th3) 0 0;
        0 0 1 0;0 0 0 1];
    T4 = [cos(th4) -sin(th4) 0 e;
        0 0 -1 -l4-l5;
        sin(th4) cos(th4) 0 0;
        0 0 0 1;];
    T5 = [-cos(th5) sin(th5) 0 0;
        0 0 1 0;
        sin(th5) cos(th5) 0 0;
        0 0 0 1;] ;
    Tt = [-1 0 0 0;
        0 0 1 l6;
        0 1 0 0;
        0 0 0 1;];

    T = T1*T2*T3*T4*T5*Tt;

    T01 = T1;
    R01 = T01(1:3,1:3);
    Z1 = T01(1:3,3);
    T1t = T2*T3*T4*T5*Tt;
    P16 = T1t(1:3,4);
    UP1 = cross(Z1,R01*P16);
    J1 = [UP1;Z1];

    T02 = T1*T2;
    R02 = T02(1:3,1:3);
    Z2 = T02(1:3,3);
    T2t = T3*T4*T5*Tt;
    P26 = T2t(1:3,4);
    UP2 = cross(Z2,R02*P26);
    J2 = [UP2;Z2];

    T03 = T1*T2*T3;
    R03 = T03(1:3,1:3);
    Z3 = T03(1:3,3);
    T3t = T4*T5*Tt;
    P36 = T3t(1:3,4);
    UP3 = cross(Z3,R03*P36);
    J3 = [UP3;Z3];

    T04 = T1*T2*T3*T4;
    R04 = T04(1:3,1:3);
    Z4 = T04(1:3,3);
    T4t = T5*Tt;
    P46 = T4t(1:3,4);
    UP4 = cross(Z4,R04*P46);
    J4 = [UP4;Z4];

    T05 = T1*T2*T3*T4*T5;
    R05 = T05(1:3,1:3);
    Z5 = T05(1:3,3);
    T5t = Tt;
    P56 = T5t(1:3,4);
    UP5 = cross(Z5,R05*P56);
    J5 = [UP5;Z5];

    J = [J1 J2 J3 J4 J5];
    K = diag([2940,2350,1956,1960,390]);
%     K = 1e5*diag([38,66,39,5.6,6.6,4.7]);
    K = (180/pi)*60*K;
    Jp = J(1:3,:);
    
    C = Jp/K*Jp';
    D = inv(C);
%     D = D'*D;
    [R, E] = eig(D*D');
    E = diag(E);
    
    norm1 = T(1:3,3);
    Force = norm1;
    deflect = orientation_stiffness(E, R, Force);
    deflect = sqrt(deflect);
%     deflect = (E(1)*E(2)*E(3)*E(4)*E(5)*E(6))^(1/12);
%     deflect = det(C*C')^(1/6);

    
    [~,s,~] = svd(J*J');
    DVal = diag(s);
    
    dexterity = sqrt(prod(DVal(1:5))^(1/5));
%     dexterity = sqrt(DVal(5)/DVal(1));

%     if dexterity <= 1
%         dexterity = 0;
%     else
%         dexterity = 1;
%     end
%     deflect = dexterity * deflect;
%     deflect = sigmoid(dexterity)*sqrt(deflect); 
    deflect = [deflect*sigmoid(dexterity), dexterity];
%     deflect = [deflect, dexterity];
end

function sig = sigmoid(x)
    x = (x - 1)*200;
    sig = 1/(1 + exp(-x));
end

%%  func for calc R matrix, using for calc the orientation stiffness about a spatial curve which a circular cone and ellipsoid intersect
    % arg: n: orientation of normal orientation force
    % return: a 3x3 rotation matrix, 
function R=R_calc(n)
    %n(1)*x+n(2)*y+n(3)*z
    i = [1, 0, (-n(1)-n(2))/n(3)];
    j = cross(n, i);     % the order of n and i means that Y = Z x X, -Y = X x Z
    i=i/norm(i);
    j=j/norm(j);
    R = [i', j', n'];
end