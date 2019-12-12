clear;clc
NUM = 50000;
joint_angle = joint_angle_generator(NUM);
high_point = zeros(1000,3);
low_point = zeros(1000,3);
j = 1;
k = 1;

avg_cond = 0;
avg_mani = 0;

for i = 1:1:NUM
    T_matrix = T_matrix_calc(joint_angle(i, :));
    T_matrix =  T_matrix(1:3,4);
    dex = cond_calc(joint_angle(i, :));
    cond = dex(1);
    mani = dex(2);
    if mani <= 1
        low_point(j,:) = T_matrix;
        j = j+1;
    else
        high_point(k,:) = T_matrix;
        k = k +1;
    end
    avg_cond = avg_cond + cond;
    avg_mani = avg_mani + mani;
end

figure(1)
scatter3(low_point(:, 1), low_point(:, 2), low_point(:, 3),'r');
figure(2)
scatter3(high_point(:, 1), high_point(:, 2), high_point(:, 3),'g');
figure(3)
scatter3(low_point(:, 1), low_point(:, 2), low_point(:, 3),'r');
hold on 
scatter3(high_point(:, 1), high_point(:, 2), high_point(:, 3),'g');


jsvd_calc(joint_angle(1, :));

%% func for jth generating
    % args: NUM: number of end points, seed: the random seed
    % return: a matrix with shape of NUMx5 
function joint_angle = joint_angle_generator(NUM)
    range = [pi*2, (145/180*pi), (230/180*pi), pi, (240/180*pi)];
    bias = [-pi, -(90/180*pi), -(70/180*pi), -pi/2, -deg2rad(120)];
%     range = [0, (145/180*pi), (230/180*pi), 0, (240/180*pi)];
%     bias = [0, -(90/180*pi), -(50/180*pi), 0, -deg2rad(120)];
    
%     rng(100, 'twister');
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
function T_matrix = T_matrix_calc(th)
    th1=th(1);th2=th(2);th3=th(3);th4=th(4);th5=th(5);

    % l1=1.0;l2=0.4;l6=0.5;e=0.37;
    l1=1.05;l2=0.5;l6=0.39;e=0.15;
    l3 = 1.3;
    % l3 = 1.3;
    % l3 = (para2/(para2+1))*2.6;
    l4 = 1.2/2;
    % l4 = (2.6-l3)/2;
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
function J = Jacobian_matrix(th)
    th1=th(1);th2=th(2);th3=th(3);th4=th(4);th5=th(5);

    % l1=1.0;l2=0.4;l6=0.5;e=0.37;
    l1=1.05;l2=0.5;l6=0.39;e=0.15;
    l3 = 1.3;
    % l3 = 1.3;
    % l3 = (para2/(para2+1))*2.6;
    l4 = 1.2/2;
    % l4 = (2.6-l3)/2;
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
end
    
    
function cond = cond_calc(th)
    J = Jacobian_matrix(th);
    Jp = J(1:3,:);
    [~,s,~] = svd(J*J');
    DVal = diag(s);
    cond = sqrt(DVal(5)/DVal(1));
    mani = sqrt(prod(DVal(1:5))^(1/5));
    cond = [cond, mani];
end

function jsvd_calc(th)
    J = Jacobian_matrix(th);
    Jp = J(1:3,:);
    [u,s,v] = svd(J*J');
end

