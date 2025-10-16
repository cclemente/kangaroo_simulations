% This script generates the equations of motion of the hopping robot using 
% MATLAB's symbolic toolbox.
% 
% Author: Friedl De Groote

%% Symbolic variables.
% xp is x position of pelvis, xy is y position of pelvis, tp is orientation 
% of pelvis, th is orientation of femur, tk is orientation of tibia, tm is orientation of foot.
% d...: velocities (states).
% d...: accelerations (slack controls).
% Th-Tk-Ta: joint torques (controls) (we might constrain Ta to be 0).
syms xt yt tm tt tf tp dxt dyt dtm dtt dtf dtp ddxt ddyt ddtm ddtt ddtf ddtp Th Tk Ta Tt Fx Fy c 'real';
% m...: masses (p - pelvis, f - femur, t - tibia, m - foot)
% l...: lengths (pelvis com - hip; hip - knee; knee - ankle; ankle - toe, i.e. contact with ground).
% I...: moments of inertia.
% lc...: distances between segment's parent joint and CoM.
% ls: attachment of spring on femur
% lso: resting length of spring
% lc: distance between ankle and heel
% ks: stiffness of spring connecting foot and femur
% g: gravity.
% c: damping coefficient
syms mp mf mt mm g lp lf lt lm Ip If It Im lcf lct lcm ls lso lc ks 'real';

%% Equations of motion.
disp('Deriving the equations of motion symbolically...')

% Positions of joints, segment COMs, and application points of forces.
T = [xt;yt;0];

A = T + lm * [cos(tm); sin(tm); 0];
G_m = T + (lm - lcm) * [cos(tm); sin(tm); 0];

K = A + lt * [cos(tt); sin(tt); 0];
G_t = A + (lt - lct) * [cos(tt); sin(tt); 0];

H = K + lf * [cos(tf); sin(tf); 0];
G_f = K + (lf - lcf) * [cos(tf); sin(tf); 0];

G_p = H + lp * [cos(tp); sin(tp); 0];

C = T +  (lm + lc) * [cos(tm); sin(tm); 0];

S = K + ls * [cos(tf); sin(tf); 0];

% Column vectors.
P = [T; A; K; H; G_p];
% G = [G_p; G_f; G_t; G_m];
q = [xt; yt; tm; tt; tf; tp];
dq = [dxt; dyt; dtm; dtt; dtf; dtp];
ddq = [ddxt; ddyt; ddtm; ddtt; ddtf; ddtp];

% Get velocities.
% dH = jacobian(H,q)*dq; 
% dK = jacobian(K,q)*dq; 
% dA = jacobian(A,q)*dq; 

dG_p = jacobian(G_p,q)*dq;
dG_f = jacobian(G_f,q)*dq;
dG_t = jacobian(G_t,q)*dq;
dG_m = jacobian(G_m,q)*dq;

% Get accelerations.
ddG_p = jacobian(dG_p,[q;dq])*[dq; ddq]; 
ddG_f = jacobian(dG_f,[q;dq])*[dq; ddq]; 
ddG_t = jacobian(dG_t,[q;dq])*[dq; ddq]; 
ddG_m = jacobian(dG_m,[q;dq])*[dq; ddq]; 

%length of the spring
dsc = sqrt( ( S(1)-C(1) )^2 + ( S(2)-C(2) )^2 );

% Joint angles and angular velocities
th = pi + tp - tf;
tk = pi - tf + tt;
ta = pi + tt - tf;

dth = dtp - dtf; %joint velocity at the hip
dtk = -dtf + dtt; %joint vel at the knee
dta = dtt - dtf;



% Equations for angular momentum balance about each joint to get five
% independent equations and construct equations of motion.

i = [1; 0; 0];
j = [0; 1; 0];
k = [0; 0; 1];


% equation of motion of pelvis with respect to hip
eq1 = [0; 0; Th] ...                                                % hip torque
    + k.*cross((G_p - H),(-mp*g*j)) ...                             %TORQUE DUE TO GRAVITY
        - ddtp*Ip*k ...                                             %ROTATIONAL INERTIA
        - k.*cross((G_p - H),(mp*ddG_p))...                         %INERTIAL TORQUE DUE TO LINEAR ACCELERATION OF THE com
        - [0; 0; c*dth];                                            % damping torque about z-axis

eq1 = eq1(3);

% equation of motion of pelvis and femur with respect to knee
% eq2 = [0; 0; Tk] ...                            % knee torque
%     + k.*cross((G_p - K),(-mp*g*j)) ...
%     + k.*cross((G_f - K),(-mf*g*j)) ...
%     + k.*cross((S-K), (ks * (dsc - lso)*((C - S)/dsc)))...
%     - (ddtp*Ip*k +ddtf*If*k + k.*cross((G_p - K),(mp*ddG_p))+ k.*cross((G_f - K),(mf*ddG_f)));

%ddtk corresponds to the angular velocity of the femur relative to the knee. If you're tracking pelvis and femur separately, and you already have:
%ddtp: angular velocity of pelvis
%ddtf: angular velocity of femur

% equation of motion of pelvis and femur with respect to knee
eq2 = [0; 0; Tk] ...                                                % knee torque
    + k.*cross((G_p - K), (-mp*g*j)) ...                            %torque due to gravity acting on the pelvis mass (mp), measured about the knee joint (K).
    + k.*cross((G_f - K), (-mf*g*j)) ...                            %Same as above, but for the femur mass (mf).
    + k.*cross((S - K), (ks * (dsc - lso) * ((C - S)/dsc))) ...     %This is the spring/muscle torque about the knee:
    - (ddtp*Ip*k + ddtf*If*k ...                                    %Torque from angular acceleration of the pelvis and femur
       + k.*cross((G_p - K), (mp*ddG_p)) ...                        %Inertial torque from linear acceleration of pelvis COM  
       + k.*cross((G_f - K), (mf*ddG_f)))...                        %Same, but for the femur
       + [0; 0; c * dtk];                                           %damping for the knee


eq2 = eq2(3);

% equation of motion of pelvis, femur, tibia with respect to ankle
eq3 = [0; 0; Ta] ...                                                % ankle torque
    + k.*cross((G_p - A),(-mp*g*j)) ...                             % torque due to gravity acting on the pelvis mass (mp), measured about the ankle joint (A).
    + k.*cross((G_f - A),(-mf*g*j)) ...                             %Same as above, but for the femur mass (mf).
    + k.*cross((G_t - A),(-mt*g*j)) ...                             %Same as above, but for the tibia mass (mt).
    + k.*cross((C - A), (ks * (dsc - lso)*((C - S)/dsc))) ...       %This is the spring/muscle torque about the knee:
    - (ddtp*Ip*k +ddtf*If*k +ddtt*It*k ...                          %Torque from angular acceleration of the pelvis, femur and tibia
        + k.*cross((G_p - A),(mp*ddG_p)) ...                        %Inertial torque from linear acceleration of pelvis COM
        + k.*cross((G_f - A),(mf*ddG_f)) ...                        %Inertial torque from linear acceleration of femur COM
        + k.*cross((G_t - A),(mt*ddG_t))) ...                       %Inertial torque from linear acceleration of tibia COM 
    - [0; 0; c * dta];                                              %damping for the knee
eq3 = eq3(3);


% Toe compliance torque (spring-damper)
%Tt = -kt * (theta_t - theta_f - theta0) - ct * (ddt_t - ddt_f);
% Parameter	Typical Range
% kt (toe stiffness)	2 – 10 N·m/rad
% ct (toe damping)	0.05 – 1 N·m·s/rad
% theta0	0 (neutral resting angle)
%Start with kt = 5, ct = 0.2, and adjust based on behavior (e.g., too stiff, oscillatory, etc.)

% equation of motion of pelvis, femur, tibia, foot with respect to toe
eq4 = [0; 0; Tt] ...
    + k.*cross((G_m - T),(-mm*g*j)) ...
    + k.*cross((G_t - T),(-mt*g*j)) ...
    + k.*cross((G_f - T),(-mf*g*j)) ...
    + k.*cross((G_p - T),(-mp*g*j)) ...
    - (ddtm*Im*k +ddtt*It*k +ddtf*If*k +ddtp*Ip*k ...
        + k.*cross((G_m - T),(mm*ddG_m))...
        + k.*cross((G_t - T),(mt*ddG_t))...
        + k.*cross((G_f - T),(mf*ddG_f))...
        + k.*cross((G_p - T),(mp*ddG_p)))...
        - [0; 0; c * dtt];
eq4 = eq4(3);

eq56 = (-mm*g*j) + (-mt*g*j) + (-mf*g*j) + (-mp*g*j) + [Fx; Fy; 0] ...
    - (mm*ddG_m + mt*ddG_t + mf*ddG_f + mp*ddG_p);

eq5 = eq56(1);
eq6 = eq56(2);

eq_systemDynamics = simplify([eq1; eq2; eq3; eq4; eq5; eq6]);
f_eq_systemDynamics = matlabFunction(eq_systemDynamics,'File','getSystemDynamics.m');

%% Joint positions and velocities as well as relative joint angles and angular velocities
% disp('Creating symbolic functions to compute joint positions and velocities as well as relative joint angles and angular velocities...')
%P = [T; A; K; H; G_p];
jointPositions = P([1 2 4 5 7 8 10 11 13 14],:);
joint_angle_Velocities = [dth, dtk, dta];

P_fcn = matlabFunction(jointPositions,'File','getJointPositions.m');
dP_fcn = matlabFunction(joint_angle_Velocities,'File','getJointAngularVelocities.m');

% % Relative joint angles.
% q_ANK = q1;
% q_stanceKNEE = q1 - q2;
% q_stanceHIP = q2 - q3;
% q_swingHIP = q4 - q3;
% q_swingKNEE = q5 - q4;
% relativeJointAngles = [q_ANK; q_stanceKNEE; q_stanceHIP; q_swingHIP; q_swingKNEE];
% Prel_fcn = matlabFunction(relativeJointAngles,'File','getRelativeJointAngles.m');
% 
% % Relative joint velocities.
% dq_ANK = dq1;
% dq_stanceKNEE = dq1 - dq2;
% dq_stanceHIP = dq2 - dq3;
% dq_swingHIP = dq4 - dq3;
% dq_swingKNEE = dq5 - dq4;
% relativeJointAngularVelocities = [dq_ANK; dq_stanceKNEE; dq_stanceHIP; dq_swingHIP; dq_swingKNEE];
% dPrel_fcn = matlabFunction(relativeJointAngularVelocities,'File','getRelativeJointAngularVelocities.m');
