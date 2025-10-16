clear; close all; clc;

% This script formulates and solves a trajectory optimization problem
% underlying a predictive simulation of a hopping (let's hope) robot. 
% A direct collocation method with a backward Euler
% integration scheme is employed to formulate a nonlinear program (NLP)
% from the continuous optimal control problem.
%
% The formulation is deeply inspired from the five-link biped example 
% described in: "Kelly, An Introduction to Trajectory Optimization:
% How to do your own direct collocation (2017), SIAM REVIEW.
% DOI. 10.1137/16M1062569".
%
% Author: Friedl De Groote, Christofer Clemente, Robin Maag

%addpath('C:\Program Files\casadi')
addpath('C:\Users\User\Desktop\3dpredsim_DwarfsAndGiants-master\casadi-windows-matlabR2016a-v3.5.5')

import casadi.*

%% Plot settings.
% You might no want to generate the animation and the figures every time
% you run the code. Feel free to adjust the variables below accordingly.
generate_animation = true;
generate_plots = true;

%% Selection of walking pattern.
% Options:
%   - nominal
%   - on_the_moon
%   - on_mars
selected_gait = 'nominal';

%% Selection of weight of the cost function terms.
% The cost function minimizes the sum of the squared joint torques. The
% contribution of each joint is weighted with a separate weight factor as
% follows: J = w1*T1^2 + w2*T2^2 + w3*T3^2

% We expose those weight factors here so that you can easily change them,
% and see what is the impact on the predicted walking pattern.
w1 = 1; % hip
w2 = 1; % knee
w3 = 1; % ankle
w4 = 1; % toes

%% Model: physical parameters.
% The model consists of a pelvis (p) and a leg. The leg has three segments, 
% femur (f), tibia (t), and foot (m). The ground reaction force is applied
% to the toes.

% Mass of the segments.
mp = 5.0;  %mp = 5.8;
% mp = 5.8;
mf = 2.8;
mt = 0.8; mm = 0.16;

% %27kg kangaroo
 % mp = 17.1; mf = 7.66;
 % mt = 2.04; mm = 0.72;

% Rotational inertia of the segments.
Ip = 0.1; If = 0.005;
It = 0.003; Im = 0.001;

% %27kg kangaroo
% Ip = 0.05; If = 0.027;
% It = 0.016; Im = 0.001;

% Length of the segments.
lp = 0.1; %was 0.1
lf = 0.12;
lt = 0.31; lm = 0.12;

% %27kg kangaroo
% lp = 0.20; lf = 0.17;
% lt = 0.431; lm = 0.17;

% resting length of spring
lso = 0.31;
% lso = 0.30;

% %27kg kangaroo
% lso = 0.431;

% ls: attachment of spring on femur
% lso: resting length of spring
% lc: distance between ankle and heel
% ks: stiffness of spring connecting foot and femur

lc = 0.04; ls = 0.04;
% ks = 10;
ks = 10;

% Distance from segment center of mass to parent joint.
lcf = 0.034;
lct = 0.12; lcm = 0.045;

% %27kg kangaroo
% lcf = 0.05;
% lct = 0.16; lcm = 0.064;
% lc = 0; ls = 0;
% ks = 10; %original value was 200



%damping coefficient
%0.5–1.0 N·m·s/rad for passive joints
%Up to 10 N·m·s/rad for damping in driven models or optimization stability
% For mechanical systems or simulations that need stability, damping may be tuned empirically:
% Start with a small value like 0.1- 1.0
% Increase if your system oscillates
% Decrease if your system feels over-damped or sluggish
c = 1;
% c = 0.1;
%c = 1;




% Gravity.
% Note that on_the_moon and on_mars might be slightly counter-intuitive,
% since there are constraints that prevent the model to 'fly' as you might
% expect from someone walking on the moon.
if strcmp(selected_gait, 'on_the_moon')
    g = 1.62;
elseif strcmp(selected_gait, 'on_mars')
    g = 3.72;
else % on_earth
    g = 9.81;
end

%% Model: dynamics.
% For the simulation to be dynamically consistent, we want those equations of
% motion to be enforced. In pratice, we do that by having path constraints
% in the problem formulation (implicit formulation of the dynamics).
%
% Here, we create a CasADi function that returns the 'model' constraint
% errors based on the model states (q, dq) and controls (ddq, T). This
% function is initialized based on the physical parameters of the model,
% such that if you change those parameters, the equations of motion get
% updated. During the actual optimization, we will impose the contraint
% errors to be null. Note that you can change physical parameters of
% the model (e.g., mass or length), but not for instance add a segment. 
% This would make the equations of motion we generated invalid.
%
% f_getModelConstraintErrors:
% Inputs:
%   - states: segment angles q (5x1)
%   - states: segment angular velocities dq (5x1)
%   - controls: segment angular accelerations ddq (5x1)
%   - controls: joint torques T (5x1)
% Outputs:
%   - model contraint errors (5x1)
f_getModelConstraintErrors = getModelConstraintErrors(...
    mp, mf, mt, mm, ...
    Ip, If, It, Im, ...
    lp, lf, lt, lm, ...
    lcf, lct, lcm, lc, ls, ...
    g, ks, lso, c);

%% Trajectory optimization problem formulation.
% Stride length and time, and mesh size.
% Those are parameters you can play with. If you use a lower mesh size,
% this should increase the accuracy of your simulation, but likely at the
% cost of higher computational time. In practice, if your solution changes
% when lowering the mesh size, it suggests that your current mesh size is 
% not low enough. It is advised to do such type of convergence analysis 
% to make sure you do not misinterpret the results.
speed = 1;         % average forward speed (m/s)
N = 100;          % Number of mesh intervals


% The NLP is formulated using Opti Stack, which is a collection of CasADi
% helper classes that provides a close correspondence between mathematical
% NLP notation and computer code.
% More info here: https://web.casadi.org/docs/#document-opti

% Create opti instance.
opti = casadi.Opti(); 

% Create design variables.
%
% Backward Euler scheme:
% x(t+1) = x(t) + u(t+1)dt
%
% We define the states at N+1 mesh points (starting at k=1).
% We define the controls at N mesh points (starting at k=2)
%
% k=1   k=2   k=3             k=N   k=N+1
% |-----|-----|-----|...|-----|-----|
%
% The dynamic contraints and equations of motion are NOT enforced in k=0.
%
% States.
% Segment angles.
xt_s = opti.variable(1,N+1);   yt_s = opti.variable(1,N+1);   
tp = opti.variable(1,N+1);   tf = opti.variable(1,N+1);   
tt = opti.variable(1,N+1);   tm = opti.variable(1,N+1);
% Segment angular velocities.
dxt_s = opti.variable(1,N+1);   dyt_s = opti.variable(1,N+1);   
dtp = opti.variable(1,N+1);   dtf = opti.variable(1,N+1);   
dtt = opti.variable(1,N+1);   dtm = opti.variable(1,N+1);
% Controls.
% Segment angular accelerations.
ddxt_s = opti.variable(1,N);   ddyt_s = opti.variable(1,N);   
ddtp = opti.variable(1,N);   ddtf = opti.variable(1,N);   
ddtt = opti.variable(1,N);   ddtm = opti.variable(1,N);
% Joint torques.
Th = opti.variable(1,N);     Tk = opti.variable(1,N);     
Ta = opti.variable(1,N);     Fx = opti.variable(1,N);     
Fy = opti.variable(1,N);     Tt = opti.variable(1,N);




% ---- Free stride time and length ----
T = opti.variable(1);                 % stride time [s]
L = opti.variable(1);                 % stride length [m]

% Sensible bounds (adjust as you like)
opti.subject_to(0.15 <= T <= 0.8);    % keep time from going to 0 or huge
opti.subject_to(0.2  <= L <= 2.0);    % avoid degenerate tiny/huge steps

% Initial guesses used also to seed other variables:
T0 = 0.3;    % numeric seeds
L0 = 0.5;

opti.set_initial(T, T0);
opti.set_initial(L, L0);

% Collocation step
dt = T / N;

% T = 0.3;
% L = 0.5;
% dt = T/N;
%time = 0:0.3/100:0.3;

xt = xt_s/100;
yt = yt_s/100;
dxt = dxt_s/100;
dyt = dyt_s/100;
ddxt = ddxt_s/100;
ddyt = ddyt_s/100;


% Set bounds on segment angles (if not otherwise specified, design
% variables will be bounded between +/- Inf).
opti.subject_to(0 <= xt_s <= 1 * 100);
opti.subject_to(-0.01 * 100 <= yt_s <= 1.0 * 100);%CJC
% opti.subject_to(-0.01 * 100 <= yt_s <= 0.5 * 100);
opti.subject_to(0 <= tp <= pi/2);
opti.subject_to(pi/2 <= tf <= pi);
opti.subject_to(0 <= tt <= pi/2);
opti.subject_to(pi/2 <= tm <= pi);

% Set physiological joint limits (no knee hyperextension).
% opti.subject_to(-pi <= q1 - q2 <= 0);
% opti.subject_to(-pi <= q5 - q4 <= 0);
opti.subject_to(0 <= tf - tt <= pi);

opti.subject_to(-50 <= Th <= 50);
opti.subject_to(-50 <= Tk <= 50);
opti.subject_to(-20 <= Ta <= 20);%CJC reduced from 50
opti.subject_to(-20 <= Tt <= 20);%CJC reduced from 50
opti.subject_to(-400 <= Fx <= 400);
opti.subject_to(0 <= Fy <= 600);

% TO DO bounds for other variables

% Set naive initial guess for the segment angles
% (linearly spaced vector between lower and upper bounds).
% When no initial guess is provided, numerical zero is assumed.
xt_init = 0;        xt_final = 2;  
yt_init = 0;        yt_final = 0;
tp_init = pi/4;    tp_final = pi/4;
tf_init = 3*pi/4;    tf_final = 3*pi/4;
tt_init = pi/4;    tt_final = pi/4;
tm_init = 3*pi/4;    tm_final = 3*pi/4;
xtguess = linspace(xt_init, xt_final, N+1);
ytguess = [zeros(1,25) 0.05*ones(1,51) zeros(1,25)];
tpguess = linspace(tp_init, tp_final, N+1);
tfguess = linspace(tf_init, tf_final, N+1);
ttguess = linspace(tt_init, tt_final, N+1);
tmguess = linspace(tm_init, tm_final, N+1);
opti.set_initial(xt_s, xtguess*100);
opti.set_initial(yt_s, ytguess*100);
opti.set_initial(tp, tpguess);
opti.set_initial(tf, tfguess);
opti.set_initial(tt, ttguess);
opti.set_initial(tm, tmguess);

opti.set_initial(dxt_s, 100*L0/T0*ones(1,N+1));
opti.set_initial(dyt_s, 0.1*ones(1,N+1));
opti.set_initial(dtp, ones(1,N+1));
opti.set_initial(dtf, ones(1,N+1));
opti.set_initial(dtt, ones(1,N+1));
opti.set_initial(dtm, ones(1,N+1));

opti.set_initial(ddxt_s, 100*ones(1,N));
opti.set_initial(ddyt_s, 100*ones(1,N));
opti.set_initial(ddtp, ones(1,N));
opti.set_initial(ddtf, ones(1,N));
opti.set_initial(ddtt, ones(1,N));
opti.set_initial(ddtm, ones(1,N));

% Initialize the cost function (J).
J = 0;

% Loop over mesh points.
% k=1   k=2   k=3             k=N   k=N+1
% |-----|-----|-----|...|-----|-----|
for k=1:N
    % States at mesh point k.
    % Segment angles.
    xtk = xt(:,k);     ytk = yt(:,k);     tpk = tp(:,k);
    tfk = tf(:,k);     ttk = tt(:,k);     tmk = tm(:,k);     
    % Segment angular velocities.
    dxtk = dxt(:,k);     dytk = dyt(:,k);     dtpk = dtp(:,k);
    dtfk = dtf(:,k);     dttk = dtt(:,k);     dtmk = dtm(:,k);
    
    % States at mesh point k+1.
    % Segment angles.
    xtk_plus = xt(:,k+1);     ytk_plus = yt(:,k+1);     tpk_plus = tp(:,k+1);
    tfk_plus = tf(:,k+1);     ttk_plus = tt(:,k+1);     tmk_plus = tm(:,k+1); 
    % Segment angular velocities.
    dxtk_plus = dxt(:,k+1);     dytk_plus = dyt(:,k+1);     dtpk_plus = dtp(:,k+1);
    dtfk_plus = dtf(:,k+1);     dttk_plus = dtt(:,k+1);     dtmk_plus = dtm(:,k+1);
    
    % Controls at mesh point k+1.
    % (Remember that controls are defined from k=2, so 'mesh point k+1 for
    % the states corresponds to element k for the controls', which is why
    % we use k and not k+1 here).
    % Segment angular accelerations.
    ddxtk_plus = ddxt(:,k);     ddytk_plus = ddyt(:,k);     ddtpk_plus = ddtp(:,k);
    ddtfk_plus = ddtf(:,k);     ddttk_plus = ddtt(:,k);     ddtmk_plus = ddtm(:,k);
    % Joint torques.
    Thk_plus = Th(:,k);     Tkk_plus = Tk(:,k);     Tak_plus = Ta(:,k);     
    Fxk_plus = Fx(:,k);     Fyk_plus = Fy(:,k);     Ttk_plus = Tt(:,k);
       
    % Stack states at mesh points k and k+1.
    Xk = [xtk; ytk; tpk; tfk; ttk; tmk; ...
          dxtk; dytk; dtpk; dtfk; dttk; dtmk];
    Xk_plus = [xtk_plus; ytk_plus; tpk_plus; tfk_plus; ttk_plus; tmk_plus; ...
               dxtk_plus; dytk_plus; dtpk_plus; dtfk_plus; dttk_plus; dtmk_plus];
    
    % Stack state derivatives at mesh points k+1.
    Uk_plus = [dxtk_plus; dytk_plus; dtpk_plus; dtfk_plus; dttk_plus; dtmk_plus; ...
               ddxtk_plus; ddytk_plus; ddtpk_plus; ddtfk_plus; ddttk_plus; ddtmk_plus];
       
    % Path constraints - dynamic constraints.
    % The function eulerIntegrator returns the error in the dynamics.
    % We impose this error to be null (i.e., dqdt* = dqdt and
    % ddqdt* = ddqdt, where * indicates the approximated state derivatives
    % computed based on the integration scheme and no * represents the
    % actual states or controls. Both should match - collocation).
    % The integration is performed using a backward Euler scheme
    % (see eulerIntegrator.m)
    opti.subject_to(eulerIntegrator(Xk, Xk_plus, Uk_plus, dt) == 0);
     
    % Path constraints - model constraints (implicit skelton dynamics).
    % We impose this error to be null (i.e., f(q, dq, ddq, T) = 0).
    modelConstraintErrors = f_getModelConstraintErrors(...
        tpk_plus,tfk_plus,ttk_plus,tmk_plus,...
        dtpk_plus,dtfk_plus,dttk_plus,dtmk_plus,...
        ddxtk_plus,ddytk_plus,ddtpk_plus,ddtfk_plus,ddttk_plus,ddtmk_plus,...
        Thk_plus,Tkk_plus,Tak_plus,Ttk_plus,Fxk_plus,Fyk_plus);
    opti.subject_to(modelConstraintErrors == 0);
    
    if k == 1
    Tkk = Tk(:,end);
    Thk = Th(:,end);
    Tak = Ta(:,end);
    else
    Tkk = Tk(:,k-1);
    Thk = Th(:,k-1);
    Tak = Ta(:,k-1);
    end

    opti.subject_to((Thk_plus - Thk) / dt <= 100); %constrain kknee torque
    opti.subject_to((Tkk_plus - Tkk) / dt <= 100); %constrain kknee torque
    opti.subject_to((Tak_plus - Tak) / dt <= 100); %constrain kknee torque

    % Path constraints
    % if k > 30 && k < 70
    %     %opti.subject_to(ytk >= 0.05);
    %     opti.subject_to(ytk >= 0.05);
    % elseif k < 10 || k > 90
    %     opti.subject_to(ytk == 0);
    % else
    %     opti.subject_to(ytk >= 0);
    % end

% Path constraints
    if k > 30 && k < 70
        %opti.subject_to(ytk >= 0.05);
        opti.subject_to(ytk >= 0.03);
    elseif k < 14 || k > 86
        opti.subject_to(ytk == 0);
    else
        opti.subject_to(ytk >= 0);
    end



    % if k == 50
    %     opti.subject_to(ytk >= 0.05);
    % end
    % 
    % Make sure that we can only apply a ground reaction force when the
    % toes are on the ground (i.e. vertical position is 0)
    opti.subject_to(ytk_plus*Fxk_plus == 0);
    opti.subject_to(ytk_plus*Fyk_plus == 0);

    %opti.subject_to(Fxk_plus - 100*dxtk_plus*Fyk_plus == 0);
    
    % --- Simple Coulomb friction cone at contact ---
    mu = 0.8;                             % friction coefficient
    opti.subject_to(Fyk_plus >= 0);       % no "pulling" from the ground
    opti.subject_to(-mu*Fyk_plus <= Fxk_plus <= mu*Fyk_plus);
    

    % Cost function.
    % Minimize the weighted sum of the squared joint torques.
    J = J + (w1*Thk_plus.^2 + w2*Tkk_plus.^2 + w3*Tak_plus.^2 + w4*Ttk_plus.^2)*dt;
    % Penalize (with low weight) segment angular accelerations for
    % regularization purposes.
    J = J + 1e-4*(ddxtk_plus.^2 + ddytk_plus.^2 + ddtpk_plus.^2 + ...
        ddtfk_plus.^2 + ddttk_plus.^2 + ddtmk_plus.^2)*dt;
    
end

% J = J/(xt(end)-xt(1));

% Boundary constraints - periodic gait.
opti.subject_to(xt(1) == 0);
% opti.subject_to(yt(1) == 0);
opti.subject_to(yt(1) - yt(end) == 0);
opti.subject_to(tp(1) - tp(end) == 0);
opti.subject_to(tf(1) - tf(end) == 0);
opti.subject_to(tt(1) - tt(end) == 0);
opti.subject_to(tm(1) - tm(end) == 0);

opti.subject_to(dxt(1) - dxt(end) == 0);
opti.subject_to(dyt(1) - dyt(end) == 0);
opti.subject_to(dtp(1) - dtp(end) == 0);
opti.subject_to(dtf(1) - dtf(end) == 0);
opti.subject_to(dtt(1) - dtt(end) == 0);
opti.subject_to(dtm(1) - dtm(end) == 0);

% Boundary constraints - gait speed.
% opti.subject_to(xt(end)-xt(1) - T*speed == 0);
opti.subject_to(xt(end)-xt(1) - L == 0);

% Set cost function.
opti.minimize(J);

% Create an NLP solver.
optionssol.ipopt.max_iter = 1e5;
optionssol.ipopt.tol = 1e-6;
optionssol.ipopt.constr_viol_tol = 1e-6;
optionssol.ipopt.nlp_scaling_method = 'none';


% Termination tolerances (looser = easier to ?declare success?)
optionssol.ipopt.tol                 = 1e-3;   % was 1e-6
optionssol.ipopt.constr_viol_tol     = 1e-3;   % was 1e-6
optionssol.ipopt.compl_inf_tol       = 1e-3;   % complementarity tolerance

% ?Acceptable? (early-exit) tolerances ? very handy when exact tol is hard
optionssol.ipopt.acceptable_tol                  = 1e-2;
optionssol.ipopt.acceptable_constr_viol_tol      = 1e-2;
optionssol.ipopt.acceptable_compl_inf_tol        = 1e-2;
optionssol.ipopt.acceptable_obj_change_tol       = 1e-4;
optionssol.ipopt.acceptable_iter                 = 10;

% Slightly relax bound handling (helps when stuck on bounds)
optionssol.ipopt.bound_relax_factor   = 1e-6;   % default 1e-8

% Diagnostic output (optional but useful)
optionssol.ipopt.print_level          = 5;      % 0?12


opti.solver('ipopt',optionssol);

% Solve the NLP.
sol = opti.solve();

% %% Extract the optimal design variables.
% % Optimal segment angles.
xt_opt = sol.value(xt);
yt_opt = sol.value(yt);
tp_opt = sol.value(tp);
tf_opt = sol.value(tf);
tt_opt = sol.value(tt);
tm_opt = sol.value(tm);

% T_opt = sol.value(T);
% t_opt = 0: T_opt/N: T_opt;
T_opt = sol.value(T);
L_opt = sol.value(L);
dt_opt = T_opt / N;

% Build time vectors numerically
time   = linspace(0, T_opt, N+1);  % for states
time_u = time(1:end-1);            % for controls/forces

Fx_opt = sol.value(Fx);
Fy_opt = sol.value(Fy);

% % Optimal segment angular velocities.
% dq1_opt = sol.value(dq1);
% dq2_opt = sol.value(dq2);
% dq3_opt = sol.value(dq3);
% dq4_opt = sol.value(dq4);
% dq5_opt = sol.value(dq5);
% % Optimal segment angular accelerations.
% ddq1_opt = sol.value(ddq1);
% ddq2_opt = sol.value(ddq2);
% ddq3_opt = sol.value(ddq3);
% ddq4_opt = sol.value(ddq4);
% ddq5_opt = sol.value(ddq5);

% % Optimal joint torques.
Th_opt = sol.value(Th);
Tk_opt = sol.value(Tk);
Ta_opt = sol.value(Ta);
Tt_opt = sol.value(Tt);

%% Generate an animation.
if generate_animation == true
    jointPositions_opt = getJointPositions(...
        lf,lm,lp,lt,tf_opt,tm_opt,tp_opt,tt_opt,xt_opt,yt_opt)';
    generateAnimation(jointPositions_opt, dt_opt, L_opt);
end
 
%% Plots.
if generate_plots == true
    % Generalized coordinates
    figure()
    subplot(231)
    plot(time,xt_opt,'k','LineWidth',2)
    xlabel('time [s]')
    ylabel('position [m]')
    title('x_T')

    subplot(232)
    plot(time,yt_opt,'k','LineWidth',2)
    xlabel('time [s]')
    ylabel('position [m]')
    title('y_T')

    subplot(233)
    plot(time,180/pi*tp_opt,'k','LineWidth',2)
    xlabel('time [s]')
    ylabel('angle [^o]')
    title('pelvis')

    subplot(234)
    plot(time,180/pi*tf_opt,'k','LineWidth',2)
    xlabel('time [s]')
    ylabel('angle [^o]')
    title('femur')

    subplot(235)
    plot(time,180/pi*tt_opt,'k','LineWidth',2)
    xlabel('time [s]')
    ylabel('angle [^o]')
    title('tibia')

    subplot(236)
    plot(time,180/pi*tm_opt,'k','LineWidth',2)
    xlabel('time [s]')
    ylabel('angle [^o]')
    title('foot')
    
    % Ground reaction forces
    figure()
    subplot(121)
    plot(time(1:end-1),Fx_opt,'k','LineWidth',2)
    xlabel('time [s]')
    ylabel('force [N]')
    title('F_x')

    subplot(122)
    plot(time(1:end-1),Fy_opt,'k','LineWidth',2)
    xlabel('time [s]')
    ylabel('force [N]')
    title('F_y')

    % Joint torques
    figure()
    subplot(221)
    plot(time(1:end-1),Th_opt,'k','LineWidth',2)
    xlabel('time [s]')
    ylabel('hip torque [Nm]')
   
    subplot(222)
    plot(time(1:end-1),Tk_opt,'k','LineWidth',2)
    xlabel('time [s]')
    ylabel('knee torque [Nm]')

    subplot(223)
    plot(time(1:end-1),Ta_opt,'k','LineWidth',2)
    xlabel('time [s]')
    ylabel('ankle torque [Nm]')
   
    subplot(224)
    plot(time(1:end-1),Tt_opt,'k','LineWidth',2)
    xlabel('time [s]')
    ylabel('mtp torque [Nm]')
end

% %% Maximum torque.
% max_torque=max(abs([T1_opt T2_opt T3_opt T4_opt T5_opt]));
% 
% disp(['The maximum torque is ', num2str(max_torque), ' Nm. '...
%       'Try to make it lower by playing with the cost term weights.'])
