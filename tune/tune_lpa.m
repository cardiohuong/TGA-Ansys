clear; close all; clc;

Q_prev_tbl  = readtable('flowrate.csv');      
P_step_tbl  = readtable('P_step_lpa.csv');    
Q_prev_t = Q_prev_tbl.CardiacPhase;
Q_prev_Q = Q_prev_tbl.MeanCurve;
P_step_t = P_step_tbl.T;
P_step_P2 = P_step_tbl.P_lpa;

t = P_step_t(:);

Q_interp = interp1(Q_prev_t, Q_prev_Q, t, 'spline', 'extrap');
Q_truth  = Q_interp*0.5;

Pin = P_step_P2(:);

x0 = [1.46, 0.0361, 0.382]; 
lb = x0 * 0.5;
ub = x0 * 2.0;

options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','iter', ...
    'MaxIterations',2000, ...
    'MaxFunctionEvaluations',5000);


[xopt, Jopt] = fmincon(@(x)objectiveWK(x, Q_truth, Pin, t), ...
                       x0, [],[],[],[], lb, ub, [], options);

C_opt  = xopt(1);
Rp_opt = xopt(2);
Rd_opt = xopt(3);
fprintf('C  = %.6f\n', C_opt);
fprintf('Rp = %.6f\n', Rp_opt);
fprintf('Rd = %.6f\n', Rd_opt);

Q_sim = simulate_WK(C_opt, Rp_opt, Rd_opt, Pin, t);

figure;
plot(t, Q_truth, 'LineWidth',2); hold on;
%plot(t, Q_1d,'Color','g', 'LineWidth',2); hold on; <- using the 1.44,
%0.036, 0.38 one
plot(t, Q_sim, '--','LineWidth',2);
legend('Q_{CFD}','Q_{Windkessel}');
xlabel('T (s)');
ylabel('Flow-lpa');
grid on;

function J = objectiveWK(x, Q_truth, Pin, t)

    C  = x(1);
    Rp = x(2);
    Rd = x(3);

    Q_sim = simulate_WK(C, Rp, Rd, Pin, t);

    J = sum((Q_truth - Q_sim).^2) / sum(Q_truth.^2);

end
function Q = simulate_WK(C, Rp, Rd, Pin, t)

    N  = length(t);
    Q  = zeros(N,1);
    Pc = zeros(N,1);

    Pc(1) = Pin(1);

    for n = 2:N

        dt = t(n) - t(n-1);

        Pc(n) = ( Pin(n)/Rp + (C/dt)*Pc(n-1) ) ...
              / ( 1/Rp + 1/Rd + C/dt );

        Q(n) = (Pin(n) - Pc(n)) / Rp;

    end

end


%%write data

timestep = 0.034483;
t_new = (t(1):timestep:t(end));
Qsim_rpa = interp1(t, Q_sim, t_new, 'spline');
output_table = table(t_new, Qsim_rpa, ...
                     'VariableNames', {'Time','Flow'});

writetable(output_table, 'Q_sim.csv');
