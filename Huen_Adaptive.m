interval = [0, 10];
v0 = [1/pi*atan(-100)+1/2; 0];
order = 2; %Here we're using Huen's Method
tolerance = 0.01;

h_initial = tolerance*(norm(v0)+1)/norm(f(v0));
h_wiggle = h_initial;
alpha = 0.9; %decay constant for rejecting timestep
beta = 0.9; %decay sconstant for accepting timestep
gamma = 2; %max factor of change for h_wiggle
time = interval(1);
v_o = v0;
v_n = v0;
points = [v0];
reject_counter = 0;
times = [time];

while time < interval(2)
    remaining = interval(2) - time;
    if h_wiggle > remaining
        h_wiggle = remaining;
    end
    if reject_counter > 10
        fprintf('Number of h-wiggle rejects exceeded 10, try adjusting your parameters.')
        return
    end
    halfstep = h_wiggle/2;
    fullstep = h_wiggle;
    v_hat = v_o+fullstep/2*f(v_o)+fullstep/2*f(v_o+fullstep*f(v_o));
    v_nhalf = v_o+halfstep/2*f(v_o)+halfstep/2*f(v_o+halfstep*f(v_o));
    v_n = v_nhalf+halfstep/2*f(v_nhalf)+halfstep/2*f(v_nhalf+halfstep*f(v_nhalf));
    error_h = norm(v_hat-v_n)/(2^order-1);
    %Do the comparison
    left_ineq = error_h/h_wiggle;
    right_ineq = tolerance/(interval(2)-interval(1));
    if left_ineq < right_ineq
        %Advance the solution
        time = time+h_wiggle;
        v_o = v_n;
        points = [points, v_n];
        times = [times, time];
        reject_counter = 0;
        %Pick a new h_wiggle
        %If Err/h is a lot smaller than tolerance/interval, pick a new
        %h_wiggle, otherwise keep it
        if (right_ineq - left_ineq)/right_ineq > 0.1
            h_wiggle = min(beta*h_wiggle*(tolerance/(interval(2)-interval(1))/(error_h/h_wiggle))^(1/order), gamma*h_wiggle);
        end
    else
        %Reject and pick a new h_wiggle
        reject_counter = reject_counter + 1;
        h_wiggle = alpha*h_wiggle*(tolerance/(interval(2)-interval(1))/(error_h/h_wiggle))^(1/order);
    end
end
plot(times, points(1,:), '-o');
hold on
plot(times, points(2,:), '-o');
legend('y1', 'y2')
title(sprintf('Huen`s Method Solution with Adaptive Timestepping'));
fprintf('The solution is y1=%s and y2=%s', num2str(v_n(1)), num2str(v_n(2)));

function [out] = f(v)
    out = [20/(pi*(1+(20*v(2)-100)^2));1];
end