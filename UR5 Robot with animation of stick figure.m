clc;
clear;
%% Define DH parameters for UR 5 robot(a, alpha, d, theta)
a = [0 0 0.5 0.5 0 0]; % Link length
alpha = [0 pi/2 0 0 pi/2 -pi/2]; % Link twist
d = [0.1 0 0 0.1 0.1 0.1]; % Link offset
theta = [0, 0, 0, 0, 0, 0]; % Joint angle
%% Performing inverse Kinematics for the UR5 robot
% 1. Set convergence criteria
max_iterations = 100;
epsilon = 1e-6;

% 2. Initialize joint angles (initial guess)

theta_guess = zeros(1, 6); 
final_joint_angles = zeros(10, 6);
% 3. Get Initial Position;
homeposition = forward_kinematics(theta_guess, a, alpha, d, theta);
% 4. Center point of the circle for State space trajectory
center = homeposition';

radius = 0.2; % Radius of the circle
num_points = 100; % Number of waypoints (adjust as needed)
t = linspace(0, 2*pi, num_points); % Angle increment to generate waypoints

% 5. Calculate waypoints for the path
x_points = repmat(center(1), 1, num_points); % All points have the same x-coordinate
y_points = center(2) + radius * cos(t);
z_points = center(3) + radius * sin(t);

% 6. Target end-effector waypoints of moving in circle:
for j=1:num_points
    target_pos = [x_points(j),y_points(j), z_points(j)]';

    % Newton's method iteration loop
    for iter = 1:max_iterations
        current_pos = forward_kinematics(theta_guess, a, alpha, d, theta); % Compute end-effector position using forward kinematics
        error = target_pos - current_pos;  % Compute error between current and target positions     
        % Check convergence
        if norm(error) < epsilon
            disp('Converged!');
            break;
        end                
        J = compute_jacobian(theta_guess, a, alpha, d, theta); % Compute Jacobian matrix       
        % Compute change in joint angles using Newton's method
        delta_theta = pinv(J) * error; % Pseudo-inverse to handle singularities                
        theta_guess = theta_guess + delta_theta'; % Update joint angles
    end
    final_joint_angles(j,:) = theta_guess;
end
disp('Final joint angles:');
disp(final_joint_angles);

%% Plot and Animation of UR5 robot
% Initialize figure
figure;
hold on;
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal;
plot3(x_points, y_points, z_points, 'b', 'LineWidth', 1); % Plot circle
view(3);

% Draw the robot continuously
for k = 1:num_points
    theta_f = final_joint_angles(k, :);
    cla; % Clear current axes
    plot3(x_points, y_points, z_points, 'b', 'LineWidth', 1); % Plot circle again to keep it visible
    draw_robot(a, alpha, d, theta_f);
    pause(0.05); % Pause for 0.01 seconds to control the animation speed
    drawnow; % Update the figure window to show the animation
end

% Function to compute forward kinematics
function pos = forward_kinematics(theta, a, alpha, d, theta_offset)
    num_links = length(theta);
    pos = eye(4);
    for i = 1:num_links
        A_i = [
            cos(theta(i) + theta_offset(i)), -sin(theta(i) + theta_offset(i)) * cos(alpha(i)), sin(theta(i) + theta_offset(i)) * sin(alpha(i)), a(i) * cos(theta(i) + theta_offset(i));
            sin(theta(i) + theta_offset(i)), cos(theta(i) + theta_offset(i)) * cos(alpha(i)), -cos(theta(i) + theta_offset(i)) * sin(alpha(i)), a(i) * sin(theta(i) + theta_offset(i));
            0, sin(alpha(i)), cos(alpha(i)), d(i);
            0, 0, 0, 1
        ];
        pos = pos * A_i;
    end
    pos = pos(1:3, 4); % Extract position
end

% Function to compute Jacobian matrix
function J = compute_jacobian(theta, a, alpha, d, theta_offset)
    num_links = length(theta);
    J = zeros(3, num_links);
    for i = 1:num_links
        % Perturb joint angle theta_i
        theta_perturbed = theta;
        theta_perturbed(i) = theta_perturbed(i) + 1e-6; % Small perturbation        
        % Compute end-effector position for perturbed and current joint angles
        pos_perturbed = forward_kinematics(theta_perturbed, a, alpha, d, theta_offset);
        pos_current = forward_kinematics(theta, a, alpha, d, theta_offset);        
        % Compute column i of the Jacobian matrix
        J(:, i) = (pos_perturbed - pos_current) / 1e-6; % Finite difference approximation
    end
end

function draw_robot(a, alpha, d, theta)
    % Initialize transformation matrix
    T = eye(4);
    positions = zeros(3, length(a) + 1);
    
    % Loop through each link
    for i = 1:length(a)
        % Compute transformation matrix for current link
        A = [cos(theta(i)), -sin(theta(i)) * cos(alpha(i)), sin(theta(i)) * sin(alpha(i)), a(i) * cos(theta(i));
             sin(theta(i)), cos(theta(i)) * cos(alpha(i)), -cos(theta(i)) * sin(alpha(i)), a(i) * sin(theta(i));
             0, sin(alpha(i)), cos(alpha(i)), d(i);
             0, 0, 0, 1];
        
        % Update overall transformation matrix
        T = T * A;
        % Extract position of the end of current link
        positions(:, i + 1) = T(1:3, 4);
    end
    
    % Plot the robot
    hold on;
    grid on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis equal;
    
    % Plot coordinate frames
    for i = 1:length(a)
        % Plot coordinate frame
        if i == length(a)
            % End effector frame
            frame_length = 0.1;
            x_axis = T(1:3, 1) * frame_length + positions(:, i + 1);
            y_axis = T(1:3, 2) * frame_length + positions(:, i + 1);
            z_axis = T(1:3, 3) * frame_length + positions(:, i + 1);
            plot3([positions(1, i + 1) x_axis(1)], [positions(2, i + 1) x_axis(2)], [positions(3, i + 1) x_axis(3)], 'r', 'LineWidth', 2);
            plot3([positions(1, i + 1) y_axis(1)], [positions(2, i + 1) y_axis(2)], [positions(3, i + 1) y_axis(3)], 'g', 'LineWidth', 2);
            plot3([positions(1, i + 1) z_axis(1)], [positions(2, i + 1) z_axis(2)], [positions(3, i + 1) z_axis(3)], 'b', 'LineWidth', 2);
        else
            % Plot coordinate frame
            frame_length = 0.05;
            x_axis = T(1:3, 1) * frame_length + positions(:, i + 1);
            y_axis = T(1:3, 2) * frame_length + positions(:, i + 1);
            z_axis = T(1:3, 3) * frame_length + positions(:, i + 1);
            plot3([positions(1, i + 1) x_axis(1)], [positions(2, i + 1) x_axis(2)], [positions(3, i + 1) x_axis(3)], 'r', 'LineWidth', 2);
            plot3([positions(1, i + 1) y_axis(1)], [positions(2, i + 1) y_axis(2)], [positions(3, i + 1) y_axis(3)], 'g', 'LineWidth', 2);
            plot3([positions(1, i + 1) z_axis(1)], [positions(2, i + 1) z_axis(2)], [positions(3, i + 1) z_axis(3)], 'b', 'LineWidth', 2);
            
            % Connect to previous joint
            plot3([positions(1, i) positions(1, i + 1)], [positions(2, i) positions(2, i + 1)], [positions(3, i) positions(3, i + 1)], 'k');
        end
    end
    
    % Adjust axis limits to ensure all frames are fully visible
    max_lim = max(positions, [], 2);
    min_lim = min(positions, [], 2);
    xlim([min_lim(1) - 0.3, max_lim(1) + 0.3]);
    ylim([min_lim(2) - 0.5, max_lim(2) + 0.5]);
    zlim([min_lim(3) - 0.3, max_lim(3) + 0.3]);
end
