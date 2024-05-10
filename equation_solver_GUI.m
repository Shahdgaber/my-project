function equation_solver_GUI
    % Create the main figure
    fig = uifigure('Name', 'Equation Solver', 'Position', [100 100 400 400]);

    % Create input fields for equations
    eq1_edit = uieditfield(fig, 'text', 'Position', [50 350 300 30], 'Value', 'Enter Equation 1:');
    eq2_edit = uieditfield(fig, 'text', 'Position', [50 300 300 30], 'Value', 'Enter Equation 2:');
    eq3_edit = uieditfield(fig, 'text', 'Position', [50 250 300 30], 'Value', 'Enter Equation 3:');
    eq4_edit = uieditfield(fig, 'text', 'Position', [50 200 300 30], 'Value', 'Enter Equation 4:');

    % Create dropdown menu for method selection
    method_dropdown = uidropdown(fig, 'Items', {'Elliptic', 'Parabolic', 'Hyperbolic', 'ODE'}, ...
        'Position', [50 150 150 30], 'Value', 1);

    % Create button to solve equations
    solve_button = uibutton(fig, 'push', 'Text', 'Solve', ...
        'Position', [230 150 120 30], 'ButtonPushedFcn', @solveEquations);

    % Create textbox to display results
    result_textbox = uitextarea(fig, 'Position', [50 50 300 80]);

    % Callback function to solve equations
    function solveEquations(~, ~)
        % Get input equations
        eq1 = sym(eq1_edit.Value);
        eq2 = sym(eq2_edit.Value);
        eq3 = sym(eq3_edit.Value);
        eq4 = sym(eq4_edit.Value);

        % Get selected method
        method = method_dropdown.Value;

        % Call corresponding solver function
        switch method
            case 'Elliptic'
                result = ellipticSolver(eq1, eq2, eq3, eq4);
            case 'Parabolic'
                result = parabolicSolver(eq1, eq2, eq3, eq4);
            case 'Hyperbolic'
                result = hyperbolicSolver(eq1, eq2, eq3, eq4);
            case 'ODE'
                result = odeSolver(eq1, eq2, eq3, eq4);
        end

        % Display result
        result_textbox.Value = result;
    end
end

% Elliptic solver function
function result = ellipticSolver(eq1, eq2, eq3, eq4)
    % Laplace's equation solver using finite difference method
    % Set grid size
    nx = 50;
    ny = 50;
    % Set tolerance
    tol = 1e-6;
    % Initialize solution matrix
    u = zeros(nx, ny);
    % Set boundary conditions
    u(1,:) = 0;
    u(end,:) = 0;
    u(:,1) = 1;
    u(:,end) = 1;
    % Iterative solution using finite difference method
    while true
        u_new = u;
        for i = 2:nx-1
            for j = 2:ny-1
                u_new(i, j) = 0.25 * (u(i-1, j) + u(i+1, j) + u(i, j-1) + u(i, j+1));
            end
        end
        if max(abs(u_new - u), [], 'all') < tol
            break;
        end
        u = u_new;
    end
    result = u;
end

% Parabolic solver function
function result = parabolicSolver(eq1, eq2, eq3, eq4)
    % Heat equation solver using explicit finite difference method
    L = 1; % Length of the domain
    nx = 100; % Number of spatial grid points
    dx = L / (nx - 1); % Grid spacing
    nt = 100; % Number of time steps
    dt = 0.01; % Time step size
    alpha = 0.1; % Diffusion coefficient
    % Initialize solution array
    u = zeros(nx, nt + 1);
    % Set initial condition
    x = linspace(0, L, nx);
    u(:, 1) = exp(-(x - 0.5).^2 / (2 * 0.1^2));
    % Implement explicit finite difference method
    for t = 1:nt
        for i = 2:nx-1
            u(i, t+1) = u(i, t) + alpha * dt / dx^2 * (u(i+1, t) - 2 * u(i, t) + u(i-1, t));
        end
    end
    result = u;
end