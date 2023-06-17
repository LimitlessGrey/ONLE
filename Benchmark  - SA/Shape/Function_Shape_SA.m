function mass_shape = Function_Shape_SA(x)

P = 10^5; % Load in N
A = (40 * 10^-3)^2 / 4 * pi; % Bars Area [m^3]
E = 68.9 * 10^9; % Young modulus in Pa
rho = 2700; % [kg/m^3]

conn = readmatrix('Connectivity.txt');
coord = [0 0.915; x(1) x(2); x(3) x(4); 0 0; x(5) x(6); 0.915*2 0];

% get the number of rows of connectivity matrix
num_rows_conn = size(conn,1);

% get the number of rows of coordinates matrix
num_rows_coord = size(coord,1);

L = zeros(num_rows_conn, 1); % Initialize bar lenght vector

k_g = zeros(2*num_rows_coord); % initialize global stiffness matrix

% Global stiffness matrix calculation
for i=1:num_rows_conn
    n1 = conn(i, 1);
    n2 = conn(i, 2);
    x_n1 = coord(n1, 1);
    y_n1 = coord(n1, 2);
    x_n2 = coord(n2, 1);
    y_n2 = coord(n2, 2);

    L(i, 1) = norm([x_n2, y_n2] - [x_n1, y_n1]);

    alfa = atan2d(y_n2 - y_n1, x_n2 - x_n1);

    k_e = E * A / L(i, 1) * [cosd(alfa)^2, cosd(alfa)*sind(alfa), -cosd(alfa)^2, -cosd(alfa)*sind(alfa); ...
        cosd(alfa)*sind(alfa), sind(alfa)^2, -cosd(alfa)*sind(alfa), -sind(alfa)^2; ...
        -cosd(alfa)^2, -cosd(alfa)*sind(alfa), cosd(alfa)^2, cosd(alfa)*sind(alfa); ...
        -cosd(alfa)*sind(alfa), -sind(alfa)^2, cosd(alfa)*sind(alfa), sind(alfa)^2];

    % Divide element stiffness matrix into 4 sections
    sections = mat2cell(k_e, [2 2], [2 2]);

    % Define the positions in the global stiffness matrix to add the element stiffness matrix
    positions = [n1, n1; n1, n2; n2, n1; n2, n2];

    % Iterate over the positions and sections and add the sections to the positions
    for j = 1:size(positions, 1)
        row_idx = positions(j, 1); 
        col_idx = positions(j, 2);
        section = sections{j};
        k_g(row_idx*2-1:row_idx*2, col_idx*2-1:col_idx*2) = k_g(row_idx*2-1:row_idx*2, col_idx*2-1:col_idx*2) + section;
    end
end

% Define the applied load vector
F_applied = zeros(num_rows_coord*2, 1);
F_applied(10) = -P; % Node 5, y-component
F_applied(12) = -P; % Node 6, y-component

% Initialize displacement vector
u = zeros(num_rows_coord*2, 1);

% Define the indices of the fixed nodes (1, 4 and 6)
fixed_nodes = [1*2-1 1*2 4*2-1 4*2 6*2-1 6*2];

% Define the indices of the free nodes
free_nodes = setdiff(1:size(k_g, 1), fixed_nodes);

% Solve for the displacement vector
u(free_nodes, :) = k_g(free_nodes, free_nodes) \ (F_applied(free_nodes, :) - k_g(free_nodes, fixed_nodes) * u(fixed_nodes, :));

% Reshape displacement vector into a 6x2 matrix
u_reshaped = reshape(u, 2, 6)';
u_reshaped = reshape(u_reshaped, 6, 2);

% Add the displacements to the coordinate matrix
new_coord = coord + u_reshaped;

% Initialize final bar lenght vector
new_L = zeros(num_rows_conn, 1);
% Initialize delta L vector
dL = zeros(num_rows_conn, 1);
% Initialize axial bar force vector
f = zeros(num_rows_conn, 1);
% Initialize axial stress vector
sigma = zeros(num_rows_conn, 1);

for i=1:num_rows_conn
    n1 = conn(i, 1);
    n2 = conn(i, 2);
    x_n1 = new_coord(n1, 1);
    y_n1 = new_coord(n1, 2);
    x_n2 = new_coord(n2, 1);
    y_n2 = new_coord(n2, 2);

    new_L(i, 1) = norm([x_n2, y_n2] - [x_n1, y_n1]);
    dL(i, 1) = new_L(i, 1) - L(i, 1);

    f(i, 1) = E * A / L(i ,1) * dL(i, 1); % Axial loads on bars
    sigma(i, 1) = f(i, 1) ./ A; % Axial stresses on bars
end

mass_shape = 0;
r = 10;
for i=1:num_rows_conn
    if abs(sigma(i, 1)) >= 276e6 || new_L(i, 1) < 0.2
        g1 = abs(sigma(i, 1)) - 276e6;
        g2 = -new_L(i, 1) + 0.2;
        mass_shape = mass_shape + rho * new_L(i, 1) * A + r*max(0, g1/276e6)^2 + r*max(0, g2/0.2)^2;
    else
        mass_shape = mass_shape + rho * new_L(i, 1) * A;
    end
end
end