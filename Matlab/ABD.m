% Material properties

E1  = 11380;
E2  = 161000;
E3  = 11380;     
nu21 = 0.32;
nu13 = 0.45;     
nu23 = 0.32;
nu12 = (E1/E2)*nu21;

G12 = 5200;
G13 = 3900;     
G23 = 5200;

% Laminate definition
theta_deg = [60, -30, -60, 30, 60, -30, -60, 30];
tply = 0.15;
nply = numel(theta_deg);

% Reduced stiffness Q

nu21_calc = nu12 * E2 / E1;
den = 1 - nu12*nu21_calc;

Q11 = E1 / den;
Q22 = E2 / den;
Q12 = nu12 * E2 / den;
Q66 = G12;

Q = [ Q11  Q12  0;
      Q12  Q22  0;
      0    0    Q66 ];

% z-coordinates through thickness

h = nply * tply;
z = linspace(-h/2, h/2, nply+1);  % ply interfaces (bottom -> top)

% Build A, B, D

A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);

for k = 1:nply
    th = theta_deg(k);
    Qbar = Qbar_from_Q(Q, th);

    zk1 = z(k);
    zk2 = z(k+1);

    A = A + Qbar * (zk2 - zk1);
    B = B + 0.5 * Qbar * (zk2^2 - zk1^2);
    D = D + (1/3) * Qbar * (zk2^3 - zk1^3);
end

ABD = [A B; B D];


% Print results
disp('A (units: modulus*length):'); disp(A);
disp('B (units: modulus*length^2):'); disp(B);
disp('D (units: modulus*length^3):'); disp(D);
disp('ABD:'); disp(ABD);


function Qbar = Qbar_from_Q(Q, theta_deg)

    m = cosd(theta_deg);
    n = sind(theta_deg);

    Q11 = Q(1,1); Q22 = Q(2,2); Q12 = Q(1,2); Q66 = Q(3,3);

    Qbar11 = Q11*m^4 + 2*(Q12+2*Q66)*m^2*n^2 + Q22*n^4;
    Qbar22 = Q11*n^4 + 2*(Q12+2*Q66)*m^2*n^2 + Q22*m^4;
    Qbar12 = (Q11+Q22-4*Q66)*m^2*n^2 + Q12*(m^4+n^4);
    Qbar16 = (Q11 - Q12 - 2*Q66)*m^3*n - (Q22 - Q12 - 2*Q66)*m*n^3;
    Qbar26 = (Q11 - Q12 - 2*Q66)*m*n^3 - (Q22 - Q12 - 2*Q66)*m^3*n;
    Qbar66 = (Q11+Q22-2*Q12-2*Q66)*m^2*n^2 + Q66*(m^4+n^4);

    Qbar = [Qbar11 Qbar12 Qbar16;
            Qbar12 Qbar22 Qbar26;
            Qbar16 Qbar26 Qbar66];
end
