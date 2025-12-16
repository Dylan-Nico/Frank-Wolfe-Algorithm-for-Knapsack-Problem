
eps = 1e-8;

%% TEST 1 – Q ben condizionata + vincolo attivo
%% già fatto
disp("======== TEST 1 : Q ben condizionata + vincolo attivo ========");

Q = [1 0; 0 2];
q = [-4; -2];

a = [1; 1];
b = 2;

l = [0; 0];
u = [3; 3];

x0 = [2; 0];   % ammissibile

[x_FW, f_FW] = frank_wolfe(Q, q, x0, a, b, l, u, eps);

disp("Risultato TEST 1");
disp("x_FW = ");
disp(x_FW);
fprintf("a'*x = %.4f\n", a'*x_FW);
fprintf("f(x)  = %.6f\n\n", f_FW);



%% TEST 2 – Q ben condizionata + vincolo NON attivo

disp("======== TEST 2 : Q ben condizionata + vincolo NON attivo ========");

Q = [1 0; 0 2];
q = [-4; -2];

a = [1; 1];
b = 0.1;   % vincolo NON attivo

l = [0; 0];
u = [3; 3];

x0 = [1; 1];

[x_FW, f_FW] = frank_wolfe(Q, q, x0, a, b, l, u, eps);

disp("Risultato TEST 2");
disp("x_FW = ");
disp(x_FW);
fprintf("a'*x = %.4f (vincolo NON attivo)\n", a'*x_FW);
fprintf("f(x)  = %.6f\n\n", f_FW);



%% TEST 3 – Q mal condizionata
% dovremmo avere x1 bassa e x2 alta
disp("======== TEST 3 : Q mal condizionata ========");

Q = [20 0; 0 1];
q = [-4; -2];

a = [1; 1];
b = 2;

l = [0; 0];
u = [3; 3];

x0 = [2; 0];

[x_FW, f_FW] = frank_wolfe(Q, q, x0, a, b, l, u, eps);

disp("Risultato TEST 3");
disp("x_FW = ");
disp(x_FW);
fprintf("a'*x = %.4f\n", a'*x_FW);
fprintf("f(x)  = %.6f\n\n", f_FW);



%% TEST 4 – Q con autovalore ZERO
disp("======== TEST 4 : Q con autovalore ZERO========");

Q = [1 0; 0 0];
q = [-4; -2];

a = [1; 1];
b = 2;

l = [0; 0];
u = [3; 3];

x0 = [2; 0];

[x_FW, f_FW] = frank_wolfe(Q, q, x0, a, b, l, u, eps);

disp("Risultato TEST 4");
disp("x_FW = ");
disp(x_FW);
fprintf("a'*x = %.4f\n", a'*x_FW);
fprintf("f(x)  = %.6f\n\n", f_FW);

%% TEST 6 – Vincolo molto stringente (quasi tutto fuori)
% x1 e x2 al bordo del box
disp("======== TEST 6 : Vincolo molto stringente ========");

Q = [1 0; 0 1];
q = [-4; -2];

a = [1; 1];
b = 5.9;    % quasi massimo possibile

l = [0; 0];
u = [3; 3];

x0 = [3; 3];

[x_FW, f_FW] = frank_wolfe(Q, q, x0, a, b, l, u, eps);

disp("Risultato TEST 6");
disp("x_FW = ");
disp(x_FW);
fprintf("a'*x = %.4f (vincolo molto attivo)\n", a'*x_FW);
fprintf("f(x)  = %.6f\n\n", f_FW);


disp("============== FINE TEST ==============");