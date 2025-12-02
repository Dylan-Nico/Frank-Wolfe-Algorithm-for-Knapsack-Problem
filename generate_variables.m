function [M, q, l, u, a, b, x0] = generate_variables(n)

 A = randn(n);
 M = A' * A;
 q = randn(n,1); %sempre valido
 LOW  = 1;
 HIGH =  10;

% generare i bound inferiori
 l = LOW + (HIGH-LOW) * rand(n,1);

% generare i bound superiori, sempre > l
 u = l + 0.5 + rand(n,1);   % così u_i > l_i garantito
 a = rand(n,1);  % tutti positivi → più stabile

 % calcola range valido
 b_min = a' * l;
 b_max = a' * u;

 % scegli b in modo "random" dentro l'intervallo possibile
 b = b_min + (b_max - b_min) * rand();
 t = (b - a'*l) / (a'*(u-l));   % scala tra 0 e 1
t = max(0, min(1, t));          % assicura 0 <= t <= 1
x0 = l + t*(u - l);             % x0 nel box, soddisfa il vincolo
 disp("n:")
 disp(n)
 disp("M:")
 disp(M)
 disp("q:")
 disp(q)
 disp("l:")
 disp(l)
 disp("u:")
 disp(u)
 disp("a:")
 disp(a)
 disp("b:")
 disp(b)
 disp("x0:")
 disp(x0)

end