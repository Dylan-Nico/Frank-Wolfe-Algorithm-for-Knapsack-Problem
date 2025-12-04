function [Q, q, l, u, a, b, x0] = generate_variables(n)

 A = randn(n);
 Q = A' * A;
 q = randn(n,1); 
 LOW  = 1;
 HIGH =  10;

% genero i bound inferiori
 l = LOW + (HIGH-LOW) * rand(n,1);

% genero i bound superiori, sempre > l
 u = l + 0.5 + rand(n,1); 
 a = rand(n,1); 

 % calcolo range valido
 b_min = a' * l;
 b_max = a' * u;

 % scegli b in modo "random" dentro l'intervallo possibile
 b = b_min + (b_max - b_min) * rand();
 t = (b - a'*l) / (a'*(u-l));  
t = max(0, min(1, t));          
x0 = l + t*(u - l);             % x0 nel box, soddisfa il vincolo
 disp("n:")
 disp(n)
 disp("Q:")
 disp(Q)
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