% Fundamentos de Robotica
% Ejemplo: Metodo del gradiente para un robot RR
% ================================================

xd = [1.2; 1.2];    % Valor deseado en el espacio Cartesiano
q  = [0; 0];        % Valor Inicial en el espacio articular
epsilon = 1e-3;
max_iter = 1000;    % Maximo numero de iteraciones
alpha = 0.5;

% Iteraciones del metodo del gradiente
iter = 0;
for i=1:max_iter
    q1=q(1); q2=q(2);
    JT = [-sin(q1)-sin(q1+q2), cos(q1)+cos(q1+q2);
                  -sin(q1+q2),         cos(q1+q2)];
    f = [cos(q1)+cos(q1+q2); sin(q1)+sin(q1+q2)];
    e = xd-f;
    q = q + alpha*JT*e;
    iter = iter+1;
    % Condicion de termino
    if (norm(e)<epsilon)
        break;
    end
end
disp(q)

l1=1; l2=1;
rr_plot(q,l1,l2); axis([0,1.4,0,1.4])
