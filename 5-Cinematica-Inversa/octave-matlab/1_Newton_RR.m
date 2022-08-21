% Fundamentos de Robotica
% Ejemplo: Metodo de Newton para un robot RR
% ============================================

xd = [1.2; 1.2];    % Valor deseado en el espacio Cartesiano
q  = [0.5; 0.5];    % Valor Inicial en el espacio articular
epsilon = 1e-3;
max_iter = 100;     % Maximo numero de iteraciones

% Iteraciones del metodo de Newton
iter = 0;
for i=1:max_iter
    q1=q(1); q2=q(2);
    Jinv = 1/sin(q2)*[cos(q1+q2), sin(q1+q2);
                      -cos(q1)-cos(q1+q2), -sin(q1)-sin(q1+q2)];
    f = [cos(q1)+cos(q1+q2); sin(q1)+sin(q1+q2)];
    e = xd-f;
    q = q + Jinv*e;
    iter = iter+1;
    % Condicion de termino
    if (norm(e)<epsilon)
        break;
    end
end
disp(q)

l1=1; l2=1;
rr_plot(q,l1,l2); axis([0,1.4,0,1.4])
