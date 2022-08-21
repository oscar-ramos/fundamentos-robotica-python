% Fundamentos de Robotica
% Ejemplo: Metodo de Newton en 1D
% ================================
% Encontrar la interseccion con cero de f(x)=3+(x-2).^3

% Valor inicial (cualquier valor)
x0 = 0;
% Iteraciones del metodo de Newton
for i=1:10
   x0 = x0 - (3+(x0-2)^3)/(3*(x0-2)^2);
end
disp(['Interseccion en: ' num2str(x0)])

% Grafico de la funcion
x = 0:0.001:1;
y = 3 + (x-2).^3;
plot(x,y), hold on, plot(x0,0,'o'), grid on
title('f(x)=3+(x-2)^3')
