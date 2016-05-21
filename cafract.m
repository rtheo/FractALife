function cafract(ord, np, time)
close all
global order len dim
order = ord;len = np;
dim = 2^order;  %construct an arbitrary arithmetic (number theoretic) fractal as the spatial background
for i=1:dim
    for j=1:dim
        matand(i, j) = bitand(i, j);
        matxor(i, j) = bitxor(i, j);
        mat(i, j) = abs(matand(i, j) - matxor(i, j));
    end
end
mov = avifile('FractalCA.avi');
inpi = floor((dim-1)*rand(1,np))+1; %initialize 'particle' positions
inpj = floor((dim-1)*rand(1,np))+1;
[inpi, inpj] = collision(inpi, inpj);
[newmat, inpi, inpj] = update(mat, inpi, inpj);
imagesc(newmat)
for k=1:time
    [newmat, inpi, inpj] = update(mat, inpi, inpj);
    imagesc(1 - newmat), colormap(gray), pause(0.1)
     F=getframe(gca);
     mov=addframe(mov,F);
end
close(mov);
end

function [new, x, y] = update(m, x, y) %update each current 'particle' position
global len dim
new = zeros(dim, dim);
for i=1:len
    for j=1:len
        [new(x(i), y(i)), newx, newy] = move(m(x(i), y(i)), x(i), y(i));
        x(i) = x(i) + newx; y(i) = y(i) + newy;
        x(i) = mod(x(i), dim+1); y(i) = mod(y(i), dim+1);   % boundary conditions
        if x(i)==0, x(i)=1; end, if y(i)==0, y(i) = 1; end
    end
end
[x, y] = collision(x, y);
for i=1:len
    for j=1:len
        new(x(i), y(i)) = 1;
    end
end
end

function [x, y] = collision(x, y)
global len dim
% Collision Rules
z = x + dim*y;
u = unique(z);
if length(u) < len, 
    for i=1:length(u)
         v = find( u(i) == z ); % find indices of unique elements
         for j=2:length(v)
             ptr = v(j);    % check if neighbor positions are free 
             x1 = x(ptr)-1; if x1<0 x1 = dim;end
             x2 = x(ptr)+1; if x2>dim, x2=1;end
             y1 = y(ptr)-1; if y1<0 y1 = dim;end
             y2 = y(ptr)+1; if y2>dim, y2=1;end            
             idx1 = length(find( x(ptr) -1== x)); 
             idx2 = length(find( x(ptr)+1== x)); 
             idx3 = length(find( y(ptr) -1 == y));
             idx4 = length(find( y(ptr)+1 == y));
             if (idx1 == 0), x(ptr) = x1;else
             if (idx2 == 0), x(ptr) = x2; end, end
             if (idx3 == 0), y(ptr) = y1;else
             if (idx4 == 0), y(ptr) = y2; end, end             
         end
    end
end
end

function [val, newx, newy] = move(m, x, y) % update functions for background and 'particle' position 
global order 
val = bitand( m, x)*bitxor(m, y);
newx = 1 - 2*bitget(val, floor(order*rand)+1);
newy = 1 - 2*bitget(val, floor(order*rand)+1);
val = 1;
end