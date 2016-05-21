function fraca(d, np, time, viz, rec)
if nargin < 5
    rec = 0;
end
if nargin < 4
    viz = 0;
end
close all
global dim  len order order2 maxdim
dim = d;len = np; 
order = nextpow2(dim)-1;
order2 = 2*order; maxdim = 2^order;
if rec, mov = avifile('FraCA.avi'); end
inpi = floor((dim-1)*rand(1,np))+1; %initialize 'particle' positions
inpj = floor((dim-1)*rand(1,np))+1;
mat = zeros(dim, dim); Sumat = mat;
[mat, inpi, inpj] = update(mat, inpi, inpj);
for k=1:time
    [mat, inpi, inpj] = update(mat, inpi, inpj);
    Sumat = Sumat + 1 - mat;
    %w = xcorr2(1-mat, 1-mat);
    %J(k) = entropy(1-mat);
    if viz,
        imagesc(1-mat), colormap(gray), 
        title(['Step = ',num2str(k)]), pause(0.1)
    end
    if rec
        F=getframe(gca);
        mov=addframe(mov,F);
    end
end
if rec, mov = close(mov);end
imagesc(Sumat*(100/time))
title(['Site occupation after ',num2str(time),' steps'])
%plot(J)
end

function [new, x, y] = update(m, x, y) %update each current 'particle' position
global len dim 
new = zeros(dim, dim);
for i=1:len
    for j=1:len
        [newx, newy] = BooleanMove(x(i), y(j));
        %[newx, newy] = VerletMove(x(i), y(j), i+maxdim*j);      
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
             x1 = x(ptr)-1; if x1<0, x1 = dim;end
             x2 = x(ptr)+1; if x2>dim, x2=1;end
             y1 = y(ptr)-1; if y1<0, y1 = dim;end
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

function [newx, newy] = BooleanMove(x, y) % update functions for background and 'particle' position via Boolean functions
global order
val = bitand(x, y)*bitxor(x, y);
newx = 1 - 2*bitget(val, floor(order*rand)+1);
newy = 1 - 2*bitget(val, floor(order*rand)+1);
end

function [newx, newy] = VerletMove(x, y, val) % update functions for background and 'particle' position a la Verlet
global order2 
sdig = sum(bitget( val, 1:order2));
newx = floor(10*sqrt(sdig)); newy = floor(10*sqrt(order2 - sdig));
end