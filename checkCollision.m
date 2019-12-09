function doesCollide = checkCollision(xi, xj, yi, yj, zi, zj, surfEq)

syms x y z d;
doesCollide = false;

x = xi + (xj - xi) * d;
y = yi + (yj - yi) * d;
z = zi + (zj - zi) * d;

surfEq = eval(surfEq);

d = double(vpasolve(surfEq));

% empty means complex/imaginary, no collision forever.
if isempty(d)
    return
end

x = double(eval(x));
y = double(eval(y));
z = double(eval(z));

if x >= min(xi,xj) && x <= max(xi,xj) && y >= min(yi,yj) && y <= max(yi,yj) && z >= min(zi,zj) && z <= max(zi,zj)
    doesCollide = true;
end