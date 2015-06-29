function fxy = gauss_2D(nx,ny,sigx,sigy,angle)
%2d gauss, checks for even or odd length so as to center it correctly and
%avoid division by a 0 sigma.  0 sigma can be entered and it will return
%the correct gauss form, i.e repeated 1d gaussians

if isempty(sigy), sigy = sigx;end
if (nargin<5) || isempty(angle), angle = 0;end

fxy=zeros(nx,ny);


[x , y]=meshgrid( -(nx-1)/2:(nx-1)/2,-(ny-1)/2:(ny-1)/2);

if angle ~= 0
    theta=angle/180e0*pi;
    xd=x*cos(theta)-y*sin(theta);
    yd=x*sin(theta)+y*cos(theta);
    x=xd;
    y=yd;
    xd=0;
    yd=0;
end

if sigx < 0
    if mod(nx,2) == 1
        if sigx == 0, gx=(x == 0);else gx=exp(0.5.*x.^2./sigx^2);end
    else gx=exp(0.5.*x.^2./sigx^2);end
else
    if mod(nx,2) == 1
        if sigx == 0, gx=(x == 0);else gx=exp(-0.5.*x.^2./sigx^2);end
    else gx=exp(-0.5.*x.^2./sigx^2);end
end

if sigy < 0
    if mod(ny,2) == 1
        if sigy == 0, gy=(y == 0);else gy=exp(0.5.*y.^2./sigy^2);end
    else gy=exp(0.5.*y.^2./sigy^2);end
else
    if mod(ny,2) == 1
        if sigy == 0, gy=(y == 0);else gy=exp(-0.5.*y.^2./sigy^2);end
    else gy=exp(-0.5.*y.^2./sigy^2);end    
end

fxy=gx.*gy;

fxy=fxy/sum(sum(fxy));

end

