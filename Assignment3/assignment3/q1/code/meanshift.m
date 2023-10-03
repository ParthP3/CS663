function im2 = meanshift(im1, ss, sr, kernel_size,dist_weight, intens_weight, thresh)
%  MEANSHIFT Brief summary of this function.
% 
% Detailed explanation of this function.
a = size(im1);
im2 = zeros(a);
im1_pad = padarray(im1,[floor(kernel_size/2), floor(kernel_size/2)],0,'both');
mesh = zeros(1,kernel_size);
for i=1:size(mesh)
    mesh(1,i)=abs(i-floor(kernel_size/2));
end
[X,Y] = meshgrid(mesh,mesh);
sp = X.^2 + Y.^2;
Gss = (1/(ss*sqrt(2*pi)))*exp(-0.5*sp/ss^2);
mesh1 = zeros(1,kernel_size);
for i=1:size(mesh)
    mesh1(1,i)=i-ceil(kernel_size/2);
end
[X1,Y1] = meshgrid(mesh1);
cents = 1;
centroids = [0,0,0,0];
for i=floor(kernel_size/2)+1:a(1)+floor(kernel_size/2)
    for j=floor(kernel_size/2)+1:a(2)+floor(kernel_size/2)
        x2 = i;
        y2 = j;
        I2 = im1_pad(x2,y2);
        x1 = x2+2;
        y1 = y2+2;
        I1 = I2+5;
        while ~(abs(x2-x1)<=1 && abs(y2-y1)<=1 && abs(I1-I2)<=3)
            x1=x2;
            y1=y2;
            I1=I2;
            intensity = im1_pad(x1-floor(kernel_size/2):x1+floor(kernel_size/2),y1-floor(kernel_size/2):y1+floor(kernel_size/2));
            Gsr = (1/(sr*sqrt(2*pi)))*exp(-0.5*(I1-intensity).^2/sr^2);
            x = X1+x1;
            y = Y1+y1;
            x2=max( min( floor(sum(x.*Gss.*Gsr,'all')/sum(Gss.*Gsr,'all')), a(1)+floor(kernel_size/2) ), floor(kernel_size/2)+1);
            y2=max( min( floor(sum(y.*Gss.*Gsr,'all')/sum(Gss.*Gsr,'all')), a(2)+floor(kernel_size/2) ), floor(kernel_size/2)+1);
            I2=sum(intensity.*Gss.*Gsr,'all')/sum(Gss.*Gsr,'all');
        end
        check = 0;
        if cents == 1
            cents = cents+1;
            centroids(cents,:) = [x2, y2, I2, 1];
        else
            
            for p=2:cents
                if dist_weight*((centroids(p,1)-x2)^2 + (centroids(p,2)-y2)^2) + intens_weight*((centroids(p,3)-I2)^2) < thresh^2
                    I2 = centroids(p,3);
                    l = centroids(p,4);
                    centroids(p,1) = (l/(l+1))*(centroids(p,1)) + (1/(l+1))*(x2);
                    centroids(p,2) = (l/(l+1))*(centroids(p,2)) + (1/(l+1))*(y2);
                    centroids(p,3) = (l/(l+1))*(centroids(p,3)) + (1/(l+1))*(I2);
                    centroids(p,4) = l+1;
                    check = 1;
                    break;
                end
            end
        end
        if check == 0
            cents = cents+1;
            centroids(cents, :) = [x2, y2, I2, 1];
        end
        im2(i-floor(kernel_size/2),j-floor(kernel_size/2))=I2;
    end
end
cents
end