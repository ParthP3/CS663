function out_img = mybilateralfilter(img, sigma_s, sigma_r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
out_img = double(zeros(size(img)));
kernel_size = 11;
mesh = zeros(1, kernel_size);
for i=1:kernel_size
    mesh(1,i) = abs(i-(ceil(kernel_size/2.0)));
end
img_padded = padarray(img, [floor(kernel_size/2.0), floor(kernel_size/2.0)], 0, "both");
[X, Y] = meshgrid(mesh, mesh);
dist = X.^2+Y.^2;
gauss_dist = exp(-dist/(2*pi*(sigma_s^2)))/(2*pi*sigma_s);
a = size(img);
for i=1:a(1)
    for j=1:a(2)
        intens = img_padded(i:i+kernel_size-1, j:j+kernel_size-1);
        intens1 = intens(ceil(kernel_size/2.0), ceil(kernel_size/2.0))-intens;
        gauss_intens = exp(-(intens1.^2)/(2*pi*(sigma_r^2)))/(2*pi*sigma_r);
        final = gauss_dist.*gauss_intens.*intens;
        W = sum(gauss_dist.*gauss_intens, "all");

        final = sum(final, "all");
        out_img(i, j) = final/W;
    end
end
end