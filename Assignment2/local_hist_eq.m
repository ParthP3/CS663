function out_img = local_hist_eq(img, kernel_size)

bins = 64;
div = 255.0/bins;
a = size(img);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
img_padded = padarray(img, [floor(kernel_size/2), floor(kernel_size/2)], 0, "both");

out_img = double(zeros(a));
K=4;
for i=1:a(1)
    for j=1:a(2)
        histogram = zeros(1, bins);
        intens = img_padded(i:i+kernel_size-1, j:j+kernel_size-1);
        % Constructing hisogram
        % =======================================================================================================================
        for p=1:kernel_size
            for q=1:kernel_size
                histogram(1, min(max(ceil(intens(p,q)/div), 1), bins)) = histogram(1, min(max(ceil(intens(p,q)/div), 1), bins)) + 1;
            end
        end
        %=========================================================================================================================
        histogram = histogram/(kernel_size^2);
        curr_int = intens(ceil(kernel_size/2), ceil(kernel_size/2));
        curr_bin = min(max(ceil(curr_int/div), 1), bins);
        final_int = 0.0;
        % Finding I'(i,j)
        for k=1:curr_bin
            if k == curr_bin
                final_int = final_int+(histogram(1,k)/K);
            else
                final_int = final_int + histogram(1, k);
            end
        end
        final_int = floor(255*final_int);
        out_img(i, j) = final_int;
    end
end
end