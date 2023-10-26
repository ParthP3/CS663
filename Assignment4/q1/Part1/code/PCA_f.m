function [output, coeff] = PCA_f(train, test, outside, ks, type)

s0 = size(ks);
output =double(zeros(size(ks)));
s1 = size(train);
s2 = size(test);
s3 = size(outside);
train_mean = mean(train, 2);
train1 = train-train_mean;

if type == 2
    cov = (train1'*train1)/(s1(2)-1);
    [V_1, D] = eig(cov);
    [d, ind] = sort(diag(D), "descend");
    V1 = V_1(:, ind);
    V1 = train1*(V1);
    norm = vecnorm(V1);
    V1 = V1./norm;
    
elseif type == 1
    svd1 = train1'*train1;
    [U, S, U1] = svd(svd1);
    V1= train1*(U);
    norm = vecnorm(V1);
    V1 = V1./norm;

elseif type == 0
    cov = (train1*train1')/(s1(2)-1);
    [V_1, D] = eig(cov);
    [d, ind] = sort(diag(D), "descend");
    V1 = V_1(:, ind);
end
figure;
recon = 89;
imshow(reshape((train(:,recon))/255, [92,112])');
t = [2, 10, 20, 50, 75, 100, 125, 150, 175];

figure;
for t_c =1:9
    V = V1(:, 1:t(1, t_c));
    coeffs = V'*train1;
    subplot(3,3,t_c);
    img = (reshape(train_mean + V*(coeffs(:,recon)), [92, 112]))';
    %min1 = min(img, [], 'all');
    max1 = max(img, [], 'all');
    img = (1/(max1))*(img);
    imshow(img);        
end

c1 = [];
c1_index = 1;
c2 = [];
c2_index = 1;
figure;
fp = 0;
fn = 0;
e = 4.5*10^4;
for j=1:s0(2)
    k = ks(j);
    V = V1(:, 1:k);
    %size = dxk
    coeffs = V'*train1;
    test_coeffs = V'*(test-train_mean);
    out_coeffs = V'*(outside-train_mean);
    correct = 0;
    if k == 170
        for i = 1:s3(2)
            errors1 = coeffs - out_coeffs(:, i);
            errors1 = errors1.^2;
            errors1 = mean(errors1, 1);
            [C1, I1] = min(errors1);
            if (C1(1) < e)
                fp = fp+1;
            end
            c2(c2_index) = C1(1);
            c2_index = c2_index+1;
        end
    end
    for i=1:s2(2)
        errors = coeffs - test_coeffs(:, i);
        errors = errors.^2;
        errors = mean(errors, 1);
        [C, I] = min(errors);
        if k ==170
            if C(1)>e
                fn = fn+1;
            end
            c1(c1_index) = C(1);
            c1_index = c1_index+1;
        end
        if(floor((I(1)-1)/6) == floor((i-1)/4))
            correct = correct+1;
        end
    end
    output(1,j) = (correct*100)/s2(2);
end
fp
fn
coeff = [c1,c2];
figure;
for i = 1:25
    subplot(5,5,i);
    img = reshape(V1(:,i), [92,112])';
    min1 = min(img, [], 'all');
    max1 = max(img, [], 'all');
    img = (1/(max1-min1))*(img-min1);
    imshow(img);
end
end