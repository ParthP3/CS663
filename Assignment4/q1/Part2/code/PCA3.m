function [output, output1] = PCA3(train, test, ks, count)
cs = size(count);
s0 = size(ks);
output =double(zeros(size(ks)));
output1 =double(zeros(size(ks)));
s1 = size(train);
s2 = size(test);
train_mean = mean(train, 2);
train1 = train-train_mean;
cov = (train1'*train1)/(s1(2)-1);
[V_1, D] = eig(cov);
[d, ind] = sort(diag(D), "descend");
V1 = V_1(:, ind);
V1 = train1*(V1);
norm = vecnorm(V1);
V1 = V1./norm;
for j=1:s0(2)
    k = ks(j);
    V = V1(:, 1:k);
    coeffs = V'*train1;
    test_coeffs = V'*(test-train_mean);

    if k>3
    test_coeffs1 = test_coeffs(4:k, :);
    end
    correct = 0;
    correct1 = 0;
    for i=1:s2(2)
        errors = coeffs - test_coeffs(:, i);
        if k>3
            errors1 = coeffs(4:k, :)- test_coeffs1(:, i);
            errors1 = errors1.^2;
            errors1 = mean(errors1, 1);
            [C1, I1] = min(errors1);
        end
        errors = errors.^2;
        errors = mean(errors, 1);
        
        [C, I] = min(errors);
        

        train_i = 0;
        person = 0;
        for p=1:cs(2)
            train_i = train_i+count(p);
            if(train_i >= i)
                break;
            end
            person = person+1;
        end
        
        if(floor((I(1)-1)/40) == person)
            correct = correct+1;
        end
        if k>3
        if(floor((I1(1)-1)/40) == person)
            correct1 = correct1+1;
        end
        end
    end
    if k<=3
        output(1,j) = (correct*100)/s2(2);
        output1(1, j) = (correct*100)/s2(2);
    else
        output(1,j) = (correct*100)/s2(2);
        output1(1, j) = (correct1*100)/s2(2);
    end

end
end