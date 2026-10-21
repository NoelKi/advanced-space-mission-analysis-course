matrix_with_zeros = zeros(10000,6);
matrix_with_ones = ones(10000,6);
matrix_with_nan = nan(10000,6); 


tA = tic;
for n = 3:35;
    disp(n)
end
measuaredTimeA = toc(tA)

