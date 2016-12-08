%collinearity
 files1 = dir('*prosocial.dat');
 num_of_subjects = length(files1);
 files2 = dir('*consistWrep1*.dat');
 i = 1;
% 
for index = 1:num_of_subjects
    if index~=36
        filename1 = files1(index).name;
        fprintf('File processing: %s\n', filename1);
        id = filename1(isstrprop(filename1,'digit'));
        regs1 = load(filename1);
        regs2 = load(files2(index).name);
        [R, P] = corrcoef(regs1(:,3),regs2(:,3));
        corr_matrix(i) = R(1,2);
        i=i+1;
    end
end