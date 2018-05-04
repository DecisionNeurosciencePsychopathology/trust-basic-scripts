function b = asterisk()
%adding asterisk to existing .dat files
data_dir_str= 'E:/data/trust/regs/';
data_dump_str = 'E:/data/trust/regs/asterisked_regs/';

if ~exist(data_dump_str,'file')
    mkdir(data_dump_str)
    fprintf('Creating specific reg folder in: %s\n\n',data_dump_str);
end

cd(data_dir_str)
files = dir('*.dat');
num_of_subjects = length(files);

parfor index = 1:num_of_subjects
    try
        filename=files(index).name;
        fprintf('File processing: %s\n', filename);
        x = load(filename);
        block1=num2cell(x(1:48,:));
        block2=num2cell(x(49:96,:));
        block3=num2cell(x(97:144,:));
        block4=num2cell(x(145:192,:));
        ast = {'*', '*', '*'};
        c = [block1; ast; block2; ast; block3; ast; block4];
        %writetable(cell2table(c), [data_dump_str filename],'Delimiter','\t');
        dlmcell([data_dump_str filename],c,'\t');
    catch
        continue
    end
end

%Move all censor files into final dir
%N.B. you may jsut want to fix the rsync funtion in the custom move reg exp
%file to dump the astrisked regs in the "regs" folder of trust_analyses,
%that way to don't have to worry about new paths, but if you end up doing
%this, you must change the 'make_baseline_model' script to refelect the new
%path change!
copyfile([data_dir_str '*censor*'], data_dump_str);

return

