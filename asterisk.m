function b = asterisk()
%adding asterisk to existing .dat files
%data_dir_str= ['/Users/polinavanyukov/Box Sync/Project Trust Game/regs/18-Jul-2016/'];
%data_dump_str = '/Users/polinavanyukov/Box Sync/Project Trust Game/regs/asterisked regs/';
%data_dir_str='/Users/polinavanyukov/Box Sync/Project Trust Game/regs/Sept 2018/';
%data_dump_str='/Users/polinavanyukov/Box Sync/Project Trust Game/regs/Sept 2018/asterisked_regs/';

%data_dir_str=['/Users/polinavanyukov/Regs/',date, '/'];
%data_dir_str='/Users/polinavanyukov/Regs/PEs_aligned_decision/';
data_dir_str=['/Users/polinavanyukov/Regs/clinical_policy_PEsValues_FixedPs/',date];
data_dump_str =['/Users/polinavanyukov/Regs/asterisked_regs/',date];

if ~exist(data_dump_str,'file')
    mkdir(data_dump_str)
    fprintf('Creating specific reg folder in: %s\n\n',data_dump_str);
end

cd(data_dir_str)
files = dir('*.dat');
num_of_subjects = length(files);

for index = 1:num_of_subjects
    filename=files(index).name;
    fprintf('File processing: %s\n', filename);
    x = load(filename);
    blocks = length(x)/48;
    ast = {'*', '*', '*'};
    if(blocks == 1)
        block1=num2cell(x(1:48,:));
        c = [block1];
    elseif(blocks == 2)
        block1=num2cell(x(1:48,:));
        block2=num2cell(x(49:96,:));
        c = [block1; ast; block2];
    elseif(blocks == 3)
        block1=num2cell(x(1:48,:));
        block2=num2cell(x(49:96,:));
        block3=num2cell(x(97:144,:));
        c = [block1; ast; block2; ast; block3];
    else
        block1=num2cell(x(1:48,:));
        block2=num2cell(x(49:96,:));
        block3=num2cell(x(97:144,:));
        block4=num2cell(x(145:192,:));
        c = [block1; ast; block2; ast; block3; ast; block4];
    end   
    dlmcell([data_dump_str filename '.dat'],c,'\t');
end


return

