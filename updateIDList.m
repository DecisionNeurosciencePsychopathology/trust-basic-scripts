function updateIDList( varargin, spss_flag )
%This funciton assumes you have access to both the W: and R: (I think)
%drives. Since this function is specific for only two files maybe it would
%be ok to hard code?

%5/13/15 So now this code will run a macro and save the query spalshdemo2
%from my remote copy in the Y:drive. It also saves a copy to the L drive as
%a backup, and creates a local backup (thikning about that now it seems
%like it's overkill?). this code will also create and save that query in a
%.mat file which will be used for the subsequent Neuropsych processing
%functions (CreateSubjList, LoadAllIds, ect). This function will now also
%update the spss demogfile as well in one go

%Example of how to run
%UpdateIDList <- this will check to see if the list needs updated and react
%accordingly

%UpdateIDList([],1) <- this will check to see if the list needs updated and react
%accordingly, then update the spss .dat file

%UpdateIDList('force',1) <- force updates list, updates spss demog .dat file

%%%%%%%%-------BEGIN CODE-------------%%%%%%%%%%%%%%


%Added this so we could overwrite if needed
if(strcmpi('force',varargin))
    run_flag = true;
else
    run_flag = false;
end

%Update spss demog file
if nargin < 2
    spss_flag=0;
end


data_dir = '/Users/polinavanyukov/Documents/tmp/';

orig=pwd;
%fpath = [pathroot 'db/master id list.xlsx'];
%fpath = 'L:/Summary Notes/Data/matlab/db/master id list.xlsx';
fpath = '/Users/polinavanyukov/Box Sync/Project Trust Game/data/splashDemo2.xlsx';

if ~exist(fpath,'file')
    error('No master list found, see help'); % Just in case
end

hardcopy = '/Users/polinavanyukov/Desktop/rcn/cpad/PROTECT/Front End (For Copying)/Jon/Protect 2 - Front End (Jon).accdb';
%Check original file's date
org_file=dir(strcat(fpath));
org_file=org_file.date;

%Check databases file date
db_file = dir(strcat(hardcopy));
db_file=db_file.date; 

%Compare dates buffer added to negate redundancy
if datenum(org_file)+.03<=datenum(db_file) || run_flag
    disp('Database is out of date! Don''t worry we can put a bird on it....');
    
    %Connect to database and export mast list
    h= actxserver('Access.Application');
    invoke(h,'OpenCurrentDatabase',hardcopy); %Protect DB
    invoke(h.DoCmd,'RunMacro','ExportDemographics'); %Macro is currently only executable on my (Jon's) computer
    h.Visible = 0;
    
    %Copy it to the L drive
    %cd('L:/Summary Notes/Data/matlab/db/')
    copyfile(fpath,'/Users/polinavanyukov/Desktop/llmd/Summary Notes/Data/matlab/db/splashDemo2.xlsx');
%     movefile('master id list.xlsx','master id list_old_bckup.xlsx'); %Create backup
%     movefile('master_id_list.xlsx','master id list.xlsx');
    disp('All DONE, press any key to continue!');
    %pause();
else
    disp('Database is up to date!');
end

%This portion of the code will replace the demog_data file in the tmp
%directory used for the older Jan code

%Load data into matlab
filename = '/Users/polinavanyukov/Desktop/llmd/Summary Notes/Data/matlab/db/splashDemo2.xlsx';
[~,~,data] = xlsread(filename);

%Strip header information
header_info = data(1,:);
data(1,:)=[];

%Map out needed info from cols index, rmv white space for setfield
index=[];
for i = 1:length(header_info)
    %index = setfield(index, regexprep(header_info{i},'[^\w'']',''),...
    %find(strcmpi(header_info,header_info(i))));
    index = setfield(index, regexprep(header_info{i},'[^\w'']',''),i);
end


%Grab all ids numbers
id_nums=[data{:,index.ID}]';
tmp = zeros(length(id_nums),1);

%Find where there are duplicates
for i = 2:length(id_nums)
    if id_nums(i)==id_nums(i-1)
        tmp(i,1)=1;
    end
end
idx=find(tmp);
count = 0;
for i=1:length(idx)
    %Not so clean workaround look into making this more robust
    idx(i)=idx(i)-count;
    
    %PROTO field
    data{idx(i)-1,index.PROTO}=strcat(data{idx(i)-1,index.PROTO},'/',data{idx(i),index.PROTO});
    
    %Combine dates adjust for SPSS UTF-16 error
    if ~strcmpi(data{idx(i)-1,index.Consent},data{idx(i),index.Consent}) && strcmpi(data{idx(i)-1,index.Consent},'NaN')
        data{idx(i)-1,index.Consent}=strcat(data{idx(i)-1,index.Consent},'--',data{idx(i),index.Consent});
    end
    
    if ~strcmpi(data{idx(i)-1,index.DateofTermination},data{idx(i),index.DateofTermination}) && strcmpi(data{idx(i)-1,index.DateofTermination},'NaN')
        data{idx(i)-1,index.DateofTermination}=strcat(data{idx(i)-1,index.DateofTermination},'--',data{idx(i),index.DateofTermination});
    end
    
    %Delete the duplicate row
    data(idx(i),:)=[];
    
    count=count+1;
end

%Clean up "other paitent" data, rmv 'NaN' entries
other_idx=strcmp(data(:,index.PATTYPE),'OTHER PATIENT');
data(other_idx,:)={'NaN'};
qnan = cellfun(@any,cellfun(@isnan,data,'UniformOutput',0));
qnan = ( qnan |  strcmp('NaN',data) ); % get rid of 'NaN' entries
data(qnan) = {''}; % replace NaN's with zeros
emptyCells = cellfun('isempty', data);
data(all(emptyCells,2),:) = [];

% Archive and save it
cd(data_dir)
copyfile('demogs_data.mat', 'demogs_data_archive.mat');
save([data_dir 'demogs_data'],'data');


if spss_flag ==1
    % open file pointer for writing DAT file
    cd(pathroot)
    cd('../SPSS/data')
    fid = fopen('demogspss.dat' ,'w');
    
    
    %Take care of headers
    for hn = 1:length(header_info)
        fprintf(fid,'%s\t',header_info{hn});
    end
    fprintf(fid, '\n');
    
    str_fmt=[]; %automate str formatting
    for i = 1:size(data,2)
        if ischar(data{1,i})
            tempstr = '%s';
        else
            tempstr = '%d';
        end
        str_fmt = strcat(str_fmt,tempstr,'\t');
    end
    str_fmt(end)='n';
    
    for i=1:size(data,1)
        fprintf(fid,str_fmt,data{i,:});
    end
    
    % Save it
    save([data_dir 'demogs_data'],'data');
end

%Just added this is so it's one less function to run really
createSubjIDlist;

% varargout
if(nargout), demogs = data; end

%kill pointer
fclose all;

%Return to origin
cd(orig);
end

