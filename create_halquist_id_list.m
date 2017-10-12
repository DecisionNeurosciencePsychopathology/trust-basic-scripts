function id_list=create_halquist_id_list(path_to_files)
%So far all ids are put into folders with id being the start of the string,
%extract that and make it a numerical list to be used to processing into
%.mat files

%Grab the files into a cell
file_list=dir(path_to_files);
file_list = {file_list.name}';

%Create regexp and pull numerical ids from file string (do not grab leading
%zeros!)
expression = '([1-9]\d+)';
ids = regexp(file_list,expression,'match','once'); %Just the first match

%Convert to double and remove zeros, there must be a vectorized way!!!
for i = 1:length(ids)
    try
        id_list(i) = str2num(cell2mat(ids(i))); 
    catch
        continue
    end
end
id_list(id_list==0)=[];
