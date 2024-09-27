% Prompts for batch conversion of ecgs, ishne to arash formats

% Prompt user for ecg folder and automatically create a folder to store outputs
inp_folder = uigetdir('.','Please select input folder...'); 
out_folder = strcat(inp_folder,'\ArashFormats');
mkdir(out_folder);

% Extract a list of .ecg files in input folder
dirList = dir(strcat(inp_folder,'\*.ecg'));

% Iterate through ecg file list and perform conversion
for i = 1:length(dirList)
    
    [~,name,~] = fileparts(dirList(i).name);
    inp_namestr  = [inp_folder '/' name '.ecg'];
    disp(inp_namestr);
    out_namestr  = [out_folder '/' name '.mat'];

    % V This conversion function can also take optional arguments scale_factor and first_leadind
    convert_ishne_2_arashfmt(inp_namestr,out_namestr); 
end
