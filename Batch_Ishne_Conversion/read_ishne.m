% Read ISHNE ecgSig file : 
%  function mat=read_ishne(fileName, startOffset, length);    
% Input--------------------------------------------------------       
%  fileName : the ishne filename including the path
%  startOffset: the start offset to read ecgSig
%  length:      length of ecgSig to read
% Ouput--------------------------------------------------------
%  ecgSig : ECG signal
%
function [ishneHeader, ecgSig]=read_ishne(fileName, startOffset, length, tmp_header)
fid=fopen(fileName,'r');
if ne(fid,-1)
    
    %Magic number
    magicNumber = fread(fid, 8, 'char');
   
    % get checksum
	checksum = fread(fid, 1, 'uint16');
	
	%read header
    Var_length_block_size = fread(fid, 1, 'int32');
    ishneHeader.Sample_Size_ECG = fread(fid, 1, 'long');	
    Offset_var_lenght_block = fread(fid, 1, 'long');
    ishneHeader.inf.Offset_var_lenght_block = Offset_var_lenght_block;
    Offset_ECG_block = fread(fid, 1, 'long');
    ishneHeader.inf.Offset_ECG_block=Offset_ECG_block;
    File_Version = fread(fid, 1, 'short');
    ishneHeader.inf.File_Version=File_Version;
    First_Name = fread(fid, 40, 'ubit8');
    ishneHeader.inf.First_Name=First_Name;
    Last_Name = fread(fid, 40, 'ubit8');
    ishneHeader.inf.Last_Name=Last_Name;
    ID = fread(fid, 20, 'ubit8');
    ishneHeader.inf.ID=ID;
    Sex = fread(fid, 1, 'short');
    ishneHeader.inf.Sex=Sex;
    Race = fread(fid, 1, 'short');
    ishneHeader.inf.Race=Race;
    Birth_Date = fread(fid, 3, 'short');
    ishneHeader.inf.Birth_Date=Birth_Date;
    Record_Date = fread(fid, 3, 'short');
    ishneHeader.inf.Record_Date=Record_Date;
    File_Date = fread(fid, 3, 'short');
    ishneHeader.inf.File_Date=File_Date;
    Start_Time = fread(fid, 3, 'short');
    ishneHeader.inf.Start_Time=Start_Time;
    ishneHeader.nbLeads = fread(fid, 1, 'short');
    Lead_Spec = fread(fid, 12, 'short');
    ishneHeader.inf.Lead_Spec=Lead_Spec;
    Lead_Qual = fread(fid, 12, 'short');
    ishneHeader.inf.Lead_Qual=Lead_Qual;
    ishneHeader.Resolution = fread(fid, 12, 'short');	
    Pacemaker = fread(fid, 1, 'short');
    ishneHeader.inf.Pacemaker=Pacemaker;
    Recorder = fread(fid, 40, 'char');
    ishneHeader.inf.Recorder=Recorder;
    ishneHeader.Sampling_Rate = fread(fid, 1, 'short');	
    Proprietary = fread(fid, 80, 'ubit8');
    Copyright = fread(fid, 80, 'ubit8');
    Reserved = fread(fid, 88, 'ubit8');
    
    % read variable_length block
    varblock = fread(fid, Var_length_block_size, 'ubit8');
    
    % get data at start
    offset = startOffset*ishneHeader.Sampling_Rate*ishneHeader.nbLeads*2; % each data has 2 bytes
    fseek(fid, Offset_ECG_block+offset, 'bof');
    
   
    % read ecgSig signal
    numSample = ishneHeader.Sample_Size_ECG*ishneHeader.Sampling_Rate;
%     numSample = 100000;% length*ishneHeader.Sampling_Rate;
    ecgSig = fread(fid, [ishneHeader.nbLeads, numSample], 'int16')';
    
    fclose(fid);
 else
     ihsneHeader = [];
     ecgSig=[];
     disp('Error:  Cannot read file');
 end
