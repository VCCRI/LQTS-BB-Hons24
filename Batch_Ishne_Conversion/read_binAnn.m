% Read ISHNE binary ann file : 
%  function [ihsneHeader, Ann, Rloc]=read_ishne(fileName);
% Input--------------------------------------------------------       
%  fileName : the ishne binary annotation filename including the path
% Ouput--------------------------------------------------------
%  ishneHeader : ishne header info
%  Ann         : annotation of the beats
%  Rloc          : R peak location
%

function [ihsneHeader, Ann, Rloc, RR]=read_ishne(fileName, endPos);

fid=fopen(fileName,'r');
if ne(fid,-1)
   
    %Magic number
    magicNumber = fread(fid, 8, 'uchar');
   
    % get checksum
	checksum = fread(fid, 1, 'uint16');
	
	%read header
    Var_length_block_size = fread(fid, 1, 'long');
    ihsneHeader.Sample_Size_ECG = fread(fid, 1, 'long');	
    Offset_var_lenght_block = fread(fid, 1, 'long');
    Offset_ECG_block = fread(fid, 1, 'long');
    File_Version = fread(fid, 1, 'short');
    First_Name = fread(fid, 40, 'uchar');  									        								
    Last_Name = fread(fid, 40, 'uchar');  									        								
    ID = fread(fid, 20, 'uchar');  									        								
    Sex = fread(fid, 1, 'short');
    Race = fread(fid, 1, 'short');
    Birth_Date = fread(fid, 3, 'short');	
    Record_Date = fread(fid, 3, 'short');	
    File_Date = fread(fid, 3, 'short');	
    Start_Time = fread(fid, 3, 'short');	
    ihsneHeader.nbLeads = fread(fid, 1, 'short');
    Lead_Spec = fread(fid, 12, 'short');	
    Lead_Qual = fread(fid, 12, 'short');	
    ihsneHeader.Resolution = fread(fid, 12, 'short');	
    Pacemaker = fread(fid, 1, 'short');	
    Recorder = fread(fid, 40, 'uchar');
    ihsneHeader.Sampling_Rate = fread(fid, 1, 'short');	
    Proprietary = fread(fid, 80, 'uchar');
    Copyright = fread(fid, 80, 'uchar');
    Reserved = fread(fid, 88, 'uchar');
    
    % read variable_length block
    if Var_length_block_size > 0
        varblock = fread(fid, Var_length_block_size, 'char');
    end
    
    % read first location
    loc = fread(fid, 1, 'int32');
    
    % get number of beats
    currentPosition = ftell(fid);
    fseek(fid, 0, 'eof');
    endPosition = ftell(fid);
    numBeat = (endPosition-currentPosition)/4;
    fseek(fid, currentPosition - endPosition, 'eof');
    
    for i = 1:numBeat
        Ann(i) = fread(fid, 1, 'uchar');
        internalUse = fread(fid, 1, 'uchar');
        RR(i) = fread(fid, 1, 'uint16');
        loc = RR(i) + loc;
        Rloc(i) = loc;
        if ((endPos ~=0) && (Rloc(i) > endPos))
            break;
        end
    end
     
    fclose(fid);
 else
     ihsneHeader = [];
     Ann = [];
     Rloc = [];
     RR = [];
 end
