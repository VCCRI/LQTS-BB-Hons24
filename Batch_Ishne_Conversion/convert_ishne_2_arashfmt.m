function [ output_args ] = convert_ishne_2_arashfmt( infnamestr, outfnamestr, scale_factor, first_leadind )
%CONVERT_ISHNE_2_ARASHFMT this function is converting the ECGs in ishne
%format to arash format

    if (nargin == 2)
        scale_factor=1;
        first_leadind=1;
    end
    
    % Read files using file readers from THEW website
    Format_Description = 'This is a new format to storing ECG information- implemented by Arash Sadrieh';
    [arash_Header, raw_ecgSig] = read_ishne(infnamestr,0,0);
    annfnamestr = strrep(infnamestr,'.ecg','.ann');
    [tmp_Header, Ann, Rloc, RR] = read_binAnn(annfnamestr,0);
    
    raw_ecgSig = raw_ecgSig .* scale_factor;

    % V Swapping 'first lead' signal data indicated by first_leadind with the signal data in the first index
    temp = raw_ecgSig(:, 1);
    raw_ecgSig(:, 1) = raw_ecgSig(:, first_leadind);
    raw_ecgSig(:, first_leadind) = temp;
    
    save(outfnamestr,'-v7.3');

end

