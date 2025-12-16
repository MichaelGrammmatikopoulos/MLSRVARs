function tcode_full = transformation_code_creator(var_mnemonics)

tcode_full = cell(size(var_mnemonics));

    for n = 1 : length(var_mnemonics)
        switch var_mnemonics{n}
            case {'M2SL'}
                tcode_full{n} = 6;
            case {'CPIAUCSL','GDP','RPI','PINCOME','PCECC96', ...
                    'PPIACO','INDPRO','PAYEMS'}
                tcode_full{n} = 5;
            case {'HOUST'}
                tcode_full{n} = 4;
            case {'PSAVERT','CUMFNS','MCUMFN','EXUSUK'}
                tcode_full{n} = 2;
            case {'UNRATE',"PREMIA",'FEDFUNDS', 'TB3MS', 'TB6MS', 'GS1', 'GS3', 'GS5', 'GS10', 'BAA', 'GS20'}
                tcode_full{n} = 1;
        end % switch
    end % for n
