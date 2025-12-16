
all_models_pretty = cellstr(all_models)';
for m_i = 1:size(all_models,2)

    if contains(all_models_pretty{m_i},'minnesota','ignorecase',true)
        all_models_pretty{m_i} = strrep(all_models_pretty{m_i},"minnesota","Minnesota");
    end

    if contains(all_models_pretty{m_i},'_stvol','ignorecase',true)
        all_models_pretty{m_i} = strrep(all_models_pretty{m_i},"_stvol","-SV");
    end

    if contains(all_models_pretty{m_i},'_ssp','ignorecase',true)
        all_models_pretty{m_i} = strrep(all_models_pretty{m_i},"_ssp","-VSLR");
    end

    if contains(all_models_pretty{m_i},'_srp','ignorecase',true)
        all_models_pretty{m_i} = strrep(all_models_pretty{m_i},"_srp","-SR");
    end

    if contains(all_models_pretty{m_i},'blasso_A','ignorecase',true)
        all_models_pretty{m_i} = strrep(all_models_pretty{m_i},"blasso_A","BLASSO");
    end

    if contains(all_models_pretty{m_i},'blasso_G','ignorecase',true)
        all_models_pretty{m_i} = strrep(all_models_pretty{m_i},"blasso_G","BLASSO(G)");
    end

    if contains(all_models_pretty{m_i},'DirLap','ignorecase',true)
        all_models_pretty{m_i} = strrep(all_models_pretty{m_i},"DirLap","DL");
    end

    if contains(all_models_pretty{m_i},'ssvs','ignorecase',true)
        all_models_pretty{m_i} = strrep(all_models_pretty{m_i},"ssvs","SSVS") ;
        if contains(all_models_pretty{m_i},'_Min','ignorecase',true)
            all_models_pretty{m_i} = strrep(all_models_pretty{m_i},"_Min","-Min");
        end
    end
end

all_models_pretty = string(all_models_pretty);
%     all_models_pretty{1}    = 'Minnesota';
% all_models_pretty{2}    = 'SSVS';
% all_models_pretty{3}    = 'BLASSO';
% all_models_pretty{4}    = 'Minnesota-SSP';
% all_models_pretty{5}    = 'SSVS-SSP';
% all_models_pretty{6}    = 'BLASSO-SSP';
% all_models_pretty{7}    = 'Minnesota-SV';
% all_models_pretty{8}    = 'SSVS-SV';
% all_models_pretty{9}    = 'BLASSO-SV';
% all_models_pretty{10}   = 'Minnesota-SSP-SV';
% all_models_pretty{11}   = 'SSVS-SSP-SV';
% all_models_pretty{12}   = 'BLASSO-SSP-SV';
% all_models_pretty{13}   = 'SR-Minnesota';
% all_models_pretty{14}   = 'SR-SSVS';
% all_models_pretty{15}   = 'SR-BLASSO';
% all_models_pretty{16}   = 'SR-Minnesota-SSP';
% all_models_pretty{17}   = 'SR-SSVS-SSP';
% all_models_pretty{18}   = 'SR-BLASSO-SSP';
% all_models_pretty{19}   = 'SR-Minnesota-SV';
% all_models_pretty{20}   = 'SR-SSVS-SV';
% all_models_pretty{21}   = 'SR-BLASSO-SV';
% all_models_pretty{22}   = 'SR-Minnesota-SSP-SV';
% all_models_pretty{23}   = 'SR-SSVS-SSP-SV';
% all_models_pretty{24}   = 'SR-BLASSO-SSP-SV';