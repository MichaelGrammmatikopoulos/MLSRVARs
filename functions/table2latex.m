% ----------------------------------------------------------------------- %
% Function table2latex(T, filename) converts a given MATLAB(R) table into %
% a plain .tex file with LaTeX formatting.                                %
%                                                                         %
%   Input parameters:                                                     %
%       - T:        MATLAB(R) table. The table should contain only the    %
%                   following data types: numeric, boolean, char or string.
%                   Avoid including structs or cells.                     %
%       - filename: (Optional) Output path, including the name of the file.
%                   If not specified, the table will be stored in a       %
%                   './table.tex' file.                                   %  
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       LastName = {'Sanchez';'Johnson';'Li';'Diaz';'Brown'};             %
%       Age = [38;43;38;40;49];                                           %
%       Smoker = logical([1;0;1;0;1]);                                    %
%       Height = [71;69;64;67;64];                                        %
%       Weight = [176;163;131;133;119];                                   %
%       T = table(Age,Smoker,Height,Weight);                              %
%       T.Properties.RowNames = LastName;                                 %
%       table2latex(T);                                                   %                                       
% ----------------------------------------------------------------------- %
%   Version: 1.1                                                          %
%   Author:  Victor Martinez Cagigal                                      %
%   Date:    09/10/2018                                                   %
%   E-mail:  vicmarcag (at) gmail (dot) com                               %
% ----------------------------------------------------------------------- %
function table2latex(T, filename, tbl_title, tbl_subtitle, table_id, all_models)
    % Error detection and default parameters
    if nargin < 2
        filename = 'table.tex';
        fprintf('Output path is not defined. The table will be written in %s.\n', filename); 
    elseif ~ischar(filename)
        error('The output file name must be a string.');
    else
        if ~strcmp(filename(end-3:end), '.tex')
            filename = [filename '.tex'];
        end
    end
    if nargin < 1, error('Not enough parameters.'); end
    if ~istable(T), error('Input must be a table.'); end
    
    % Parameters
    n_col = size(T,2);
    col_spec = [];
    for c = 1:n_col, col_spec = [col_spec 'l']; end
    col_names = strjoin(T.Properties.VariableNames, ' & ');
    row_names = T.Properties.RowNames;
    if ~isempty(row_names)
        col_spec = ['l' col_spec]; 
        col_names = ['& ' col_names];
    end
    
    % Writing header
    fileID = fopen(filename, 'w');
    fprintf(fileID, '\\begin{table} \n', col_spec);
    fprintf(fileID, '\\caption{%s} \n', tbl_title);
    fprintf(fileID, '\\begin{adjustbox}{tabular = %s, scale = {0.7}, center }\n', col_spec);
    if table_id==2
        fprintf(fileID, '\\  & \\multicolumn{4}{c}{MSFE} & & \\multicolumn{4}{c}{MAFE} \\\\ \n');
        fprintf(fileID, '\\cline{2-5} \n');
        fprintf(fileID, '\\cline{7-10} \n');
        col_names = 'Model & RGDP & Unempl. & CPI & 3M-Tbill & & RGDP & Unempl. & CPI & 3M-Tbill';
    end
    fprintf(fileID, '%s \\\\ \n', col_names);
    fprintf(fileID, '\\hline \n');
    
    % Writing the data
    try
        for row = 1:size(T,1)
            if(contains(all_models(row+1),"srp"))
                fprintf(fileID, '\\rowcolor{LightCyan} \n');
            end
            temp = cell(1, n_col);
            for col = 1:n_col
                value = T{row,col};
                if isstruct(value), error('Table must not contain structs.'); end
                while iscell(value), value = value{1,1}; end
                if isinf(value), value = '$\infty$'; end
                % Clean non-printable characters
                value = regexprep(value, '[^\x20-\x7E]', '');
                temp{1,col} = num2str(value);
            end
            if ~isempty(row_names)
                temp = [{row_names{row}}, temp];
            end
            fprintf(fileID, '%s \\\\ \n', strjoin(temp, ' & '));
        end
    catch ME
        fprintf('Error processing row %d, column %d: %s\n', row, col, ME.message);
        fprintf('Table size: %d rows, %d columns\n', size(T,1), size(T,2));
        fprintf('Row data: %s\n', strjoin(cellstr(T{row,:}), ', '));
        error('Unknown error. Make sure that table only contains chars, strings or numeric values.');
    end
    
    % Closing the file
    fprintf(fileID, '\\hline \n');   
    fprintf(fileID, '\\end{adjustbox} \n');
    fprintf(fileID, '\\vspace{0.5cm} \n');    
    fprintf(fileID, '\\caption*{%s} \n',tbl_subtitle);
    fprintf(fileID, '\\end{table}');
    fclose(fileID);
end
