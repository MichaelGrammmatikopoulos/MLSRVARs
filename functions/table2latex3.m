function [latexTableContent, actualFilename] = table2latex3(T, filename, tbl_title, tbl_subtitle, table_id, all_models)
    % Initialize the output string
    latexTableContent = '';
    actualFilename = ''; % Initialize filename output

    % Error detection and default parameters
    write_to_file = true; % Flag to control writing to file

    if nargin < 2 || isempty(filename)
        filename = ''; % No file output if filename is empty or not provided
        write_to_file = false;
    elseif ~ischar(filename) && ~isstring(filename)
        error('The output file name must be a string or char array.');
    else
        filename = char(filename); % Ensure it's a char array for consistency
        if ~strcmp(filename(max(1,end-3):end), '.tex')
            filename = [filename '.tex'];
        end
        actualFilename = filename; % Store the resolved filename for output
    end

    if nargin < 1, error('Not enough parameters.'); end
    if ~istable(T), error('Input must be a table.'); end

    % Parameters
    n_col = (size(T,2)-1)/3+1;
    col_spec = ['l'];
    horizon_row = [' '];
    for c = 2:n_col
        col_spec = [col_spec 'lll|'];
        horizon_row = [horizon_row '& h = 1 & h = 4 & h = 8'];
    end
    unique_vars = strrep(strrep(strrep(string(T.Properties.VariableNames),'h1',''),'h4',''),'h8','');
    unique_vars=[(unique_vars(1)) (unique_vars(2:3:end))];
    col_names = append(strjoin(unique_vars(1:2),' & & '),' & & & ',strjoin(unique_vars(3:end), ' & & & '));
    row_names = T.Properties.RowNames;
    if ~isempty(row_names)
        col_spec = ['l' col_spec];
        col_names = ['& ' col_names];
    end

    % Start building the LaTeX string
    % Use sprintf to format lines and concatenate
    latexTableContent = [latexTableContent, sprintf('\\begin{flushleft} \n')];
    latexTableContent = [latexTableContent, sprintf('\\begin{table} \n')];
    latexTableContent = [latexTableContent, sprintf('\\caption{%s} \n', tbl_title)];
    latexTableContent = [latexTableContent, sprintf('\\begin{adjustbox}{tabular = %s, scale = {0.6}, center }\n', col_spec)];
    if table_id==2
        latexTableContent = [latexTableContent, sprintf('\\  & \\multicolumn{4}{c}{MSFE} & & \\multicolumn{4}{c}{MAFE} \\\\ \n')];
        latexTableContent = [latexTableContent, sprintf('\\cline{2-5} \n')];
        latexTableContent = [latexTableContent, sprintf('\\cline{7-10} \n')];
        col_names = 'Model & RGDP & Unempl. & CPI & 3M-Tbill & & RGDP & Unempl. & CPI & 3M-Tbill';
    end
    latexTableContent = [latexTableContent, sprintf('%s \\\\ \n', col_names)];
    latexTableContent = [latexTableContent, sprintf('%s \\\\ \n', horizon_row)];
    latexTableContent = [latexTableContent, sprintf('\\hline \n')];

    % Writing the data to the string
    n_col = size(T,2);
    try
        for row = 1:size(T,1)
            % if(contains(all_models(row+1),"SR"))
            %     latexTableContent = [latexTableContent, sprintf('\\rowcolor{LightCyan} \n')];
            % end
            temp = cell(1, n_col);
            for col = 1:n_col
                value = T{row,col};
                if isstruct(value), error('Table must not contain structs.'); end
                while iscell(value), value = value{1,1}; end
                if isinf(value), value = '$\\infty$'; end % Corrected backslash for LaTeX
                % Clean non-printable characters
                value = regexprep(string(value), '[^\x20-\x7E]', ''); % Ensure value is string for regexprep
                temp{1,col} = char(string(value)); % Convert value to string, then char array
            end
            if ~isempty(row_names)
                temp = [{row_names{row}}, temp];
            end
            latexTableContent = [latexTableContent, sprintf('%s \\\\ \n', strjoin(temp, ' & '))];
        end
    catch ME
        warning('Error processing row %d, column %d: %s\n', row, col, ME.message);
        warning('Table size: %d rows, %d columns\n', size(T,1), size(T,2));
        % Attempt to convert row data to cellstr safely for display
        rowDataForWarning = T{row,:};
        if iscell(rowDataForWarning)
            warning('Row data: %s\n', strjoin(rowDataForWarning, ', '));
        elseif isnumeric(rowDataForWarning) || islogical(rowDataForWarning)
             warning('Row data: %s\n', strjoin(cellstr(num2str(rowDataForWarning(:))), ', '));
        else
            warning('Row data: (Unable to display, unexpected type)');
        end
        error('Unknown error. Make sure that table only contains chars, strings or numeric values.');
    end

    % Closing the LaTeX string content
    latexTableContent = [latexTableContent, sprintf('\\hline \n')];
    latexTableContent = [latexTableContent, sprintf('\\end{adjustbox} \n')];
    latexTableContent = [latexTableContent, sprintf('\\vspace{0.5cm} \n')];
    latexTableContent = [latexTableContent, sprintf('\\caption*{%s} \n',tbl_subtitle)];
    latexTableContent = [latexTableContent, sprintf('\\end{table}\n')]; % Added newline for consistency
    latexTableContent = [latexTableContent, sprintf('\\end{flushleft}\n')]; % Added newline for consistency

    % Optional: Write to file if filename was provided
    if write_to_file
        fileID = fopen(filename, 'w');
        if fileID == -1
            error('Could not open file for writing: %s', filename);
        end
        fprintf(fileID, '%s', latexTableContent);
        fclose(fileID);
        % fprintf('LaTeX table saved to: %s\n', filename);
    end
end