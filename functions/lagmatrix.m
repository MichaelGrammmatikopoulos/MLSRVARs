function varargout = lagmatrix(varargin)
%LAGMATRIX Create lagged time series data
%
% Syntax:
%
%   [YLag,TLag] = lagmatrix(Y,lags)
%   LagTbl = lagmatrix(Tbl,lags)
%   [...] = lagmatrix(...,param,val,...)
%
% Description:
%
%   Shift regular series in time by specified lags or leads.
%
% Input Arguments:
%
%   Y - Time series data, specified as a numObs-by-numVars numeric matrix.
%
%   Tbl - Time series data, specified as a table or timetable. Specify
%       series for Y using the 'DataVariables' parameter. Timetable data
%       must be sampled with a regular datetime time step.
%
%   lags - Data shifts, specified as an integer-valued vector of length
%       numShifts. Lags are positive (shift Y forward over the time base),
%       leads are negative (shift Y backward over the time base). Each
%       shift in lags is applied, in order, to each series in Y. Shifts of 
%       regular time series have units of one time step.
%
% Optional Input Parameter Name/Value Arguments:
%
%   NAME    VALUE
%
%   'Y0'	Presample data to back fill lagged series, specified as a
%           matrix, table or timetable. Matrix Y0 must have numVars
%           columns. Tabular Y0 selects variables with the 'DataVariables'
%           parameter. Timetables must have regular sample times preceding
%           times in Tbl. The default for presample data is NaN.
%
%   'YF'    Postsample data to front fill led series, specified as a
%           matrix, table or timetable. Matrix YF must have numVars
%           columns. Tabular YF selects variables with the 'DataVariables'
%           parameter. Timetables must have regular sample times following
%           times in Tbl. The default for postsample data is NaN.
%
%	'DataVariables' Variables to use for Y, Y0 and YF, specified as names
%           in Tbl.Properties.VariableNames. Variable names are cell
%           vectors	of character vectors, string vectors, integer vectors
%           or logical vectors. The default is all variables.
%
%   'Shape'	Part of the shifted series to appear in the outputs, specified
%           as a string or character vector. Values are 'full', 'same' and
%           'valid'. When 'Shape' is 'full', outputs contain all values in
%           Y and all presample or postsample values on an expanded time
%           base. When 'Shape' is 'same', outputs contain only values on
%           the original time base. When 'Shape' is 'valid', outputs
%           contain values for times at which all series have specified
%           (non-NaN) values. The default is 'same'.
%
% Output Argument:
%
%   YLag - Shifted Y variables, returned as a matrix. Columns are, in
%          order, all series in Y shifted by the first value in lags, all
%          series in Y shifted by the second value in lags, ..., all series
%          in Y shifted by the last value in lags. Rows depend on the value
%          of the 'Shape' parameter.
%
%   TLag - Common time base for the shifted series relative to the original
%          time base of 1,2,3,...,numObs, returned as a vector of length
%          equal to the number of observations in YLag. Series with lags
%          have higher indices; series with leads have lower indices.
%
%   LagTbl - When input is Tbl, outputs YLag and TLag are returned in
%          tabular LagTbl, the same type as Tbl.
%
% Example:
%
% 	Y = [0.1 0.4 -0.2 0.1 0.2]';
% 	lags = [2 0 -2];
% 	[YLag,TLag] = lagmatrix(Y,lags,...
%                           'Y0',[0.50; 0.75]*Y(1),...
%                           'YF',[0.75; 0.50]*Y(end),...
%                           'Shape','full')
%
% See also LINEARMODEL, FILTER.

% Copyright 2022 The MathWorks, Inc.

% Parse inputs and set defaults:

parseObj = inputParser;

parseObj.addRequired('Data',...
	@(x)validateattributes(x,{'double','table','timetable'},{'nonempty','2d'}))
    
parseObj.addRequired('lags',...
	@(x)validateattributes(x,{'double'},{'integer'})) % Do not parse for nonempty vector:
                                                      % ARIMA calls with empty lags
   
parseObj.addParameter('Y0',NaN,...
	@(x)validateattributes(x,{'double','table','timetable'},{'nonempty','2d'}))
   
parseObj.addParameter('YF',NaN,...
	@(x)validateattributes(x,{'double','table','timetable'},{'nonempty','2d'}))

parseObj.addParameter('DataVariables',[],...
    @(x)validateattributes(x,{'double','logical','cell','string'},{'vector'}));

parseObj.addParameter('Shape','same',...
	@(x)validateattributes(x,{'char','string'},{'vector'}))

try
    
  parseObj.parse(varargin{:});
  
catch ME
    
  throwAsCaller(ME)
  
end

Data = parseObj.Results.Data;
lags = parseObj.Results.lags;
Y0 = parseObj.Results.Y0;
YF = parseObj.Results.YF;
varSpec = parseObj.Results.DataVariables;
shape = validatestring(parseObj.Results.Shape,{'full','same','valid'});

% Flag output type:

TableOutput = istable(Data);
TimeTableOutput = istimetable(Data);
isTabular = TableOutput || TimeTableOutput;

% Select Y, Y0, YF with 'DataVariables':

if isnumeric(Data)
    
   	Y = Data;
    
    if ~isempty(varSpec)
        
        warning(message('econ:lagmatrix:DataVariablesUnused'))
        
    end
    
else % Tabular data

    if ~isempty(varSpec)

        try

            Y = Data(:,varSpec);
            
            if ~isnumeric(Y0) && ~ismember('Y0',parseObj.UsingDefaults)
            
                Y0 = Y0(:,varSpec);
                
            end
            
            if ~isnumeric(YF) && ~ismember('YF',parseObj.UsingDefaults)
            
                YF = YF(:,varSpec);
                
            end
            
        catch ME

            throwAsCaller(ME)

        end

    else

        Y = Data; % Default

    end

end

% Validate table and timetable inputs:

try
    
    internal.econ.TableAndTimeTableUtilities.ensureTimeSeriesTypeConsistency(Y,Y0,YF,{'Y','Y0','YF'},parseObj.UsingDefaults)
    
    internal.econ.TableAndTimeTableUtilities.isTabularFormatValid(Y,'Y')
    internal.econ.TableAndTimeTableUtilities.isTabularFormatValid(Y0,'Y0')
    internal.econ.TableAndTimeTableUtilities.isTabularFormatValid(YF,'YF')
    
    internal.econ.TableAndTimeTableUtilities.isTabularDataSinglePath(Y,'Y')
    internal.econ.TableAndTimeTableUtilities.isTabularDataSinglePath(Y0,'Y0')
    internal.econ.TableAndTimeTableUtilities.isTabularDataSinglePath(YF,'YF')
    
    internal.econ.TableAndTimeTableUtilities.isTimeTableRegular(Y,'Y')
    internal.econ.TableAndTimeTableUtilities.isTimeTableRegular(Y0,'Y0')
    internal.econ.TableAndTimeTableUtilities.isTimeTableRegular(YF,'YF')
    
    internal.econ.TableAndTimeTableUtilities.areTimeTableDateTimesConsistent(Y0,Y,'Sequential',{'Y0','Y'})
    internal.econ.TableAndTimeTableUtilities.areTimeTableDateTimesConsistent(Y,YF,'Sequential',{'Y','YF'})
    
    internal.econ.TableAndTimeTableUtilities.areTabularVariableNamesConsistent(Y0,Y,{'Y0','Y'})
    internal.econ.TableAndTimeTableUtilities.areTabularVariableNamesConsistent(Y,YF,{'Y','YF'}) 
  
catch exception
    
	exception.throwAsCaller()
  
end

% Convert tabular data:

if isTabular
    
    MetaData = Y.Properties;
    Y = table2array(Y);
    Y = double(Y);
    
   	Y0 = internal.econ.TableAndTimeTableUtilities.perVariableTabular2PerPathMatrix(Y0);
    Y0 = double(Y0);

    YF = internal.econ.TableAndTimeTableUtilities.perVariableTabular2PerPathMatrix(YF);
    YF = double(YF);
    
end
    
if isvector(Y)
    
    Y = Y(:);
    
end

lags = lags(:)';

% Get sizes:

[numObs,numVars] = size(Y);
numLags = length(lags);
maxLag = max([lags(lags >= 0),0]);
maxLead = abs(min([lags(lags <= 0),0]));
numPre = size(Y0,1);
numPost = size(YF,1);

% Check size of Y0, YF:

if iscolumn(Y0)
    
    Y0 = repmat(Y0,1,numVars);
    
elseif size(Y0,2) ~= numVars
    
    error(message('econ:lagmatrix:SizeMismatch_Y0_Y'))
    
end

if iscolumn(YF)
    
    YF = repmat(YF,1,numVars);
    
elseif size(YF,2) ~= numVars
    
    error(message('econ:lagmatrix:SizeMismatch_Y_YF'))
    
end

% Preallocate outputs:

YLag = NaN(numObs+numPre+numPost+maxLag+maxLead,numVars*numLags);
TLag = [(1-(numPre+maxLead):0)';(1:numObs)';(numObs+1:numObs+numPost+maxLag)'];
idxShift = numPre+maxLead;

% Write lagged series to YLag:

Y0YYF = [Y0;Y;YF];

for lagNum = 1:numLags

    L = lags(lagNum);
    YLag(((1-numPre+L):(numObs+numPost+L))+idxShift,...
         (1+(lagNum-1)*numVars):(lagNum*numVars)) = Y0YYF;

end

% Trim rows according to 'Shape':

switch shape
    
    case 'full'
        
        rowIdx = ((1-maxLead):(numObs+maxLag))+idxShift;      
        
    case 'same'
        
        rowIdx = (1:numObs)+idxShift;
        
    case 'valid'
        
        rowIdx = ~any(isnan(YLag),2); 
        
end

% Numeric outputs:

YLag = YLag(rowIdx,:);
TLag = TLag(rowIdx);

% Table and timetable outputs:

if isTabular
    
    % Create lagged variable names:
    
    varNames = repmat(MetaData.VariableNames,1,numLags);
    
    for i = 1:numLags
        
        if lags(i) >= 0
            
            varNames((((i-1)*numVars)+1):(i*numVars)) = ...
              strcat(['Lag',num2str(lags(i),'%-u')],varNames((((i-1)*numVars)+1):(i*numVars)));
         
        else
            
            varNames((((i-1)*numVars)+1):(i*numVars)) = ...
              strcat(['Lead',num2str(-lags(i),'%-u')],varNames((((i-1)*numVars)+1):(i*numVars)));
         
        end
        
    end    
    
    % Create LagTbl:
    
    if TableOutput
               
        % For TableOutput, RowNames of the input table are removed in the
        % output table, for consistency with TimeTableOutput.

        LagTbl = array2table([TLag,YLag],'VariableNames',[{'TLag'},varNames]);
        
        if ~isempty(MetaData.VariableDescriptions)
            LagTbl.Properties.VariableDescriptions = [{'Time Lag'},repmat(MetaData.VariableDescriptions,1,numLags)];
        else
            LagTbl.Properties.VariableDescriptions = [{'Time Lag'},repmat("",1,numVars*numLags)];
        end
        
        if ~isempty(MetaData.VariableUnits)
            LagTbl.Properties.VariableUnits = [{'Time Step'},repmat(MetaData.VariableUnits,1,numLags)];
        else
            LagTbl.Properties.VariableUnits = [{'Time Step'},repmat("",1,numVars*numLags)];
        end
        
        if ~isempty(MetaData.VariableContinuity)
            LagTbl.Properties.VariableContinuity = [{'Step'},repmat(MetaData.VariableContinuity,1,numLags)];
        else
            LagTbl.Properties.VariableContinuity = [{'Step'},repmat("unset",1,numVars*numLags)];
        end

    end
    
    if TimeTableOutput

        startTime = MetaData.StartTime+(TLag(1)-1)*MetaData.TimeStep;
        LagTbl = array2timetable(YLag,'VariableNames',varNames,...
                                 'StartTime',startTime,...
                                 'TimeStep',MetaData.TimeStep);
        
        if ~isempty(MetaData.VariableDescriptions)
            LagTbl.Properties.VariableDescriptions = repmat(MetaData.VariableDescriptions,1,numLags);
        else
            LagTbl.Properties.VariableDescriptions = cellstr(repmat("",1,numVars*numLags));
        end
        
        if ~isempty(MetaData.VariableUnits)
            LagTbl.Properties.VariableUnits = repmat(MetaData.VariableUnits,1,numLags);
        else
            LagTbl.Properties.VariableUnits = cellstr(repmat("",1,numVars*numLags));
        end
        
        if ~isempty(MetaData.VariableContinuity)
            LagTbl.Properties.VariableContinuity = repmat(MetaData.VariableContinuity,1,numLags);
        else
            LagTbl.Properties.VariableContinuity = cellstr(repmat("unset",1,numVars*numLags));
        end

    end
        
    LagTbl.Properties.Description = 'Lagged Variables';
    LagTbl.Properties.DimensionNames = MetaData.DimensionNames;

end

% Assign outputs:

if isTabular
    
    nargoutchk(0,1)
    
    varargout{1} = LagTbl;
    
else
    
    nargoutchk(0,2)
    
    varargout{1} = YLag;
    varargout{2} = TLag;

end