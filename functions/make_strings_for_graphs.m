% if CHOOSE_VAR == 2
%     which_VAR = 'financial_VAR';
% elseif CHOOSE_VAR == 4
    which_VAR = 'full_VAR';
% end

% 
% if CHOOSE_VAR==1 
%     title_VAR = 'Small VAR ';
% end
% 
% if CHOOSE_VAR==2
%     pretty_names_short(5)="6M";
%     pretty_names_short(6)="1Y";
%     pretty_names_short(7)="5Y";
%     pretty_names_short(8)="10Y";
%     title_VAR = 'Financial VAR ';
% end
% 
% if CHOOSE_VAR==3 
%     title_VAR = 'Macro VAR ';
% end

% if CHOOSE_VAR==4 
    pretty_names_short(14)="3M";
    pretty_names_short(15)="6M";
    pretty_names_short(16)="1Y";
    pretty_names_short(17)="3Y";
    pretty_names_short(18)="5Y";
    pretty_names_short(19)="10Y";
    title_VAR = 'Full VAR ';
% end