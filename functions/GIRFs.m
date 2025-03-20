% forecast simulation
theseShocks_1 = (L\eyeN) * nushocks_irf(:,:,nn);
fcstX0      = Xjumpoff;

for fhn = 1 : irf_forecast_horizon
    fcstXdraw    = fcstA * fcstX0 + fcstB * theseShocks_1(:,fhn);
    ydraw        = fcstXdraw(ndxfcstY);

    % collect draw
    fcstYdraws_irf(:,fhn,nn,thisdraw) = ydraw;
    
    % prepare next iteration
    fcstX0                   = fcstXdraw;
    % censor shadow rates only for SR models
    if contains(prior_type,'srp','ignorecase',true)
        fcstX0(ndxfcstActual)    = max(fcstX0(ndxfcstShadow), ELBbound');
    end
end

% plus simulation
nushocks_irf_2 = nushocks_irf;
nushocks_irf_2(ndx_variable_to_shock,1,nn) = IRF1scale;
theseShocks_2      = (L\eyeN) * nushocks_irf_2(:,:,nn);
fcstX0           = Xjumpoff;

for fhn = 1 : irf_forecast_horizon
    fcstXdraw    = fcstA * fcstX0 + fcstB * theseShocks_2(:,fhn);
    ydraw        = fcstXdraw(ndxfcstY);

    % collect draw
    fcstYdraws1plus(:,fhn,nn,thisdraw) = ydraw;

    % prepare next iteration
    fcstX0                   = fcstXdraw;
    % censor shadow rates only for SR models
    if contains(prior_type,'srp','ignorecase',true)
        fcstX0(ndxfcstActual)    = max(fcstX0(ndxfcstShadow), ELBbound');
    end
end

% minus simulation
nushocks_irf_3 = nushocks_irf;
nushocks_irf_3(ndx_variable_to_shock,1,nn) = -1 * IRF1scale;
theseShocks_3      = (L\eyeN) * nushocks_irf_3(:,:,nn);
fcstX0           = Xjumpoff;

for fhn = 1 : irf_forecast_horizon
    fcstXdraw    = fcstA * fcstX0 + fcstB * theseShocks_3(:,fhn);
    ydraw        = fcstXdraw(ndxfcstY);

    % collect draw
    fcstYdraws1minus(:,fhn,nn,thisdraw) = ydraw;

    % prepare next iteration
    fcstX0                   = fcstXdraw;
    % censor shadow rates only for SR models
    if contains(prior_type,'srp','ignorecase',true)
        fcstX0(ndxfcstActual)    = max(fcstX0(ndxfcstShadow), ELBbound');
    end
end