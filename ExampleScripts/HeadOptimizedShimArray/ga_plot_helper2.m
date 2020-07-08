function [ state ] = ga_plot_helper2( optinfo, state, flag, Mopt, D1comp_vol, Btarg_vol, ROImsk_or_vol, cm_rb )
%     output = 1;
%     if strcmp(flag,'init')
%         return;
%     end
    [vmin imin] = min(state.Score);
    Mplot = state.Population(imin,:);
    plot_bfield_Btarg(Mplot', optinfo, state, D1comp_vol*Mopt(Mplot)'+Btarg_vol, ROImsk_or_vol, cm_rb, 1 );

end