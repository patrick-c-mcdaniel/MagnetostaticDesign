function [ state ] = ga_plot_helper( optinfo, state, flag, Mdel, D1comp_vol, Btarg_vol, ROImsk_or_vol, cm_rb )
%     output = 1;
%     if strcmp(flag,'init')
%         return;
%     end
    [vmin imin] = min(state.Score);
    Mplot = state.Population(imin,:);
    plot_bfield_Btarg(Mplot', optinfo, state, Mdel*D1comp_vol*(Mplot)'+Btarg_vol, ROImsk_or_vol, cm_rb, 1 );

end