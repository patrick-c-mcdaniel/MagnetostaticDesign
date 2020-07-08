function [ output ] = proc_A_comp( fn_MagROI, fn_A_comp_save)
    load(fn_MagROI);

    Ndof = 4*numel(inds_pp);

    Mred_null = zeros(Ndof, 1);

    [ Mxyz_null, Arf, A_comp ] = Mred2Mxyz_A_comp( Mred_null, inds_pp, inds_pn, inds_nn, inds_np, [3 5] );

    save(fn_A_comp_save,'A_comp');
end