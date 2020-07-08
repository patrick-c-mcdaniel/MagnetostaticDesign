function [ output ] = proc_A_comp( fn_MagROI, fn_A_comp_save)
    load(fn_MagROI);

    Ndof = 2*numel(inds_mid) + 3*numel(inds_pos);

    Mred_null = zeros(Ndof, 1);

    [ Mxyz_null, Arf, A_comp ] = Mred2Mxyz_A_comp( Mred_null, inds_mid, inds_pos, inds_neg, [3 5] );

    save(fn_A_comp_save,'A_comp');
end