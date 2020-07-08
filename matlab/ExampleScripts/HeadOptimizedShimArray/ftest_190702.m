function [ output ] = ftest_190702( Mred, Mopt, D )
    tmp = D*Mopt(Mred');
    output = 1;
end