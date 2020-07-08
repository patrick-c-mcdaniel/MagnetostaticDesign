function [ xx, yy, zz, Bx, By, Bz, Bm ] = BSdatread( fname )
%UNTITLED Summary of this function goes here
%   usage: [ xx, yy, zz, Bx, By, Bz, Bm ] = BSdatread( fname )
    tmp = dlmread( fname );
    xx = tmp(:,1);
    yy = tmp(:,2);
    zz = tmp(:,3);
    Bx = tmp(:,4);
    By = tmp(:,5);
    Bz = tmp(:,6);
    Bm = tmp(:,7);
end

