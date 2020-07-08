%% create/save rotation matrix for L=5 spherical harmonics

% J5anal = [  0               0           0               0           0           3*sqrt(14)/16   0               -sqrt(30)/8     0               sqrt(10)/16     0               ; ...
%         0               -sqrt(3)/2  0               1/2         0           0               0               0               0               0               0               ; ...
%         0               0           0               0           0           sqrt(70)/16     0               sqrt(6)/8       0               -9*sqrt(2)/16   0               ; ...
%         0               1/2         0               sqrt(3)/2   0           0               0               0               0               0               0               ; ...
%         0               0           0               0           0           sqrt(15)/8      0               sqrt(7)/4       0               sqrt(21)/8      0               ; ...
%         3*sqrt(14)/16   0           sqrt(70)/16     0           sqrt(15)/8  0               0               0               0               0               0               ; ...     
%         0               0           0               0           0           0               1/8             0               sqrt(42)/16     0               sqrt(210)/16    ; ...
%         -sqrt(30)/8     0           sqrt(6)/8       0           sqrt(7)/4   0               0               0               0               0               0               ; ...
%         0               0           0               0           0           0               sqrt(42)/16     0               13/16           0               -3*sqrt(5)/16   ; ...
%         sqrt(10)/16     0           -9*sqrt(2)/16   0           sqrt(21)/8  0               0               0               0               0               0               ; ...
%         0               0           0               0           0           0               sqrt(210)/16    0               -3*sqrt(5)/16   0               1/16            ];
%     
% J5anal*J5anal'
    
%     save('J5.mat','J5')
    
%% create/save rotation matrix for L=3 spherical harmonics

J3 = [  0               0           0               sqrt(10)/4  0               -sqrt(6)/4      0               ; ...
        0               1           0               0           0               0               0               ; ...
        0               0           0               sqrt(6)/4   0               sqrt(10)/4      0               ; ...
        sqrt(10)/4      0           sqrt(6)/4       0           0               0               0               ; ...
        0               0           0               0           -1/4            0               -sqrt(15)/4     ; ...
        -sqrt(6)/4      0           sqrt(10)/4      0           0               0               0               ; ...     
        0               0           0               0           -sqrt(15)/4     0               1/4             ];
    
J3*J3'
    
%     save('J3.mat','J3')    
    
%% G-matrices

kk3 = 1:7;
kk4 = 1:9;

G3z = genGLz(3);
G4z = genGLz(4);

G3y = genGLy(3);
G4y = genGLy(4);

J4a = G3y * J3 * inv( G3z(2:(end-1),:) );
J4 = zeros( size(J4a,1), size(J4a,1) );
J4( 1:end, 2:(end-1) ) = J4a;
J4( 2:(end-1), 1:end ) = J4a';
J4(end, end) = 1/8;

J5a = G4y * J4 * inv( G4z(2:(end-1),:) );
J5 = zeros( size(J5a,1), size(J5a,1) );
J5( 1:end, 2:(end-1) ) = J5a;
J5( 2:(end-1), 1:end ) = J5a';
J5(end, end) = 1/16;

save('./precomp/J3.mat','J3')  ;  
save('./precomp/J5.mat','J5') ;   