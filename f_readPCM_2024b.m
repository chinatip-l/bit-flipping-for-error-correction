function [n,m,maxVNd,maxCNd,VNd, CNd, VNlink, CNlink, PCM, PCMsparse] = f_readPCM_2024b(fname)
% reads binary parity check matrix in "alist" format from file FNAME and
% converts it to sparse matrix used in MATLAB routines.
% This is an interface to matrices at http://wol.ra.phy.cam.ac.uk/mackay/codes/

% Copyright (c) 1999 by Igor Kozintsev igor@ifp.uiuc.edu
% $Revision: 1.1 $ $Date: 2000/03/23 $ Bug fixed by Hatim Behairy
% $Modified version: 2024.b $Date: 2024/05/09 $modified by T. C.-Y. Chang

fid = fopen(fname);
n = fscanf(fid,'%d',1);
m = fscanf(fid,'%d',1);
maxVNd = fscanf(fid,'%d',1);
maxCNd = fscanf(fid,'%d',1); 
VNd = fscanf(fid,'%d',[1 n]);
CNd = fscanf(fid,'%d',[1 m]); 
position = zeros(n,maxVNd);
VNlink = position;
CNlink = zeros(m,maxCNd);
PCM = zeros(m,n);
for i=1:n
    
    for j=1:VNd(i)
        tmp=fscanf(fid,'%d',1);
        if j<=VNd(i)
            
            position(i,j) = tmp;
            VNlink(i,j) = position(i,j);
        end
    end
end
for i=1:m
    for j=1:CNd(i)
        tmp = fscanf(fid,'%d',1);
         
        if j<=CNd(i)
            
            position(j,i) = tmp; 
            CNlink(i,j)=position(j,i);
            PCM(i,CNlink(i,j)) = 1;
        end
    end
end

PCMsparse =sparse(PCM);
fclose(fid);