function results=RandomSpheres(varargin)
% 
% RandomSpheres(dimensions of x and y,sphere diameter, diameter probability,sphere number, pixel_distance)
%--------------------------------------------------------------------------
% Date: 2017-07-12
% Author: Goran Lovric
%--------------------------------------------------------------------------
% clear;clc;close all;

%--------------------------------------------------------------------------
% 1.) Constants
%--------------------------------------------------------------------------

if nargin ~= 5
    error('Incorrect input.');
    return;
end
dimensions=varargin{1};     % velicina povrsine u kojoj su uzorci
d_balls=varargin{2};  % promjeri sfera [m]
p_balls=varargin{3};    % vjerojatnosti pojedinih promjera
n_balls=varargin{4};  % ukupan broj sfera koje zelimo imati unutar volumena
delta=varargin{5};      % razmak medu cesticama

R_balls  = round((d_balls/2)/delta(2));       % [px] polje radijusa kuglica
resolution = [dimensions(1),dimensions(2),max(R_balls*2)];       % [px] broj piksela u svakom 3D smjeru
dataname = 'fullValidation';


%--------------------------------------------------------------------------
% 2.) Random distribution
%--------------------------------------------------------------------------
B_mat = zeros(resolution(1),resolution(2),resolution(3));
kk = 1

while kk < n_balls+1;
    MC_fac = rand(1);
    
    R_ind = 1+floor( rand(1)*length(R_balls) );    % bira nasumicno indeks radijusa kuglice koju ce sljedecu napraviti
    R = R_balls(R_ind);     % pridaje R nasumicno odabrani radijus
    p_1 = p_balls(R_ind)/100;   % postotak zastupljenosti kuglica radijusa R
    
    if MC_fac > p_1     % kontrolira distribuciju
        continue
    end
    
    r_pos = round (1 + rand(1,3) .* (resolution-1) );       % nasumicno daje mjesto sredista kuglice u 3D
    
    quit = false;       % flag varijabla
    
    R_mod = R + 4;   % for increasing distance between balls and borders    % uveca radijus kuglice za cetiri
    
    % (1) We test whether the ball is within the volume
    if r_pos(1) - R_mod < 1 || r_pos(2) - R_mod < 1 || r_pos(3) - R_mod < 1
        continue
    end
	
	if r_pos(1) + R_mod > resolution(1) || r_pos(2) + R_mod > resolution(2) || r_pos(3) + R_mod > resolution(3)
        continue
    end
    
    % (2) We test whether we have space to place the ball in the matrix    
    for ii = -R_mod:R_mod
        for jj = -R_mod:R_mod
            for ll = -R_mod:R_mod
                xx = r_pos(1)+ii;
                yy = r_pos(2)+jj;
                zz = r_pos(3)+ll;
                if sqrt (ii.^2 + jj.^2 + ll.^2) > R_mod;
                    continue
                end
                if B_mat(xx,yy,zz) == 1;
                    quit = true;
                    disp('Overlapping balls!')
                    break;
                end
            end
            if quit
                break;
            end
        end
        if quit
            break;
        end
    end
    
    % (3) If we have space, then we place the balls it
    if quit == false;
        for ii = -R:R%(R-1)
            for jj = -R:R%(R-1)
                for ll = -R:R%(R-1)
                xx = r_pos(1)+ii;
                yy = r_pos(2)+jj;
                zz = r_pos(3)+ll;
                if sqrt (ii.^2 + jj.^2 + ll.^2) < R;
                    B_mat(xx,yy,zz) = 1;
                end
                end
            end
        end
        R_vec(kk) = R;
        kk = kk + 1
    end
end
results=zeros(resolution(1),resolution(2));
% for ii=1:resolution(1)
%     for jj=1:resolution(2)
%         for ll=1:resolution(3)
%             results(ii,jj)=results(ii,jj)+B_mat(ii,jj,ll);
%         end
%     end
%     ii
% end
results=sum(B_mat,3);
results=results*delta(2);


%--------------------------------------------------------------------------
% 3.) Folder manipulation
%--------------------------------------------------------------------------
target_dir = dataname;

if isdir(target_dir)
    rmdir(target_dir,'s'); % delete previously created files
end
if exist(target_dir,'dir') ~= 7
    mkdir(target_dir);
end


%--------------------------------------------------------------------------
% 4.) Statistical analysis
%--------------------------------------------------------------------------
% xvalues = min(R_vec):max(R_vec);
% [nelements,centers] = hist(R_vec,xvalues);
% bar(centers,nelements)
% 
% A = [floor(centers); nelements];
% 
% fileID = fopen('quant_data.txt','w');
% fprintf(fileID,'%6s %12s\n','R [px]','Amount');
% fprintf(fileID,'%6.2f %12.0f\n',A);
% fclose(fileID);


%--------------------------------------------------------------------------
% 5.) Save images
%--------------------------------------------------------------------------
zspace = '0000';
for ind = 1:size(B_mat,3)
    zer=zspace(1:end-length(num2str(ind)));
	file=sprintf('%s/%s%d.tif',target_dir,zer,ind);
	img = uint8( (2^8-1).*(B_mat(:,:,ind)) ); % save as 8bit
	imwrite(img,file,'tif');
end
end