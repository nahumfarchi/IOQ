%% ========================================================================
%  Time the approx resistance vs the true resistance
%  ========================================================================
out_folder = 'results';
mkdir(out_folder);
EPS = 0.5;
JLFAC = 24;
USE_GPU = true;
SEED = 112;

disp('Loading mesh...')
fp = '../../../data/bunny.off';
%fp = '../../../data/bunnies/bunny_99k_faces.off';
%fp = '../../../data/ashish_nob/santa.off';
[~, meshname, ~] = fileparts(fp);
m = Mesh(fp); ne = m.nE;
m.nV

if USE_GPU
    gd = gpuDevice();
else
    gd = [];
end
%
disp('Approximating resistance distance...')
tic
Rtilde = resistance_distance(m, 1, JLFAC, USE_GPU); 
if ~isempty(gd), wait(gd); end; elapsed_Rtilde = toc
clear Rtilde

%%
tic
R = resistance_distance(m, 0, JLFAC, USE_GPU, false);
if ~isempty(gd), wait(gd); end; elapsed_R = toc
%fprintf('|R - Rtilde| = %.4g\n', norm(R - squareform(Rtilde), 'fro'));

tic
R_pdist = resistance_distance(m, 0, JLFAC, USE_GPU, true);
if ~isempty(gd), wait(gd); end; elapsed_R_pdist = toc