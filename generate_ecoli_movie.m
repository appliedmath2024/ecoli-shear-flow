% GENERATE_ECOLI_MOVIE Render E. coli simulation outputs as an MP4 or GIF.
% Run this script after generating the Fortran output files in ecoli_shear_flow.
% Required files: status.m, coord_pt.m, triad1.m, triad2.m, body_pt.m,
% AD.m, and body_face.m. These generated files are not bundled with the code.
% Paths are relative to this script. Output is written to A_Movies and an
% existing movie with the same name is overwritten. Rendering needs a display.
% MP4 playback uses 150 frames/s; the title shows simulation time, not playback
% time. GIF output uses the original zero-delay setting (viewer dependent).

clearvars
close all
clc

%% Input settings
project_dir = fileparts(mfilename('fullpath'));

folder_list = "ecoli_shear_flow";
legend_name = {}; % Optional labels; leave empty to hide the legend.

nFolder = numel(folder_list);

xz_switch = 1; % Swap the x and z coordinates for the displayed view.

%% Load simulation geometry and metadata
for ff = 1:nFolder
    
    data_dir = fullfile(project_dir, folder_list(ff));
    metadata = read_metadata(data_dir);

    S(ff).coord_pt = load(fullfile(data_dir, 'coord_pt.m'));
    S(ff).triad1   = load(fullfile(data_dir, 'triad1.m'));
    S(ff).triad2   = load(fullfile(data_dir, 'triad2.m'));
    S(ff).body_pt  = load(fullfile(data_dir, 'body_pt.m'));
    S(ff).AD       = load(fullfile(data_dir, 'AD.m'));
    S(ff).body_face= load(fullfile(data_dir, 'body_face.m'));
    S(ff).nflag = metadata.nflag;
    S(ff).nPt   = metadata.nPt;
    S(ff).nBd   = metadata.nBd;
    S(ff).nhook = metadata.nhook;
    S(ff).dt    = metadata.dt;
    S(ff).nSkip = metadata.nSkip;

    if(xz_switch==1)
        for k=1:S(ff).nflag
            S(ff).coord_pt = switch_cols(S(ff).coord_pt,3*(k-1)+1,3*k);
            S(ff).triad1   = switch_cols(S(ff).triad1,3*(k-1)+1,3*k);
            S(ff).triad2   = switch_cols(S(ff).triad2,3*(k-1)+1,3*k);
        end
        S(ff).body_pt = switch_cols(S(ff).body_pt,1,3);
        S(ff).AD      = switch_cols(S(ff).AD,5,7);
    end

    % Fortran writes the face-index array in column-major order.
    Bnum = size(S(ff).body_face,1)/12;
    S(ff).body_face = reshape(S(ff).body_face,Bnum*4,3);
end
for ff = 1:nFolder
    fprintf('folder = %s\n', folder_list(ff));
    fprintf('  nPt    = %d\n', S(ff).nPt);
    fprintf('  nBd    = %d\n', S(ff).nBd);
    fprintf('  nflag  = %d\n', S(ff).nflag);
    fprintf('  nhook  = %d\n', S(ff).nhook);
    fprintf('  frames = %d\n\n', floor(size(S(ff).coord_pt,1)/S(ff).nPt));
end
%% Fixed axis limits across all frames
maxx=-inf; maxy=-inf; maxz=-inf;
minx= inf; miny= inf; minz= inf;

for ff = 1:nFolder
    maxx=max([maxx, max(S(ff).coord_pt(:,1:3:end),[],'all'), max(S(ff).body_pt(:,1))]);
    maxy=max([maxy, max(S(ff).coord_pt(:,2:3:end),[],'all'), max(S(ff).body_pt(:,2))]);
    maxz=max([maxz, max(S(ff).coord_pt(:,3:3:end),[],'all'), max(S(ff).body_pt(:,3))]);

    minx=min([minx, min(S(ff).coord_pt(:,1:3:end),[],'all'), min(S(ff).body_pt(:,1))]);
    miny=min([miny, min(S(ff).coord_pt(:,2:3:end),[],'all'), min(S(ff).body_pt(:,2))]);
    minz=min([minz, min(S(ff).coord_pt(:,3:3:end),[],'all'), min(S(ff).body_pt(:,3))]);
end

miny=-1;

Fontsize=30;
Fontname='Times New Roman';
Linewidth=2;
close all

%% Movie settings
movie_type = 2; % 1: animated GIF; 2: MPEG-4 (requires encoder support).
frame_stride = 2; % Render every second saved simulation frame.
output_dir = fullfile(project_dir, 'A_Movies');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

figure(1)
set(gcf,'color','w');
set(gcf,'position',[20 20 1032 1252])

if (movie_type==2)
    M = VideoWriter(fullfile(output_dir, 'Movie_allFolders.mp4'), 'MPEG-4');
    M.FrameRate = 150;
    M.Quality = 75;
    open(M);
end

colors = lines(nFolder);

movie_num=0;
numFrames = inf;
for ff = 1:nFolder
    numFrames = min(numFrames, floor(size(S(ff).coord_pt,1)/S(ff).nPt));
end
%% Render flagella, cell body, and the trajectory projected onto the wall
for kk = 1:frame_stride:numFrames
    clf
    hold on

    movie_num = movie_num + 1;
    h = gobjects(nFolder,1);

    for ff = 1:nFolder

        nflag = S(ff).nflag;
        nPt   = S(ff).nPt;
        nBd   = S(ff).nBd;
        nhook = S(ff).nhook;

        X  = zeros(nPt,3,nflag);
        D1 = zeros(nPt,3,nflag);
        D2 = zeros(nPt,3,nflag);

        for k=1:nflag
            X(:,:,k)  = S(ff).coord_pt((kk-1)*nPt+1:kk*nPt,3*(k-1)+1:3*k);
            D1(:,:,k) = S(ff).triad1((kk-1)*nPt+1:kk*nPt,3*(k-1)+1:3*k);
            D2(:,:,k) = S(ff).triad2((kk-1)*nPt+1:kk*nPt,3*(k-1)+1:3*k);
        end

        B = S(ff).body_pt((kk-1)*nBd+1:kk*nBd,1:3);
        % Reconstruct a tube from the centerline and two material directors.
        m = 40;
        R = 0.04;
        dd = 2*pi/m;
        s = 0:dd:2*pi;

        x1 = zeros(m+1,nPt,nflag);
        x2 = zeros(m+1,nPt,nflag);
        x3 = zeros(m+1,nPt,nflag);

        for k=1:nflag
            for i=1:m+1
                x1(i,:,k) = X(:,1,k) + R*(cos(s(i))*D1(:,1,k) + sin(s(i))*D2(:,1,k));
                x2(i,:,k) = X(:,2,k) + R*(cos(s(i))*D1(:,2,k) + sin(s(i))*D2(:,2,k));
                x3(i,:,k) = X(:,3,k) + R*(cos(s(i))*D1(:,3,k) + sin(s(i))*D2(:,3,k));
            end
        end
        firstSurf = true;
        for k=1:nflag
            hs = surf(x1(:,nhook:nPt,k),x2(:,nhook:nPt,k),x3(:,nhook:nPt,k), ...
                'FaceColor',colors(ff,:), ...
                'EdgeColor','none', ...
                'FaceAlpha',0.9);

            if firstSurf
                h(ff) = hs;   % legend handle
                firstSurf = false;
            end

            if (nhook~=1)
                surf(x1(:,1:nhook,k),x2(:,1:nhook,k),x3(:,1:nhook,k), ...
                    'FaceColor',colors(ff,:), ...
                    'EdgeColor','none', ...
                    'FaceAlpha',0.2);
            end
        end
        patch('faces',S(ff).body_face,'vertices',B, ...
              'facecolor',colors(ff,:), ...
              'linestyle','none', ...
              'FaceAlpha',0.50);
        plot3(S(ff).AD(1:kk,5), S(ff).AD(1:kk,6)*0, S(ff).AD(1:kk,7), ...
              '-', 'Color',colors(ff,:), 'LineWidth',Linewidth);
    end
    % The wall is the y = 0 plane in display coordinates.
    patch([minx maxx maxx minx],[0 0 0 0],[minz minz maxz maxz], ...
          [230 230 230]/255,'FaceAlpha',1,'LineWidth',0.1)

    v1=180; 
    v2=0;
    view(v1,v2)

    light('Position',[cos(v2*pi/180)*sin(v1*pi/180), ...
                            -cos(v2*pi/180)*cos(v1*pi/180), ...
                             sin(v2*pi/180)]);

    lighting gouraud

    myfig = gcf;
    myfig.RendererMode = 'manual';

    material dull

    axis equal
    axis([minx maxx miny maxy minz maxz])

    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'ztick',[])

    title(sprintf('t = %6.4fs',S(1).dt*S(1).nSkip*(kk-1)), ...
        'fontsize',Fontsize,'fontname',Fontname);

    if(xz_switch==1)
        
        set(gca, 'ZDir', 'reverse')
    end

    set(gca,'fontsize',Fontsize,'fontname',Fontname)
    if ~isempty(legend_name)
        legend(h, legend_name, 'NumColumns', 2, 'Location', 'north', ...
            'Orientation', 'horizontal', 'Interpreter', 'none');
    end

    drawnow

    if (movie_type==1)
        f = getframe(gcf);
        imind = frame2im(f);
        [imind,cm] = rgb2ind(imind,256);

        if (movie_num==1)
            imwrite(imind,cm,fullfile(output_dir, 'Movie_allFolders.gif'), ...
                'Loopcount',inf,'DelayTime',0.0);
        else
            imwrite(imind,cm,fullfile(output_dir, 'Movie_allFolders.gif'), ...
                'WriteMode','append','DelayTime',0.0);
        end

    elseif (movie_type==2)
        f = getframe(gcf);
        writeVideo(M, f);
    end
end

if (movie_type==2)
    close(M);
end
function A = switch_cols(A, col1, col2)
    A(:, [col1 col2]) = A(:, [col2 col1]);
end

function metadata = read_metadata(data_dir)
% Load the Fortran-generated status script in a separate function workspace.
    required_files = {'status.m', 'coord_pt.m', 'triad1.m', 'triad2.m', ...
        'body_pt.m', 'AD.m', 'body_face.m'};
    for ii = 1:numel(required_files)
        input_path = fullfile(data_dir, required_files{ii});
        if ~isfile(input_path)
            error('generate_ecoli_movie:MissingInput', ...
                'Missing input: %s. Generate the Fortran outputs first.', input_path);
        end
    end
    run(fullfile(data_dir, 'status.m'));
    metadata = struct('nflag', nflag, 'nPt', nPt, 'nBd', nBd, ...
        'nhook', nhook, 'dt', dt, 'nSkip', nSkip);
end
