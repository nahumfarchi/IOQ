classdef MESH_IO
    %MESH_IO Summary of this class goes here
    %   Detailed explanation goes here
    % TODO: add woff, wvtk_maya
    % TODO: check wavi, wmpeg4
    
    properties
    end
    
    methods (Static)
        
        %%%%%%%%
        % Read %
        %%%%%%%%

        function [X,T] = robj(filename)
            OBJ = read_wobj( filename );
            X = OBJ.vertices;
            T = OBJ.objects(1).data.vertices;
        end
        
        function [X,T] = roff(filename)
            
            fid = fopen(filename,'r');
            if( fid==-1 )
                error('Cannot open the file.');
                return;
            end
            
            str = fgets(fid);   % -1 if eof
            
            if strcmp(str(1:4), 'COFF')
                [X,T,~] = readCoff(filename,4); % assume 4 color channels
                return;
            end
            
            if ~strcmp(str(1:3), 'OFF')
                error('The file is not a valid OFF one.');
            end
           
            str = fgets(fid);
            sizes = sscanf(str, '%d %d', 2);
            while length(sizes) ~= 2
                str = fgets(fid);
                sizes = sscanf(str, '%d %d', 2);
            end
            nv = sizes(1);
            nf = sizes(2);
            
            % Read vertices
            [X,cnt] = fscanf(fid,'%lf %lf %lf\n', [3,nv]);
            if cnt~=3*nv
                warning('Problem in reading vertices.');
            end
            X = X';
            
            [T,cnt] = fscanf(fid,'3 %ld %ld %ld\n', [3,inf]);
            T = T'+1;
            
            fclose(fid);
        end
        
        function [ V ] = rnrosy(filename)
            fid = fopen(filename,'r');
            if( fid==-1 )
                error('Cannot open the file.');
                return;
            end
            
            str = fgets(fid);
            sizes = sscanf(str, '%d %d', 2);
            nv = sizes(1);
            
            [V,~] = fscanf(fid,'%lf %lf %lf', [3,nv]);
            
            fclose(fid);
        end
        
        function [ XF ] = rffield(filename)
            fid = fopen( [filename '.ffield'], 'r' );
            if( fid == -1 )
                error('Cannot open the file.');
                return;
            end
            
            str = fgets(fid); % Thresholds
            str = fgets(fid); % ??
            str = fgets(fid); % number of faces
            nf = sscanf(str,'%d');
            str = fgets(fid); % k1 k2 ...
            
            % read cross field
            [XF,cnt] = fscanf(fid,'%lf %lf %lf %lf %lf %lf %lf %lf\n', [8,nf]);
            XF = XF'; XF = XF(:,3:end);
            
            fclose( fid );
        end
        %%%%%%%%%
        % Write %
        %%%%%%%%%

        % image files
        function wfigs(filename,mesh,F,varargin)
        
            p = inputParser;
            addOptional(p,'Titles',[]);
            addOptional(p,'Colormap','jet');
            addOptional(p,'CAxis','auto');
            addOptional(p,'CAxisV',[]);
            addOptional(p,'Plot','func');
            addOptional(p,'Montage',1);
            addOptional(p,'View',[]);
            addOptional(p,'Zoom',1);
            addOptional(p,'Resolution',1024);
            addOptional(p,'Camera',[]);
            addOptional(p,'OpenGL',0);
            addOptional(p,'LightPos',[0 0 1]);
            addOptional(p,'UVt',[]);
            addOptional(p,'FUVt',[]);
            addOptional(p,'SingularPts',0);
            addOptional(p,'SingularPtsR',20);
            addOptional(p,'AspectRatio',1);
            addOptional(p,'EdgeColor','k');
            parse(p,varargin{:});
            
            if strcmp(p.Results.Plot,'mesh') == 1, sz = 1; 
            elseif strcmp(p.Results.Plot,'func') == 1, sz = size(F,2); end;

            if size(F,1) == 1
                F = repmat(F,mesh.nv,1);
            end
            for i = 1:sz
                figure;
                if strcmp(p.Results.Plot,'mesh') == 1
                    
                    MESH_VIS.mesh( mesh, ...
                                   'CW',F,...
                                   'UVt',p.Results.UVt,...
                                   'FUVt',p.Results.FUVt,...
                                   'SingularPts',p.Results.SingularPts,...
                                   'SingularPtsR',p.Results.SingularPtsR,...
                                   'EdgeColor',p.Results.EdgeColor);
                    
                elseif strcmp(p.Results.Plot,'func') == 1; 
                    
                    MESH_VIS.func( mesh, F(:,i),...
                                   'Colormap', p.Results.Colormap );
                    colorbar off;
                    if strcmp(p.Results.CAxis,'scaled') == 1
                        ca = [min(min(F)),max(max(F))];
                        caxis( ca );
                    elseif strcmp(p.Results.CAxis,'manual') == 1
                        caxis(p.Results.CAxisV);
                    end
                    
                end
                
                if isempty(p.Results.Camera) == 0
                    MESH_VIS.set_camera(gca,p.Results.Camera);
                end
                
                if p.Results.OpenGL == 1
                    lighting phong;
                    material dull;
%                     light('Position',p.Results.LightPos,...
%                           'Style','local','Color',[.8,.8,.8]);
                    camlight('headlight');
                end
                
                if numel(p.Results.Titles) > 0
                    title(p.Results.Titles{i},'FontSize',25,...
                          'interpreter','latex');    
                end
                
%                 zoom(p.Results.Zoom);
                
%                pngname = sprintf('%s_%04d.png',filename,i);
                pngname = filename;

                dpi = get(0, 'ScreenPixelsPerInch');
                in = p.Results.Resolution/dpi;
                
                set(gcf,'color','w');
                fig = gcf;
                fig.PaperUnits = 'inches';
                fig.PaperPosition = [0 0 p.Results.AspectRatio*in in];
                fig.PaperPositionMode = 'manual';
                set(gcf, 'InvertHardCopy', 'off');
                print(pngname,'-dpng','-r0')
                RemoveWhiteSpace([], 'file', [pngname '.png']);
                close all;
            end
            
            % montage
            if p.Results.Montage == 1
                rows = 5;
                cols = size(F,2) / rows;

                A = [];
                for i = 1:rows
                    b = [];
                    for j = 1:cols
                        pngname = sprintf('%s_%04d.png',filename,cols*(i-1)+j);
                        a = imread(pngname);
                        b = cat(2,b,a);
                    end
                    A = cat(1,A,b);
                end
                imwrite(A,[filename '.png']);
                delete( sprintf('%s_*.png',filename) );
                close all;
            end
%             set(gca,'FontSize',15)
        end
        
        % geometry files
        function wobj(filename,X,T)
            OBJ.vertices = X;
            OBJ.objects(1).data.vertices = T;
            OBJ.objects(1).type = 'f';
            write_wobj( OBJ, filename );
        end
        
        function woff(filename,X,T,varargin)
            fid = fopen(filename,'w');
            if( fid==-1 )
                error('Cannot open the file.');
                return;
            end

            nv = size(X,1);
            nf = size(T,1);

            fprintf(fid, 'OFF\r\n');
            fprintf(fid, '%d %d 0\r\n', nv, nf);

            fprintf(fid, '%.9f %.9f %.9f\r\n', X');
            fprintf(fid, '3 %d %d %d\r\n', (T-1)');

            fclose(fid);
        end
        
        function wnrosy(filename,M,V,varargin)
            p = inputParser;
            addOptional(p,'n',1);
            parse(p,varargin{:});
            
%             fid = fopen(filename,'w');
%             if( fid==-1 )
%                 error('Cannot open the file.');
%                 return;
%             end
%             
%             nv = size(V,1);
%             fprintf(fid,'%d %d\n', nv, p.Results.n);
%             fprintf(fid,'%.9f %.9f %.9f\r\n', V');
%             
%             fclose(fid);
            % Based on Danielle's code
            export_to_lic(filename,M,V);
        end
        
        function wvtk(filename,X,T,varargin)
            
            p = inputParser;

            addOptional(p,'Xn',[]);
            addOptional(p,'func',[]);
            addOptional(p,'func_str','vc');
            addOptional(p,'timestep',0);
            addOptional(p,'CM',[]);

            parse(p,varargin{:});
            
            nv = length(X);
            nf = length(T);

            fid = fopen(filename,'w');

            % standard header
            fprintf(fid,'# vtk DataFile Version 3.0\n');
            fprintf(fid,'vtk output\n');
            fprintf(fid,'ASCII\n');
            fprintf(fid,'DATASET POLYDATA\n');

            % if available, add time annotation
            if isempty( p.Results.timestep ) == 0
                fprintf(fid,'FIELD FieldData 1\n');
                fprintf(fid,'TIME 1 1 double\n');
                fprintf(fid,'%g\n',p.Results.timestep);
            end

            % geometry data
            fprintf(fid,'POINTS %d float\n', nv);
            fprintf(fid,'%g %g %g\n', X');
            fprintf(fid,'POLYGONS %d %d\n', nf, 4*nf);
            fprintf(fid,'3 %d %d %d\n', T'-1);
            fprintf(fid,'\n');

            % add point data
            if isempty( p.Results.func ) == 0 || isempty( p.Results.Xn ) == 0
                fprintf(fid,'POINT_DATA %d\n', nv);

                if isempty( p.Results.func ) == 0
                    fprintf(fid,'SCALARS %s float 1\n', p.Results.func_str);
                    lt = 'default'; 
                    if isempty( p.Results.CM ) == 0; lt = 'lt'; end;
                    fprintf(fid,'LOOKUP_TABLE %s\n', lt);
                    fprintf(fid,'%g\n', p.Results.func);
                end

                if isempty( p.Results.Xn ) == 0
                    fprintf(fid,'NORMALS point_normals float\n');
                    fprintf(fid,'%g %g %g\n', p.Results.Xn');
                end

                fprintf(fid,'\n');
            end

            if isempty( p.Results.CM ) == 0
                fprintf(fid,'LOOKUP_TABLE lt %d\n', length(p.Results.CM));
                fprintf(fid,'%g %g %g 1\n', p.Results.CM');
                fprintf(fid,'\n');
            end

            fclose(fid);
        end
        
        % video files
        function wavi(dirname,filename,varargin)
            
            p = inputParser;
            addOptional(p,'FrameRate',24);
            addOptional(p,'Quality',100);
            addOptional(p,'ext','png');
            parse(p,varargin{:});
            
            outputVideo = VideoWriter( [dirname filename(1:end-1) '.avi'] );
            
            outputVideo.FrameRate = p.Results.FrameRate;
            outputVideo.Quality = p.Results.Quality;

            open(outputVideo);
            
            D = dir([dirname filename '*.' p.Results.ext]);
            N = length(D);
            for i = 1:N
                img = imread( [dirname D(i).name], ...
                              'BackgroundColor', [1,1,1,] );
                writeVideo(outputVideo,img);
            end
            
            close(outputVideo);
        end
        
        function wmpeg4(dirname,filename,varargin)
            p = inputParser;
            addOptional(p,'FrameRate',24);
            addOptional(p,'Quality',100);
            addOptional(p,'ext','png');
            parse(p,varargin{:});
            
            outputVideo = VideoWriter( [dirname filename(1:end-1) '.mp4'], ...
                                       'MPEG-4' );
            
            outputVideo.FrameRate = p.Results.FrameRate;
            outputVideo.Quality = p.Results.Quality;

            open(outputVideo);
            
            D = dir([dirname filename '*.' p.Results.ext]);
            N = length(D);
            for i = 1:N
                img = imread( [dirname D(i).name],'BackgroundColor','none');
                
                % manual padding to a multiple of 2
                sz = size(img);
                if mod(sz(1),2) == 1; img = [img; 255*ones(1,sz(2),3)]; end;
                
                sz = size(img);
                if mod(sz(2),2) == 1; img = [img 255*ones(sz(1),1,3)]; end;
                
                writeVideo(outputVideo,img);
            end
            
            close(outputVideo);
        end
       
        function montage( filename, imgnames )
            rows = 1; cols = 4;

            A = [];
            for i = 1:rows
                b = [];
                for j = 1:cols
                    a = imread( [imgnames{(i-1)*rows+j} '.png'] );
                    b = cat(2,b,a);
                end
                A = cat(1,A,b);
            end
            imwrite(A,[filename '.png']);
            delete( sprintf('%s_*.png',filename) );
            close all;
        end
    end
    
end

