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

        function [X,T,UV,UVT] = robj_t(filename)
            OBJ = read_wobj( filename );
            X = OBJ.vertices;
            T = OBJ.objects(end).data.vertices;
            UV = OBJ.vertices_texture;
            UVT = OBJ.objects(end).data.texture;
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
            T = double(T);
                
            fclose(fid);
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
            addOptional(p,'Plot',0);
            addOptional(p,'Xaxis',1:size(F(:,1),1));
            addOptional(p,'Montage',1);
            parse(p,varargin{:});
            
            for i = 1:size(F,2)
                figure;
                if p.Results.Plot == 1; 
                    plot(p.Results.Xaxis,F(:,i),'LineWidth',3);
                    set(gca,'fontsize',40);
                else
                    MESH_VIS.func(mesh,F(:,i),'Colormap',p.Results.Colormap);
                end
                colorbar off;
                
                if strcmp(p.Results.CAxis,'scaled') == 1
                    ca = [min(min(F)),max(max(F))];
                    caxis( ca );
                end
                
                pngname = sprintf('%s_%04d.png',filename,i);
                
                fig = gcf;
                fig.PaperUnits = 'points';
%                 fig.PaperPosition = [0 0 1024 1024];
                fig.PaperPosition = [0 0 512 512];
                fig.PaperPositionMode = 'manual';
                print(pngname,'-dpng','-r0')
                close all;
            end
            
            % montage
            if p.Results.Montage == 1
                rows = 2;
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
            
%             title(titles{i},'FontSize',25,'interpreter','latex');    
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
        
        function wnrosy(filename,V,varargin)
            
            fid = fopen(filename,'w');
            if( fid==-1 )
                error('Cannot open the file.');
                return;
            end
            
            nv = size(V,1);
            fprintf(fid,'%d %d\n', nv, 1);
            fprintf(fid,'%.9f %.9f %.9f\r\n', V');
            
            fclose(fid);
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
                img = imread( [dirname D(i).name] );
                writeVideo(outputVideo,img);
            end
            
            close(outputVideo);
        end
    end
    
end

