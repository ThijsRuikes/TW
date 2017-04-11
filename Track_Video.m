function mouse = Track_Video(Fname,Pname)
%%               %%%%%%%%% SELECT FILE FOR TRACKING%%%%%%%%%%%
user_def = UserDef;
pause('on')
if user_def.multiselect == 1
    FileName = Fname; PathName = Pname;
elseif user_def.multiselect == 0
    [FileName, PathName] = uigetfile([pwd,filesep,'*.mp4' '*.avi']);
end
addpath(PathName);
v = VideoReader(FileName);
tic
disp(['Video selected: ' FileName])

%%
%%%%%%%%% SETUP VARIABLES %%%%%%%%%%%
v.CurrentTime = 0;
n_steps = user_def.n_steps;
imsum = zeros(v.Height,v.Width);
n_frames = v.Duration*v.FrameRate;
stepsize = (n_frames/v.FrameRate)/n_steps;
mouse.nose = zeros(ceil(n_frames),2);
mouse.headvector = zeros(ceil(n_frames),2);
mouse.whiskers_x{1} = -999;
mouse.whiskers_y{1} = -999;
mouse.imobject = zeros(v.Height,v.Width);
mouse.gapnoise = zeros(v.Height,v.Width);
mouse.gaplocations = zeros(1,2);
mouse.speed = -999;
mouse.fps = -999;

%%%%%%%%% PARSE CONSTANT OBJECTS %%%%%%%%%%%
disp('Parse constant objects')
for i = 1:n_steps
    v.CurrentTime = (i-1)*stepsize;
    frame = readFrame(v);
    frame = im2double(frame);
    frame = rgb2gray(frame);
    if user_def.view_tracking_data == 1
        imagesc(frame);colormap gray;pause(0.1)
    end
    
    frame = ~imbinarize(frame,0.25);
    imsum = imsum +(frame./(n_steps));
    
end

imsum = imbinarize(imsum,0.99);
imobject=  imsum;

if user_def.view_tracking_data == 1
    imagesc(imobject);colormap gray
    title('Found objects')
    hold on
end

mouse.imobject = imobject;



%%%%%%%%% FIND MOVEMENT DIRECTION %%%%%%%%%%%

disp('Find general movement direction')
if user_def.find_movement_direction == 1
    mid_x = -999*ones(n_steps,1);
    mid_y = -999*ones(n_steps,1);
    for j = 1:n_steps
        v.CurrentTime = (j-1)*stepsize;
        frame = readFrame(v);
        frame = im2double(frame);
        frame = rgb2gray(frame);
        frame = ~imbinarize(frame,0.1);
        if user_def.subtrack_imobject == 1
            frame(imobject) = 0;
        end
        frame = bwareafilt(frame,1);
        if numel(find(frame)) > 12000
            stats = regionprops(frame,'Centroid');
            mid_x(j) = stats.Centroid(1);
            mid_y(j) = stats.Centroid(2);
        else
            mid_x(j) = NaN;
            mid_y(j) = NaN;
        end
    end
    mov_dir =[mean(diff(mid_x),'omitnan'),mean(diff(mid_y),'omitnan')];
    if user_def.view_tracking_data == 1
        title('Found objects with gap start & end, major movement direction')
        quiver(round(v.Width/2),round(v.Height/2),mov_dir(1)*10,mov_dir(2)*10,'g')
    end
    % dir(1,~) - hor, dir(2,~) - ver, dir(~,1/2) - dir
    if mov_dir(1)>mov_dir(2)
        dir(1,1) = 1;
        if mov_dir(1) >= 0
            dir(1,2) = 1;
        elseif mov_dir(1) < 0
            dir(1,2) = 2;
        end
        
    elseif mov_dir(2)>= mov_dir(1)
        dir(1,1) = 2;
        if mov_dir(2) >= 0
            dir(1,2) = 1;
        elseif mov_dir(2) < 0
            dir(1,2) =2;
        end
    end
    
    
end


mouse.imobject = imobject;


%%
%%%%%%%%% FIND GAP START AND END %%%%%%%%%%%
disp('Find gap start and end')
if user_def.find_gaplocations == 1
    if dir(1,1) == 2
        frame_x = sum(imobject,1);
        frame_y = sum(imobject,2);
        gap_y_idx = find(frame_y<50); gap_y_idx = [min(gap_y_idx),max(gap_y_idx)];
        gap_x_idx = find(frame_x<50); gap_x_idx = [min(gap_x_idx),max(gap_x_idx)];
        imobject(:,gap_x_idx(1):gap_x_idx(2)) = 0;
    elseif dir(1,1) == 1
        frame_x = sum(imobject,2);
        frame_y = sum(imobject,1);
        gap_x_idx = find(frame_y<50); gap_x_idx = [min(gap_x_idx),max(gap_x_idx)];
        gap_y_idx = find(frame_x<50); gap_y_idx = [min(gap_y_idx),max(gap_y_idx)];
        imobject(gap_x_idx(1):gap_x_idx(2),:) = 0;
    end
    
    
    if user_def.view_tracking_data == 1
        imagesc(imobject);colormap gray
        title('Found objects')
        hold on
    end
    
    %%%%%%%%%% Find local noise %%%%%%%%%%
    if user_def.gapnoise == 1
        
        imobject(gap_y_idx(1):gap_y_idx(2),:) = 0;
        gapnoise_simple = mouse.imobject; gapnoise_simple(imobject) = 0;
        gapnoise = imsum;
        gapnoise(1:gap_y_idx(1),1:gap_x_idx(1)) = 0; gapnoise(1:gap_y_idx(1),gap_x_idx(2):end) = 0;
        gapnoise(gap_y_idx(2):end, 1:gap_x_idx(1)) =0; gapnoise(gap_y_idx(2):end, gap_x_idx(2):end) =0;
        gapnoise_area = bwareafilt(gapnoise,4);
        gapnoise_noise = gapnoise;
        gapnoise_noise(gapnoise_area) = 0;
        gapnoise_area = imdilate(gapnoise_area,strel('diamond',3));
        
        if user_def.dirty_background == 1
            gapnoise_dil = gapnoise_noise+gapnoise_area;
        else
            gapnoise_dil = gapnoise_area;
        end
        mouse.gapnoise = gapnoise;
    end
    gap_locations = find(diff(frame_y)>10);
    if ~isempty(find(diff(gap_locations)>100, 1))
        gap_locations = [gap_locations(diff(gap_locations)>100) , gap_locations(find(diff(gap_locations)>100)+1)];
    elseif ~isempty(gap_locations) && gap_locations(1) >= v.Height/2
        gap_locations = [1,gap_locations(1)];
    elseif ~isempty(gap_locations) && gap_locations(1) < v.Height/2
        gap_locations = [gap_locations(1), v.Height];
    else
        gap_locations = [1,v.Height];
    end
    if user_def.view_tracking_data == 1
        title('Found objects with gap start and end')
        if dir(1,1) == 2
            scatter([round(v.Width/2) round(v.Width/2)],gap_locations,50,'r','filled')
        elseif dir(1,1) ==  1
            scatter(gap_locations,[round(v.Width/2) round(v.Width/2)],50,'r','filled')
        end
    end
end

mouse.gaplocations = gap_locations;


pause(0.1)
%%
%%%%%%%%% TRACK VIDEO %%%%%%%%%%%
ref = [ 0 1];
if user_def.record == 1
    disp('Start nose tracking and recording')
else
    disp('Start nose tracking')
end

v.CurrentTime = 0;
pause('on')
frame_nr =0;

user_def = UserDef;
if user_def.record == 1
    video = VideoWriter([FileName(1:end-4) '_tracked']);
    video.FrameRate = v.FrameRate;
    video.Quality = user_def.save_quality;
    open(video)
end
%
%%%%%%%%%%%%% LOOP OVER FRAMES %%%%%%%%%%%%%

while hasFrame(v)
    
    %%%%%%%%% lOAD AND PROCESS FRAME %%%%%%%%%%%
    if user_def.single_frame == 1
        t = input('Enter time of frame :');
        v.CurrentTime = t;
    end
        
    frame_nr = frame_nr+1;
    if ceil((frame_nr/round((n_frames/20)))) == (frame_nr/round((n_frames/20)))
        disp( [num2str(round((frame_nr/n_frames)*100)) ' % Processed'])
    end
    frame_load = readFrame(v);
    frame = im2double(frame_load);
    frame_gray = rgb2gray(frame);
    frame_sill = ~imbinarize(frame_gray,0.1);
    if user_def.subtrack_imobject == 1
        frame_sill(imobject) = 0;
    end
    if user_def.cut_right == 1
        if numel(find(frame_sill)) < 30000
            frame_sill(:,v.Width/2:end) =0;
        else
            frame_sill(gapnoise_simple) = 0;
        end
    end
    
    frame_sill_nose = bwareafilt(frame_sill,1);
    frame_sill_dilated_full = imdilate(frame_sill,strel('diamond',user_def.dilate_size));
    frame_sill_dilated_full = bwareafilt(frame_sill_dilated_full,1);
    frame_sill_dilated = edge(frame_sill_dilated_full);
    
    if user_def.intensity_filter == 1
        if dir(1,1) == 2
            frame_x  = sum(frame_sill_nose,1);
            frame_x_idx = frame_x <= (max(frame_x)*user_def.I_x);
            frame_y = sum(frame_sill_nose,2);
            frame_y_idx = frame_y <= (max(frame_y)*user_def.I_y);
            frame_sill_nose(:,frame_x_idx) =0;
            frame_sill_nose(frame_y_idx,:) = 0;
            frame_sill_nose =  bwareafilt(frame_sill_nose,1);
        elseif dir(1,1) == 1
            frame_y = sum(frame_sill_nose,2);
            frame_y_idx = frame_y <= (max(frame_y)*user_def.I_y);
            frame_x = sum(frame_sill_nose,1);
            frame_x_idx = frame_x <= (max(frame_x)*user_def.I_x);
            frame_sill_nose(:,frame_x_idx) = 0;
            frame_sill_nose(frame_y_idx,:) =0;
            frame_sill_nose =  bwareafilt(frame_sill_nose,1);
        end
    end
    
    imshow(frame_sill_nose)
    
    %%%%%%%%% DETECT MOUSE PRESENCE %%%%%%%%%%%
    
    
    if numel(find(frame_sill_nose))>1300 %%%%%%%%%%%
        
        %%%%%%%%% TRACK NOSE POSITION  %%%%%%%%%%%
        
        
        stats = regionprops(frame_sill_nose,'Centroid');
        base(2) = round(stats.Centroid(1));
        base(1) = round(stats.Centroid(2));
        if dir == 1
            frame_sill_nose(1:base(1),:) = 0;
        elseif dir == 2
            frame_sill_nose(base(1):v.Height,:) = 0;
        elseif dir == 3
            frame_sill_nose(:,1:base(2)) = 0;
        elseif dir == 4
            frame_sill_nose(:,base(2):v.Width) = 0;
        end
        
        edge_idx = find(edge(frame_sill_nose));
        e_idx = -999*ones(length(edge_idx),2);
        [e_idx(:,2),e_idx(:,1)] = ind2sub(size(frame_sill_nose),edge_idx);
        d1 = (e_idx(:,1)- base(2)); d1 = d1.^2;
        d2 = (e_idx(:,2) - base(1)); d2 = d2.^2;
        distance = sqrt(d1+d2);
        nose = e_idx(distance == max(distance),:);
        nose = nose(1,:);
        
        % MODULE
        if user_def.nose_as_lowest == 1
            if dir(1,1) == 2 && dir(1,2) == 2 && find(sum(frame_sill_nose,2),1,'first') < round(v.Height*0.4)
                nose(2) =  find(sum(frame_sill_nose,2),1,'first');
                nose(1) = round(mean(find(frame_sill_nose(nose(2),:))));
            elseif dir(1,1) == 2 && dir(1,2) == 1 && find(sum(frame_sill_nose,2),1,'last') < round(v.Height*0.6)
                nose(2) =  find(sum(frame_sill_nose,2),1,'last');
                nose(1) = round(mean(find(frame_sill_nose(nose(2),:))));
            elseif dir(1,1) == 1 && dir(1,2) == 1 && find(sum(frame_sill_nose,1),1,'last') < round(v.Width*0.6)
                nose(1) = find(sum(frame_sill_nose),1,'last');
                nose(2) = round(mean(find(frame_sill_nose(:,nose(1)))));
            elseif dir(1,1) == 1 && dir(1,2) == 2 && find(sum(frame_sill_nose,1),1,'first') < round(v.Width*0.4)
                nose(1) = find(sum(frame_sill_nose),1,'first');
                nose(2) = round(mean(find(frame_sill_nose(:,nose(1)))));
            end
        end
        
        
        if dir(1,1) == 2 && dir(1,2) == 1 && nose(2) >= v.Height
            nose(1) = round(mean(find(frame_sill_nose(v.Height-10,:))));
        elseif dir(1,1) == 2 && dir(1,2) == 2 && nose(2)<= 2
            nose(1) = round(mean(find(frame_sill_nose(1+10,:))));
        elseif dir(1,1) == 1 && dir(1,2) == 1 && nose(1) >= v.Width
            nose(1) = round(mean(find(frame_sill_nose(:,v.Width-10))));
        elseif dir(1,1) == 1 && dir(1,2) == 2 && nose(1) <= 2
            nose(1) = round(mean(find(frame_sill_nose(1+10))));
        end
        
        
        
        %%%%%%%%% TRACK ANGLE OF THE HEAD %%%%%%%%%%%
        
        % MODULE
        if user_def.track_angle == 1
            headvector = [0 0];
            nose_edge = edge(frame_sill_nose);
            
            theta = 1:1:360;
            c_y = round(nose(2) +  user_def.angle_r*sind(theta)  );
            c_x = round( nose(1) +  user_def.angle_r*cosd(theta)  );
            idd = find(c_x > 1 & c_x < v.Width & c_y>1 & c_y<v.Height);
            c_ind = sub2ind(size(nose_edge),c_y(idd),c_x(idd));
            
            points = find(nose_edge(c_ind));
            if ~isempty(points)
                p = [points(1),points(find(diff(points)>50)+1)];
                if size(p,2)>1
                    [a,b] = ind2sub(size(nose_edge),c_ind(p));
                    points_sub = [a',b'];
                    if dir(1,1) == 1
                        if b(2)>b(1)
                            points_sub = flip(points_sub,1);
                        end
                    elseif dir(1,1) == 2
                        if a(2)<a(1)
                            points_sub = flip(points_sub,1);
                        end
                    end
                    point_vec = points_sub(1,:)-points_sub(2,:);
                    headvector = [-point_vec(1) point_vec(2)];
                    headvector = -headvector./sqrt(headvector(1)^2+headvector(2)^2);
                    mouse.headvector(frame_nr,:) = headvector;
                end
            end
            
        end
        
        
        % MODULE
        if user_def.track_whiskers == 1
            if user_def.roi_limit == 1
                if dir(1,1) == 2 && dir(1,2) == 1 && nose(2) > v.Height-5
                    figure(2)
                    clf
                    imagesc(frame_gray);colormap gray
                    drawnow
                    continue
                elseif dir(1,1) == 2 && dir(1,2) == 2 && nose(2) < 1+5
                    figure(2)
                    clf
                    imagesc(frame_gray);colormap gray
                    drawnow
                    continue
                elseif dir(1,1) == 1 && dir(1,1) == 1 && nose(1) > v.Width-5
                    figure(2);clf;imagesc(frame_gray);colormap gray;drawnow
                    continue
                elseif dir(1,1) == 1 && dir(1,2) == 2 && nose(1) < 1+5
                    figure(2);clf;imagesc(frame_gray);colormap gray;drawnow
                    continue
                end
            end
            
            
            %%%%%%%%% FIND WHISKER ORIGINS %%%%%%%%%%%
            r = user_def.roi_size;
            roi(1,1) = nose(1)-r; if roi(1,1)<1;roi(1,1) = 1;end
            roi(1,2) = nose(2)-r; if roi(1,2)<1;roi(1,2) = 1;end
            roi(2,1) = nose(1)+r; if roi(2,1)>v.Width;roi(2,1) = v.Width;end
            roi(2,2) = nose(2)+r; if roi(2,2)>v.Height;roi(2,2) = v.Height;end
            sill_roi = frame_sill_dilated; sill = sill_roi;
            sill_roi(roi(1,2):roi(2,2),roi(1,1):roi(2,1)) = 0;
            sill_roi = sill - sill_roi;
            sill_roi(gapnoise) = 0;
            sill_roi(imobject) = 0;
            w_roi = find(sill_roi);
            roi(3,1) = roi(1,1); roi(3,2) = roi(2,2);
            roi(4,1) = roi(2,1); roi(4,2) = roi(1,2);
            
            if length(w_roi)> 5
                frame_gray(imobject) = 0;
                i_profile = 1 - frame_gray(w_roi);
                i_profile = i_profile - min(i_profile);
                [~,locs]  = findpeaks(i_profile,'MinPeakDistance',user_def.Peak_Distance,'MinPeakWidth',user_def.Peak_Width,'MinPeakProminence',user_def.Peak_Prominence);
                start = w_roi(locs);
                wstart = -999*ones(length(start),2);
                [wstart(:,1),wstart(:,2)] = ind2sub(size(frame_gray),start);
            end
            
            %%%%%%%%% TRACK WHISKERS FORWARD %%%%%%%%%%%
            
            whiskers_x = ones(length(wstart),1)*-999;
            whiskers_y = ones(length(wstart),1)*-999;
            for w_nr = 1:size(wstart,1)
                
                p_nr = 1;
                frame_gray2 = frame_gray;
                theta = 1:1:360;
                w_points= [-999 -999];
                w_points(1,:) = wstart(w_nr,:);
                c_y     = round( w_points(1,2) +  user_def.r_step * sind(theta) );
                c_x     = round( w_points(1,1) +  user_def.r_step * cosd(theta) );
                ind = find(c_y > 1 & c_y < v.Width & c_x>1 & c_x < v.Height);
                c_ind = sub2ind(size(frame_load),c_x(ind),c_y(ind));
                c_ind = c_ind(gapnoise_dil(c_ind)==0);
                c_ind = c_ind(frame_sill_dilated_full(c_ind) == 0);
                
                i_profile = frame_gray2(c_ind);
                i_min = min(i_profile);
                point_idx = frame_gray2(c_ind)==i_min;
                if isempty(point_idx);continue;end
                point_idx = round(mean(find(point_idx)));
                point_new = c_ind(point_idx);
                i_profile = i_profile-i_min;
                i_sum(p_nr) = sum(i_profile); %#ok<AGROW>
                p_nr = 2;
                countdown = -1;
                crossed  =0;
                while i_sum(p_nr-1) >= user_def.stop_threshold && p_nr<50
                    
                    [w_points(p_nr,1),w_points(p_nr,2)] = ind2sub(size(frame_gray),point_new);
                    w_direction = w_points(p_nr,:)-w_points(p_nr-1,:);
                    if countdown > 0
                        w_direction = w_points(p_nr_first,:)-w_points(p_nr_first-1,:);
                        w_points(p_nr,:) = w_points(p_nr-1,:)+ w_direction;
                        if w_points(p_nr,1)>v.Height || w_points(p_nr,2) < v.Width
                            break
                        end
                    end
                    
                    if user_def.cross_gap == 1
                        if gapnoise_dil(w_points(p_nr,1),w_points(p_nr,2)) == 1 || countdown >=0
                            if p_nr <=5
                                % break
                            end
                            if gapnoise_dil(w_points(p_nr,1),w_points(p_nr,2)) == 0
                                countdown = -1;
                                if user_def.cross_once == 1
                                    crossed = 1;
                                end
                            elseif countdown == -1 && gapnoise_dil(w_points(p_nr,1),w_points(p_nr,2)) == 1
                                w_points(p_nr,:) = w_points(p_nr-1,:) +  (w_points(p_nr-1,:)-w_points(p_nr-2,:));
                                countdown = 7;
                                p_nr_first = p_nr-1;
                            elseif countdown > 0 &&  gapnoise_dil(w_points(p_nr,1),w_points(p_nr,2)) == 1
                                countdown = countdown - 1;
                            elseif countdown == 0 &&  gapnoise_dil(w_points(p_nr,1),w_points(p_nr,2)) == 1
                                % break
                            end
                        end
                        
                    elseif user_def.cross_gap == 0
                        if gapnoise_dil(w_points(p_nr,1),w_points(p_nr,2)) == 1
                            % break
                        end
                    end
                    
                    if imobject(w_points(p_nr,1),w_points(p_nr,2)) == 1
                        % break
                    end
                    angle  = atan2d(w_direction(1),w_direction(2)) - atan2d(ref(1),ref(2));
                    theta = (30:10:150) - angle;
                    c_y = round( w_points(p_nr,2) +  user_def.r_step*sind(theta)  );
                    c_x = round( w_points(p_nr,1) +  user_def.r_step*cosd(theta)  );
                    ind = find(c_y > 1 & c_y < v.Width & c_x>1 & c_x < v.Height);
                    c_ind = sub2ind(size(frame_load),c_x(ind),c_y(ind));
                    if crossed == 1
                        c_ind = c_ind(gapnoise_dil(c_ind)==0);
                    end
                    c_ind = c_ind(frame_sill_dilated_full(c_ind) == 0);
                    i_profile = frame_gray2(c_ind);
                    i_min = min(i_profile);
                    point_idx = frame_gray2(c_ind)==i_min;
                    
                    if isempty(point_idx);break;end
                    point_idx = round(mean(find(point_idx)));
                    point_new = c_ind(point_idx);
                    i_profile = i_profile-i_min;
                    i_sum(p_nr) = sum(i_profile);
                    p_nr = p_nr +1;
                end
                whiskers_x(w_nr,1:size(w_points,1)) = w_points(:,1);
                whiskers_y(w_nr,1:size(w_points,1)) = w_points(:,2);
                if user_def.track_back == 0
                    whisk_temp_x = whiskers_x;
                    whisk_temp_y = whiskers_y;
                end
            end
        end
        
        %%%%%%%%% TRACK WHISKERS BACK %%%%%%%%%%%
        if user_def.track_back == 1
            whiskers_x_back = ones(size(whiskers_x,1),1)*-999;
            whiskers_y_back = ones(size(whiskers_x,1),1)*-999;
            for w_nr = 1:size(wstart,1)
                
                p_nr = 2;
                w_points = [-999 -999];
                if size(whiskers_x,2)<2 | find(whiskers_x(w_nr,1:2))== 0 %#ok<OR2>
                    break
                end
                w_points(1,:) = [whiskers_x(w_nr,2),whiskers_y(w_nr,2)];
                w_points(2,:) = [whiskers_x(w_nr,1),whiskers_y(w_nr,1)];
                w_direction =   w_points(2,:)-w_points(1,:);
                angle  = atan2d(w_direction(1),w_direction(2)) - atan2d(ref(1),ref(2));
                theta = (30:10:150) - angle;
                c_y = round( w_points(p_nr,2) +  user_def.r_step*sind(theta)  );
                c_x = round( w_points(p_nr,1) +  user_def.r_step*cosd(theta)  );
                ind = find(c_y > 1 & c_y < v.Width & c_x>1 & c_x < v.Height);
                c_ind = sub2ind(size(frame_load),c_x(ind),c_y(ind));
                i_profile = frame_gray2(c_ind);
                i_min = min(i_profile);
                point_idx = frame_gray2(c_ind)==i_min;
                
                if isempty(point_idx);break;end
                point_idx = round(mean(find(point_idx)));
                point_new = c_ind(point_idx);
                i_profile = i_profile-i_min;
                i_sum(p_nr) = sum(i_profile);
                p_nr = p_nr+1;
                while i_sum(p_nr-1) >= user_def.stop_threshold
                    
                    [w_points(p_nr,1),w_points(p_nr,2)] = ind2sub(size(frame_gray),point_new);
                    
                    
                    if gapnoise_dil(w_points(p_nr,1),w_points(p_nr,2)) == 1
                        break
                    end
                    
                    if imobject(w_points(p_nr,1),w_points(p_nr,2)) == 1
                        break
                    end
                    
                    if frame_sill(w_points(p_nr,1),w_points(p_nr,2)) == 1
                        
                        break
                    end
                    w_direction = w_points(p_nr,:)-w_points(p_nr-1,:);
                    angle  = atan2d(w_direction(1),w_direction(2)) - atan2d(ref(1),ref(2));
                    theta = (30:10:150) - angle;
                    c_y = round( w_points(p_nr,2) +  user_def.r_step*sind(theta)  );
                    c_x = round( w_points(p_nr,1) +  user_def.r_step*cosd(theta)  );
                    ind = find(c_y > 1 & c_y < v.Width & c_x>1 & c_x < v.Height);
                    c_ind = sub2ind(size(frame_load),c_x(ind),c_y(ind));
                    i_profile = frame_gray2(c_ind);
                    i_min = min(i_profile);
                    point_idx = frame_gray2(c_ind)==i_min;
                    
                    if isempty(point_idx);break;end
                    point_idx = round(mean(find(point_idx)));
                    point_new = c_ind(point_idx);
                    i_profile = i_profile-i_min;
                    i_sum(p_nr) = sum(i_profile);
                    p_nr = p_nr +1;
                end
                
                
                whiskers_x_back(w_nr,1:size(w_points,1)) = w_points(:,1);
                whiskers_y_back(w_nr,1:size(w_points,1)) = w_points(:,2);
                
                
            end
            n_points_found = find(sum(whiskers_x_back,1), 1, 'last' );
            n_points_tip = size(whiskers_x,2);
            
            if ~isempty(n_points_found)
                l = n_points_found+ n_points_tip;
                tag = 1;
            elseif ~isempty(n_points_tip)
                l = n_points_tip;
                tag = 2;
            else
                %break %?
            end
            if tag == 1
                whisk_temp_x = zeros(size(whiskers_x,1),l);
                whisk_temp_x(:,1:n_points_found) = flip(whiskers_x_back);
                whisk_temp_x(:,n_points_found+1:end) = whiskers_x;
                whisk_temp_y = zeros(size(whiskers_y,1),l);
                whisk_temp_y(:,1:n_points_found) = flip(whiskers_y_back);
                whisk_temp_y(:,n_points_found+1:end) = whiskers_y;
            elseif tag == 2
                whisk_temp_x = zeros(size(whiskers_x,1),l);
                whisk_temp_y = zeros(size(whiskers_y,1),l);
                whisk_temp_x(:,1:end) = whiskers_x;
                whisk_temp_y(:,1:end) = whiskers_y;
            else
                %break
            end
        end
        
        if user_def.track_whiskers == 1
            whiskers_x = -999*ones(1,size(whisk_temp_x,2));
            whiskers_y = -999*ones(1,size(whisk_temp_x,2));
            w_tick = 1;
            for idd = 1:size(whisk_temp_x,1)
                if numel(find(whisk_temp_x(idd,:))) > user_def.min_whisker_points
                    whiskers_x(w_tick,:) = whisk_temp_x(idd,:);
                    whiskers_y(w_tick,:) = whisk_temp_y(idd,:);
                    w_tick = w_tick + 1;
                end
            end
            
            mouse.whiskers_y{frame_nr} = whiskers_y;
            mouse.whiskers_x{frame_nr} = whiskers_x;
        end
        
        
        
        %%%%%%%%% SHOW RESULTS %%%%%%%%%%%
        
        
        if user_def.view_tracking_data == 1 || user_def.record == 1
            figure(2);clf; imagesc(frame_load);hold on;scatter(nose(1),nose(2),'b','filled')
            if user_def.track_angle == 1
                quiver(nose(1),nose(2),headvector(1)*30,headvector(2)*30,1,'m')
            end
            if user_def.track_whiskers == 1
                if exist('wstart','var') && sum(sum(wstart))>0
                    for kk = 1:size(whiskers_y,1)
                        scatter(whiskers_y(kk,:),whiskers_x(kk,:),2,'g','filled')
                        scatter(roi(:,1),roi(:,2),10,'r','filled','square')
                    end
                end
            end
        end
    elseif user_def.skip_empty == 1
        continue
    elseif user_def.view_tracking_data == 1 || user_def.record == 1
        figure(2);clf;imagesc(frame_load)
    end
    drawnow limitrate
    if user_def.record == 1
        frame_video = getframe;
        writeVideo(video,frame_video);
    end
    if exist('nose','var')
        mouse.nose(frame_nr,:) = nose(1,:);
    end
    
    if user_def.single_frame == 1
        cc = input('Again? (0 if no, 1 if yes): ');
        if cc == 1
            continue
        elseif cc == 0
            break
        end               
    end
    
end



% Wrapping things up
if user_def.record == 1
    close(video)
end

mouse.speed = toc;
mouse.fps = ceil(n_frames)/toc;
if user_def.single_frame == 0
    close all
    disp(['Finished, time elapsed: ' num2str(toc)])
    disp([num2str(ceil(n_frames)) ' frames analysed, at ' num2str(ceil(n_frames)/toc) ' frames per second'])
    save([FileName(end-4) '_tracked'], 'mouse' )
end
