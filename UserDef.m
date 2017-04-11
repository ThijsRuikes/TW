function [user_def] = UserDef()
%V2
%% Setup variables for tracking
% These variables will control tracking output. Set 'single_frame' to 1to
% run the tracking algorithm for one frame, to test the current settings.
user_def.single_frame = 0;
% output video quality: 0-100, lower quality for smaller file size.
user_def.save_quality = 100;     
% To use the function 'Track_more' for batch processing set multiselect 
% to 1.
user_def.multiselect = 0; 

% To view the tracking live, set 'view_tracking_data' to 1. To save the
% results to a video file, set 'record'to 1. This will slowdown the
% tracking significantly.
user_def.view_tracking_data = 1; 
user_def.record = 0;  
% Do not show empty frames
user_def.skip_empty = 1;         

% In videos where mice dont move, detection of environmental objects is not
% needed.
user_def.subtrack_imobject = 1;  

% NOISE PARAMETERS
% The edges in a gapcrossing setup produce a noise pattern which is
% indistinguiable from whiskers. 'gapnoise' will detect this noise:
user_def.gapnoise = 1; 
% 'cross_gap' will let tracking pass the noise'.
user_def.cross_gap = 1;   
% 'cross_once'  will prevent the tracking from converging on noise again.
user_def.cross_once = 1;
% Gapnoise is dilated to improve interference detection,
% 'dirty_background' will prevent dilation of other noise.
user_def.dirty_background = 1;  


user_def.cut_right = 1;

% INITIALIZE
% Find begining and end of gap
user_def.find_gaplocations = 1; 
% Find direction of movement of mouse:
user_def.find_movement_direction = 1;  
% Number of frames for finding cons. obj.
user_def.n_steps = 20;                 

% NOSE TRACKING
% When electronic cables interfere with nosetracking, intensity filtering
% is the solution:
user_def.intensity_filter = 1;
% Minimum percentage the intinsity of a row has to be from the maximum
% intensity in the profile:
user_def.I_x = 0.6;
% Minimum percentage the intensity of a column has to be from the maximum
% intensity in the profile:
user_def.I_y = 0.001;

% To stop tracking if the detected nose position reaches the border of the
% frame (there are cases when the nose is out of the frame, but whiskers
% are still detected):
user_def.roi_limit = 0;    
% In some cases, the nose tracking algorithm does not work, because the
% nose is only partly in the screen. This is another way to find the nose
% position in those cases:
user_def.nose_as_lowest = 1;   
                              
% ANGLE TRACKING
% Track angle of the head. Whisker tracking itself is not dependend on
% this, but post processing will require an angle.
user_def.track_angle = 1;     
% Radius around head for tracking, larger for less error, but headstage 
% can interfere for larger radii:
user_def.angle_r = 40;         
                              

% WHISKER TRACKING
% to turn of, set 'track_whiskers' to 0;
user_def.track_whiskers = 1;  
% Traces back the whiskers from origin to snout, a bit noisy still.
user_def.track_back = 1;        

% Region size for tracking around nose(mark with red dots in video):
user_def.roi_size = 90;
% Find origins for whiskers at a distance from the snout equal to dilation
% size:
user_def.dilate_size = 25;
% Whisker origins are tracked with a peakdetection function on a intensity
% profile (around the snout). For details see the documentation of
% 'findpeaks'.  Adjusting these values will change sensitivity towards
% detection of peaks:
user_def.Peak_Width = 1;           
user_def.Peak_Distance = 1;      
user_def.Peak_Prominence = 0.02;

% Distance the iteration takes for detecting new points on a whisker:
user_def.r_step = 5; 
% A primary filter for noisy origins is detecting when the algorithm does
% not converge on new points. Setting the minimum found points on a whisker
% required to call it a whisker:
user_def.min_whisker_points = 3;

% The final criterium to stop tracking: 'stop_threshold'. This is the
% change of intensity between the found pixel and its neighbouring pixels.
% The lower the value, the longer the algorithm will converge on noise.
user_def.stop_threshold = 0.1; 

