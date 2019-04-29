function result = windowed_average_with_freq(data,points_per_win,points_per_slide)
    %This computes the moving window average of data which is channel X freq X
    %time
    %INPUTS:
        %data: channel X freq X time
        %points_per_win: the number of points per window
        %points_per_slide: the number of points that the window is slid
        %along by
    %OUTPUT:
        %result: n_channels X n_freqs X n_windows matrix of data after
        %windowed averaging is applied
    n_channels = size(data,1);
    n_freqs = size(data,2);
    total_points = size(data,3);
    points_per_overlap = points_per_win - points_per_slide;
    n_windows = floor((total_points - points_per_overlap)/points_per_slide);
    result = zeros(n_channels,n_freqs,n_windows);

    for i = 1:n_windows
        start = (i-1)*points_per_slide + 1;
        result(:,:,i) = mean(data(:,:,start:start+points_per_win-1),3);
    end
end