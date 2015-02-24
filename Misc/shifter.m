function result = shifter(data,shift)
% function result = shifter(data,shift)
%
% illustration of the time shit in time
% 

    close
    subplot 211, plot(abs(data)) , 
	subplot 212, plot(angle(data)),
    
    % make sure the input is a row
    N=max(size(data));
    data=reshape(data,1,N);
     
	% make a vector to go from -pi tp pi
    w = [-N/2+1 : N/2] *(2*pi/N);
    w = fftshift(w);
    
    fd= fft(data);
    
    fd2 = fd  .* exp(-j.*w*shift);
    	   
    result = real(ifft(fd2));
    
	subplot 211, hold on, plot(abs(result),'r'), grid
    subplot 212, hold on, plot(angle(result),'r')
    hold on, plot(w,'y')

    return
