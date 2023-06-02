function [stp,Hend,Gend,nfev,info] = armijo(n,m,x,d,f0,Jf,theta,ftol,A,b)

% Counters

nfev = 0;

% Compute g0

g0 = Jf * d;

% Define ftest

ftest = - ftol * abs(theta);	

% Test the sufficient descent condition (sdc) at stp = 1

stp = 1;

sdc = true;     
for i = 1:m
		
	[h] = evalh(n,x+stp*d,i);
    [g,flag] = evalg(n,x+stp*d,i,A,b);
    
    f = g + h;
	nfev = nfev + 1;
	
	Hend(i) = h;
    Gend(i) = g;
	
	if ( f > f0(i) + ftest * stp )
		sdc = false;
		iA  = i;
		fiA = f;
		break
    end			
end

if ( sdc )
	info = 0;
	return
end

[stp,Hend,Gend,nfevbt,info] = backtrackingMO(stp,n,m,x,d,f0,g0,iA,fiA,ftest,theta,A,b);
nfev = nfev + nfevbt;

end

%***********************************************************************
%***********************************************************************	

function [stp,Hend,Gend,nfev,info] = backtrackingMO(stp,n,m,x,d,f0,g0,iA,fiA,ftest,theta,A,b)	

% Parameters

stpmin = 10^(-15);
sigma1 = 0.5d-1; 
sigma2 = 9.5d-1;

% Counters

outiter = 0;
nfev    = 0;
	
%-------------------------------------------------------------------
%     Main loop
%-------------------------------------------------------------------    

while (1)

	% Test the vector Armijo condition
	
	if ( outiter == 0 )
		sdc = false;
		ind = iA;
		f   = fiA;
	elseif ( infoBT == 0 )
		sdc = true;
		for i = 1:m				

			if ( i == ind ) 
                continue
            end
		
			[h] = evalh(n,x+stp*d,i);
            [g,flag] = evalg(n,x+stp*d,i,A,b);
            
            f = g + h;
	        nfev = nfev + 1;
	        
	        Hend(i) = h;
            Gend(i) = g;

			if ( f > f0(i) + ftest * stp )
				sdc = false;
				ind = i;
				break
            end
        end
		
    end
	
	% Finish backtracking with the current point
	
	if ( sdc )
		info = 0;
		return
    end
	
	% Test if stp is too small
	
	if ( stp <= stpmin )
		stp = stpmin;
		info = 1;
		
		for i = 1:m				
			
			if ( i == ind ) 
                continue
            end
			
			[h] = evalh(n,x+stp*d,i);
            [g,flag] = evalg(n,x+stp*d,i,A,b);
            
            f = g + h;
	        nfev = nfev + 1;
            
            Hend(i) = h;
            Gend(i) = g;
        end
		return
    end

	outiter = outiter + 1;
	
	% Compute new trial stepsize based on f_ind
	
	while (1)
	
		% Test Armijo condition for f_ind
		
		if ( f <= f0(ind) + ftest * stp )
			infoBT = 0;
			break
		end 
		
		if ( stp <= stpmin )
            infoBT = 1;
            break
        end
		
	if ( g0(ind) < 0 )

            stpt = ( (g0(ind) / ( (f0(ind)-f) / stp + g0(ind) ) ) / 2 ) * stp;

            if ( stpt >= sigma1 * stp && stpt <= sigma2 * stp )
                stp = stpt;
            else
                stp = stp / 2;
            end
        else
            stp = stp / 2;
        end
		
	[h] = evalh(n,x+stp*d,ind);
        [g,flag] = evalg(n,x+stp*d,ind,A,b);
        
        f = g + h;
        nfev = nfev + 1;
	        
    end
	
	Hend(ind) = h;
    Gend(ind) = g;

end

%--------------------------------------------------------------------- 
%     End of main loop
%---------------------------------------------------------------------

end
