function [n,m,l,u,x] = inip;

global problem

% ----------------------------------------------------------------------

	% AP1
	% Example 1 of "A modified Quasi-Newton method for vector optimization problem"

	if ( strcmp(problem,'AP1' )  ) 
	
		% Number of variables
		
		n = 2;
		
		% Number of objectives
		
		m = 3;

		% Box constraints
	
		l(1:n,1) = - 1.0d1;
		u(1:n,1) =   1.0d1;
		
		% Initial point
		
        for i = 1:n
            x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
        end 
		
		return
    end
	
% ----------------------------------------------------------------------

	% AP2
	% Example 2 of "A modified Quasi-Newton method for vector optimization problem"

	if ( strcmp(problem,'AP2' )  ) 
	
		% Number of variables
		
		n = 1;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints
	
		l(1:n,1) = - 1.0d2;
		u(1:n,1) =   1.0d2;
		
		% Initial point
		
		a = - 1.0d2;
		b =   1.0d2;
		
        for i = 1:n
            x(i,1) = a + ( b - a ) *  rand;
        end 

        return
    end

% ----------------------------------------------------------------------

	% AP3
	% Example 3 of "A modified Quasi-Newton method for vector optimization problem"

	if ( strcmp(problem,'AP3' )  ) 
	
		% Number of variables
		
		n = 2;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints
		
		l(1:n,1) = - 1.0d2;
		u(1:n,1) =   1.0d2;
		
		% Initial point
		
		a = - 1.0d2;
		b =   1.0d2;
		
        for i = 1:n
            x(i,1) = a + ( b - a ) *  rand;
        end 
		
        return
    end

% ----------------------------------------------------------------------

	% AP4
	% Exemple 4 of "A modified Quasi-Newton method for vector optimization problem"
	
	if ( strcmp(problem,'AP4' )  ) 
	
		% Number of variables
		
		n = 3;
		
		% Number of objectives
		
		m = 3;
		
		% Box constraints
	
		l(1:n,1) = - 1.0d1;
		u(1:n,1) =   1.0d1;
		
		% Initial point
		
		a = - 1.0d1;
		b =   1.0d1;
		
        for i = 1:n
            x(i,1) = a + ( b - a ) *  rand;
        end 
		
		return
    end
	 
	
% ----------------------------------------------------------------------

	%  BK1
	%  A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'BK1' )  ) 
	
		% Number of variables
		
		n = 2;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints
	
		l(1:n,1) = - 5.0d0;
		u(1:n,1) =   1.0d1;
		
		% Initial point
		
		a = - 5.0d0;
		b =   1.0d1;
		
		 for i = 1:n
			x(i,1) = a + ( b - a ) *  rand;
		 end 

		return
    end
	 

% ----------------------------------------------------------------------

	%   DD1
	%   I. Das and J. E. Dennis. Normal-boundary intersection: A new method for generating
	%		the Pareto surface in nonlinear multicriteria optimization problems. SIAM J. Optim.,
	%		8(3):631???657, 1998.
	
	if ( strcmp(problem,'DD1' )  ) 

		% Number of variables

		n = 5;
		
		% Number of objectives
		
		m = 2;

		% Box constraints
	
		l(1:n,1) = - 2.0d1;
		u(1:n,1) =   2.0d1;
		
		% Initial point

		a = - 1.0d0;
		b =   1.0d0;
		
		 for i = 1:n
			x(i,1) = a + ( b - a ) *  rand;
		 end 
		
		x = 1.0d1 * rand * x / norm(x);
		
		return
    end
	 			
	
% ----------------------------------------------------------------------

	%   DGO1
	%   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'DGO1' )  ) 
	
		% Number of variables

		n = 1;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints
	
		l(1:n,1) = - 1.0d1;
		u(1:n,1) =   1.3d1;
		
		% Initial point
		
		a = - 1.0d1;
		b =   1.3d1;
		
		 for i = 1:n
			x(i,1) = a + ( b - a ) *  rand;
		 end 
		
		return
    end
	
% ----------------------------------------------------------------------

	%   DGO2
	%   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'DGO2' )  ) 

		% Number of variables

		n = 1;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints
	
		l(1:n,1) = - 9.0d0;;
		u(1:n,1) =   9.0d0
		
		% Initial point
		
		a = - 9.0d0;
		b =   9.0d0;
		
		 for i = 1:n
			x(i,1) = a + ( b - a ) *  rand;
		 end 
	
		return

    end
	 		
	
% ----------------------------------------------------------------------


	%   DTLZ1
	%   Scalable test problems for evolutionary multi-objective optimization
	
	if ( strcmp(problem,'DTLZ1' )  ) 
	
		% Number of objectives
		
		m = 3;
		
		dimm = m;

		% Number of variables
		
		dimk = 5;

		n = dimk + m - 1;
		
%		dimk = n - m + 1;
		
		% Box constraints
	
		l(1:n,1) = 0.0d0;
		u(1:n,1) = 1.0d0;
		
		% Initial point
		
		a = 0.0d0;
		b = 1.0d0;
		
		 for i = 1:n
			x(i,1) = a + ( b - a ) *  rand;
		 end 
		 
         return

    end
	 	
	
% ----------------------------------------------------------------------


	%   DTLZ2
	%   Scalable test problems for evolutionary multi-objective optimization
	
	if ( strcmp(problem,'DTLZ2' )  ) 
	
		% Number of objectives
		
		m = 3;
		
		dimm = m;

		% Number of variables
		
		dimk = 5;

		n = dimk + m - 1;
		
%		dimk = n - m + 1
		
		% Box constraints
	
		l(1:n,1) = 0.0d0;
		u(1:n,1) = 1.0d0;
		
		% Initial point
		
		a = 0.0d0;
		b = 1.0d0;
		
		 for i = 1:n
			x(i,1) = a + ( b - a ) *  rand;
		 end 

		return

    end
	 
	
% ----------------------------------------------------------------------


	%   DTLZ3
	%   Scalable test problems for evolutionary multi-objective optimization
	
	if ( strcmp(problem,'DTLZ3' )  ) 
	
		% Number of objectives
		
		m = 3;
		
		dimm = m;

		% Number of variables
		
		dimk = 5;

		n = dimk + m - 1;
		
%		dimk = n - m + 1
		
		% Box constraints
	
		l(1:n,1) = 0.0d0;
		u(1:n,1) = 1.0d0;
		
		% Initial point
		
		a = 0.0d0;
		b = 1.0d0;
		
		 for i = 1:n
			x(i,1) = a + ( b - a ) *  rand;
         end 

		return

    end
	 	

% ----------------------------------------------------------------------


	%   DTLZ4
	%   Scalable test problems for evolutionary multi-objective optimization
	
	if ( strcmp(problem,'DTLZ4' )  ) 
	
		% Number of objectives
		
		m = 3;
		
		dimm = m;

		% Number of variables
		
		dimk = 5;

		n = dimk + m - 1;
		
%		dimk = n - m + 1
		
		% Define parameter alpha
		
		alpha = 2.0d0;
		
		% Box constraints
	
		l(1:n,1) = 0.0d0;
		u(1:n,1) = 1.0d0;
		
		% Initial point
		
		a = 0.0d0;
		b = 1.0d0;
		
		 for i = 1:n
			x(i,1) = a + ( b - a ) *  rand;
         end
			
		return

    end

% ----------------------------------------------------------------------		
	

	%   FA1
	%   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'FA1' )  ) 
	
		% Number of variables

		n = 3;
		
		% Number of objectives
		
		m = 3;

		% Box constraints
	
		l(1:n,1) = 0.01d0;
		u(1:n,1) = 1.0d0;
		
		% Initial point
		
		a = 0.01d0;
		b = 1.0d0;
		
		 for i = 1:n
			x(i,1) = a + ( b - a ) *  rand;
		 end 
		
		return

    end
	 			
	
% ----------------------------------------------------------------------

	%   Far1
	%   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'Far1' )  ) 
	
		% Number of variables

		n = 2;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints
	
		l(1:n,1) = - 1.0d0;
		u(1:n,1) =   1.0d0;
		
		% Initial point

		a = - 1.0d0;
		b =   1.0d0;
		
		 for i = 1:n
			x(i,1) = a + ( b - a ) *  rand;
         end	 
		
		return

    end
	 	
	
% ----------------------------------------------------------------------
	
	% FDS
	% NEWTONS METHOD FOR MULTIOBJECTIVE OPTIMIZATION
	
	if ( strcmp(problem,'FDS' )  ) 
	
		% Number of variables
		
		n = 5;
		
		% Number of objectives
		
		m = 3;

		% Box constraints
	
		l(1:n,1) = - 2.0d0;
		u(1:n,1) =   2.0d0;
		
		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return
    end
	
% ----------------------------------------------------------------------
	
	% FF1
	% C. M. Fonseca and P. J. Fleming, ???An overview of evolutionary algorithms in multiobjective optimization
	
	if ( strcmp(problem,'FF1' )  ) 
	
		% Number of variables
		
		n = 2;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = - 1.0d0;
		u(1:n,1) =   1.0d0;
		
		% Initial point
		
		a = - 1.0d0;
		b =   1.0d0;
		
		 for i = 1:n
			x(i,1) = a + ( b - a ) *  rand;
		 end 
		
		return

    end
	 		
				
	% ----------------------------------------------------------------------	

	%   Hil1
	%   Generalized Homotopy Approach to Multiobjective Optimization
	
	if ( strcmp(problem,'Hil1' )  ) 
	
		% Number of variables

		n = 2;

		% Number of objectives

		m = 2;	
		
		% Box constraints

		l(1:n,1) = 0.0d0;
		u(1:n,1) = 1.0d0;

		% Initial point

		a =  0.0d0;
		b =  1.0d0;

		 for i = 1:n
			x(i,1) = a + ( b - a ) *  rand;
		 end 
		
		return
    end
	
% ----------------------------------------------------------------------	

	%   IKK1
	%   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit

	
	if ( strcmp(problem,'IKK1' )  ) 
	
		% Number of variables

		n = 2;

		% Number of objectives

		m = 3;
		
		% Box constraints

		l(1:n,1) = - 5.0d1;
		u(1:n,1) =   5.0d1;

		% Initial point

		a = - 5.0d1;
		b =   5.0d1;

		 for i = 1:n
			x(i,1) = a + ( b - a ) *  rand;
		 end 
		
		return

    end
	 	
	
% ----------------------------------------------------------------------	

	%   IM1
	%   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit

	
	if ( strcmp(problem,'IM1' )  ) 
	
		% Number of variables

		n = 2;

		% Number of objectives

		m = 2;
		
		% Box constraints

		l(1) = 1.0d0;
		u(1) = 4.0d0;
		
		l(2) = 1.0d0;
		u(2) = 2.0d0;
        
        l = l';
        u = u';

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
         end

		return
    end
	 

% ----------------------------------------------------------------------	

	% JOS1
	% Dynamic Weighted Aggregation for Evolutionary Multi-Objetive Optimization: Why  fores It Work and How?
	
	if ( strcmp(problem,'JOS1' )  ) 

		% Number of variables

		n = 2;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = - 1.0d2;
		u(1:n,1) =   1.0d2;
		
		% Initial point
		
		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return
    end
	
% ----------------------------------------------------------------------	

	% JOS4
	% Dynamic Weighted Aggregation for Evolutionary Multi-Objetive Optimization: Why  fores It Work and How?
	% See also: NEWTON???S METHOD FOR MULTIOBJECTIVE OPTIMIZATION
	
	if ( strcmp(problem,'JOS4' )  ) 

		% Number of variables

		n = 20;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = 1.0d-2;
		u(1:n,1) = 1.0d0;
		
		% Initial point
		
		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end

% ----------------------------------------------------------------------

	% 	KW2
	%   I.Y. Kim, O.L. de Weck, Adaptive weighted-sum method for bi-objective optimization: Pareto front generation
	
	if ( strcmp(problem,'KW2' )  ) 
	
		% Number of variables

		n = 2;
		
		% Number of objectives
		
		m = 2;

		% Box constraints

		l(1:n,1) = -3.0d0;
		u(1:n,1) =  3.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return
    end
	
% ----------------------------------------------------------------------	

	%   LE1
	%   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit

	
	if ( strcmp(problem,'LE1' )  ) 
	
		% Number of variables

		n = 2;

		% Number of objectives

		m = 2;

		% Box constraints

		l(1:n,1) =  1.0d0;
		u(1:n,1) =  1.0d1;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return
    end
	
% ----------------------------------------------------------------------
	
	% Lov1
	% "Singular Continuation: Generating Piecewise Linear Approximations to Pareto Sets via Global Analysis"
	
	if ( strcmp(problem,'Lov1' )  ) 
	
		% Number of variables
		
		n = 2;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = -1.0d1;
		u(1:n,1) =  1.0d1;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
	
		return
    end
	 
	
% ----------------------------------------------------------------------
	
	% Lov2
	% "Singular Continuation: Generating Piecewise Linear Approximations to Pareto Sets via Global Analysis"
	
	if ( strcmp(problem,'Lov2' )  ) 

		% Number of variables
		
		n = 2;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = -0.75d0;
		u(1:n,1) =  0.75d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return
    end
	 	
	
% ----------------------------------------------------------------------
	
	% Lov3
	% "Singular Continuation: Generating Piecewise Linear Approximations to Pareto Sets via Global Analysis"
	
	if ( strcmp(problem,'Lov3' )  ) 
	
		% Number of variables
		
		n = 2;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = -2.0d1;
		u(1:n,1) =  2.0d1;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end
	
% ----------------------------------------------------------------------
	
	% Lov4
	% "Singular Continuation: Generating Piecewise Linear Approximations to Pareto Sets via Global Analysis"
	
	if ( strcmp(problem,'Lov4' )  ) 
	
		% Number of variables
		
		n = 2;
		
		% Number of objectives
		
		m = 2;

		% Box constraints

		l(1:n,1) = -2.0d1;
		u(1:n,1) =  2.0d1

		% Initial point;

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return
    end
	
% ----------------------------------------------------------------------
	
	% Lov5
	% "Singular Continuation: Generating Piecewise Linear Approximations to Pareto Sets via Global Analysis"
	
	if ( strcmp(problem,'Lov5' )  ) 
	
		% Number of variables
		
		n = 3;
		
		% Number of objectives
		
		m = 2;

		% Box constraints

		l(1:n,1) = -2.0d0;
		u(1:n,1) =  2.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end
	
% ----------------------------------------------------------------------
	
	% Lov6
	% "Singular Continuation: Generating Piecewise Linear Approximations to Pareto Sets via Global Analysis"
	
	if ( strcmp(problem,'Lov6' )  ) 
	
		% Number of variables
		
		n = 6;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1) = 0.1d0;
		u(1) = 0.425;
		
		l(2:6) = - 0.16d0;
		u(2:6) =   0.16d0;
        
        l = l';
        u = u';

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return
    end

% ----------------------------------------------------------------------

	%  LTDZ
	%	 Combining convergence and diversity in evolutionary multiobjective optimization
	
	if ( strcmp(problem,'LTDZ' )  ) 
	
		% Number of variables

		n = 3;
		
		% Number of objectives
		
		m = 3;

		% Box constraints

		l(1:n,1) = 0.0d0;
		u(1:n,1) = 1.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end
	
% ----------------------------------------------------------------------	

	% 	MGH9 (Gaussian)

	% More, J.J., Garbow, B.S., Hillstrom, K.E.: Testing unconstrained optimization software. ACM T. Math.
	% Softw. 7(1), 17???41 (1981)
	% See also: Mita,K. ,Fukuda,E.H. ,Yamashita,N.: Nonmonotone linesearches for unconstrained multiobjective 
	% optimization problems. J. Glob. Optim. 75(1), 63???90 (2019)

	if ( strcmp(problem,'MGH9' )  ) 
	
		% Number of variables

		n = 3;
		
		% Number of objectives
		
		m = 15;

		% Box constraints

		l(1:n,1) = - 2.0d0;
		u(1:n,1) =   2.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
	
		return
	 		
    end
	
% ----------------------------------------------------------------------

	% 	MGH16 (Brown and Dennis)

	% More, J.J., Garbow, B.S., Hillstrom, K.E.: Testing unconstrained optimization software. ACM T. Math.
	% Softw. 7(1), 17???41 (1981)
	% See also: Mita,K. ,Fukuda,E.H. ,Yamashita,N.: Nonmonotone linesearches for unconstrained multiobjective 
	% optimization problems. J. Glob. Optim. 75(1), 63???90 (2019)
	
	if ( strcmp(problem,'MGH16' )  ) 
	
		% Number of variables

		n = 4;
		
		% Number of objectives
		
		m = 5;
		
		% Box constraints

		l(1) = - 2.5d1;
		u(1) =   2.5d1;
		
		l(2:3) = - 5.0d0;
		u(2:3) =   5.0d0;

		l(4) = - 1.0d0;
		u(4) =   1.0d0;
        
        l = l';
        u = u';

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end
	
% ----------------------------------------------------------------------

	% 	MGH26 (Trigonometric)

	% More, J.J., Garbow, B.S., Hillstrom, K.E.: Testing unconstrained optimization software. ACM T. Math.
	% Softw. 7(1), 17???41 (1981)
	% See also: Mita,K. ,Fukuda,E.H. ,Yamashita,N.: Nonmonotone linesearches for unconstrained multiobjective 
	% optimization problems. J. Glob. Optim. 75(1), 63???90 (2019)

	
	if ( strcmp(problem,'MGH26' )  ) 
	
		% Number of variables

		n = 4;
														
		% Number of objectives
		
		m = n;
		
		% Box constraints

		l(1:n,1) = - 1.0d0;
		u(1:n,1) =   1.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return
    end
	
% ----------------------------------------------------------------------

	% 	MGH33 (Linear function - rank 1)

	% More, J.J., Garbow, B.S., Hillstrom, K.E.: Testing unconstrained optimization software. ACM T. Math.
	% Softw. 7(1), 17???41 (1981)
	% See also: Mita,K. ,Fukuda,E.H. ,Yamashita,N.: Nonmonotone linesearches for unconstrained multiobjective 
	% optimization problems. J. Glob. Optim. 75(1), 63???90 (2019)
	
	if ( strcmp(problem,'MGH33' )  ) 
	
		% Number of variables

		n = 10;
		
		% Number of objectives (m>=n)
		
		m = n;
		
		% Box constraints

		l(1:n,1) = - 1.0d0;
		u(1:n,1) =   1.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end
	
% ----------------------------------------------------------------------	

	%   MHHM2
	%   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit

	
	if ( strcmp(problem,'MHHM2' )  ) 
	
		% Number of variables

		n = 2;

		% Number of objectives

		m = 3;
		
		% Box constraints

		l(1:n,1) = 0.0d0;
		u(1:n,1) = 1.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return
    end
	
	
% ----------------------------------------------------------------------

	%   MLF1
	%   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'MLF1' )  ) 
	
		% Number of variables

		n = 1;
		
		% Number of objectives
		
		m = 2;

		% Box constraints

		l(1:n,1) = 0.0d0;
		u(1:n,1) = 2.0d1;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end
	 	
	
% ----------------------------------------------------------------------

	%   MLF2
	%   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'MLF2' )  ) 

		% Number of variables
		
		n = 2;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = - 1.0d2;
		u(1:n,1) =   1.0d2;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end
	 	
	
% ----------------------------------------------------------------------
	
	%  MMR1
	%  Box-constrained multi-objective optimization: A gradient-like method without ??????a priori?????? scalarization	
	
	if ( strcmp(problem,'MMR1' )  ) 
	
		% Number of variables

		n = 2;
		
		% Number of objectives
		
		m = 2;

		% Box constraints

		l(1) = 0.1d0;
		u(1) = 1.0d0;
		
		l(2) = 0.0d0;
		u(2) = 1.0d0;
        
        l = l';
        u = u';

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return

    end
	 
	
% ----------------------------------------------------------------------
	
	%  MMR2
	%  Box-constrained multi-objective optimization: A gradient-like method without ??????a priori?????? scalarization	
	
	if ( strcmp(problem,'MMR2' )  ) 
	
		% Number of variables

		n = 2;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = 0.0d0;
		u(1:n,1) = 1.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return

    end
	 
	
% ----------------------------------------------------------------------
	
	%  MMR3
	%  Box-constrained multi-objective optimization: A gradient-like method without ??????a priori?????? scalarization	
	
	if ( strcmp(problem,'MMR3' )  ) 
	
		% Number of variables

		n = 2;
		
		% Number of objectives
		
		m = 2;

		% Box constraints

		l(1:n,1) = -1.0d0;
		u(1:n,1) =  1.0d0;
		
		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return
    end
	
% ----------------------------------------------------------------------
	
	%  MMR4
	%  Box-constrained multi-objective optimization: A gradient-like method without ??????a priori?????? scalarization	
	
	if ( strcmp(problem,'MMR4' )  ) 

		% Number of variables

		n = 3;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = 0.0d0;
		u(1:n,1) = 4.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return

    end
	 		
	
% ----------------------------------------------------------------------

	%  MOP2
	%  A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit		
	
	if ( strcmp(problem,'MOP2' )  ) 
	
		% Number of variables

		n = 2;
		
		% Number of objectives
		
		m = 2;

		% Box constraints

		l(1:n,1) = - 1.0d0;
		u(1:n,1) =   1.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return

    end
	 		
	
% ----------------------------------------------------------------------

	%  MOP3	
	%  A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit	
	
	if ( strcmp(problem,'MOP3' )  ) 
	
		% Number of variables
	
		n = 2;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = - pi;
		u(1:n,1) =   pi;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return
    end
	
% ----------------------------------------------------------------------

	%  MOP5
	%  A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'MOP5' )  ) 

		% Number of variables

		n = 2;
		
		% Number of objectives
		
		m = 3;

		% Box constraints

		l(1:n,1) = - 1.0d0;
		u(1:n,1) =   1.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return
    end
	
% ----------------------------------------------------------------------

	%  MOP6
	%  A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'MOP6' )  ) 
	
		% Number of variables

		n = 2;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = 0.0d0;
		u(1:n,1) = 1.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return
    end
	
% ----------------------------------------------------------------------

	%  MOP7
	%  A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'MOP7' )  ) 
	
		% Number of variables
		
		n = 2;
		
		% Number of objectives
		
		m = 3;
		
		% Box constraints

		l(1:n,1) = - 4.0d2;
		u(1:n,1) =   4.0d2;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return
    end
	
% ----------------------------------------------------------------------

	%  PNR
	%	 M. Preuss, B. Naujoks, and G. Ru forlph, Pareto set and EMOA behaviour for simple
	%  multimodal multiobjective functions, In: Parallel Problem O Solving from Nature-PPSN IX.
	%  Springer. Berlin, 2006, pp. 513???522.
	
	if ( strcmp(problem,'PNR' )  ) 

		% Number of variables

		n = 2;
		
		% Number of objectives
		
		m = 2;

		% Box constraints

		l(1:n,1) = - 2.0d0;
		u(1:n,1) =   2.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end		
	
% ----------------------------------------------------------------------	
	
	%   QV1
	%   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'QV1' )  ) 
	
		% Number of variables
		
		n = 10;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = 1.0d-2;
		u(1:n,1) = 5.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end
	
% ----------------------------------------------------------------------	

	% SD
	% Stadler, W., Dauer, J.: Multicriteria optimization in engineering: a tutorial and survey. In: Kamat, M.P.
    % (ed.) Progress in Aeronautics and Astronautics: Structural Optimization: Status and Promise, vol. 150,
    % pp. 209???249. American Institute of Aeronautics and Astronautics, Reston (1992)
 
	
	if ( strcmp(problem,'SD' )  ) 
	
		% Number of variables
		
		n = 4;
		
		% Number of objectives
		
		m = 2;

		% Box constraints

		l(1) = 1.0d0;
		u(1) = 3.0d0;
		
		l(2:3) = sqrt(2.0d0);
		u(2:3) = 3.0d0;
		
		l(4) = 1.0d0;
		u(4) = 3.0d0;
        
        l = l';
        u = u';

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		 return
    end

% ----------------------------------------------------------------------	

	% SLCDT1
	% Convergence of stochastic search algorithms to finite size pareto set approximations
 
	
	if ( strcmp(problem,'SLCDT1' )  ) 
	
		% Number of variables
		
		n = 2;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = - 1.5d0;
		u(1:n,1) =   1.5d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
         end 

         return
    end

% ----------------------------------------------------------------------

	%  SLCDT2
	%	 Convergence of stochastic search algorithms to finite size pareto set approximations

	
	if ( strcmp(problem,'SLCDT2' )  ) 
	
		% Number of variables

		n = 10;
		
		% Number of objectives
		
		m = 3;

		% Box constraints

		l(1:n,1) = - 1.0d0;
		u(1:n,1) =   1.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
         return
    end
	
% ----------------------------------------------------------------------

	%   SP1
	%  	A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'SP1' )  ) 
	
		% Number of variables

		n = 2;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = - 1.0d2;
		u(1:n,1) =   1.0d2;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return

    end
	 		
	
% ----------------------------------------------------------------------

	%   SSFYY2
	%   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'SSFYY2' )  ) 
	
		% Number of variables
	
		n = 1;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = - 1.0d2;
		u(1:n,1) =   1.0d2;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return

    end
	 		
	
% ----------------------------------------------------------------------

	%   SK1
	%   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'SK1' )  ) 

		% Number of variables
	
		n = 1;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = - 1.0d2;
		u(1:n,1) =   1.0d2;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end
	
% ----------------------------------------------------------------------

	%   SK2
	%  	A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'SK2' )  ) 
	
		% Number of variables

		n = 4;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = - 1.0d1;
		u(1:n,1) =   1.0d1;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return
    end
	
% ----------------------------------------------------------------------

	%   TKLY1
	%  	A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'TKLY1' )  ) 
	
		% Number of variables

		n = 4;
		
		% Number of objectives
		
		m = 2;

		% Box constraints

		l(1) = 0.1d0;
		u(1) = 1.0d0;
		
		l(2:4) = 0.0d0;
		u(2:4) = 1.0d0;
        
        l = l';
        u = u';

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end
	
% ----------------------------------------------------------------------

	% 	Toi4
	
	% Toint, P.L.: Test problems for partially separable optimization and results for the routine PSPMIN.
	% Technical Report, The University of Namur, Department of Mathematics, Belgium (1983)
	% See also: Mita,K. ,Fukuda,E.H. ,Yamashita,N.: Nonmonotone linesearches for unconstrained multiobjective 
	% optimization problems. J. Glob. Optim. 75(1), 63???90 (2019)
	
	if ( strcmp(problem,'Toi4' )  ) 
	
		% Number of variables

		n = 4;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = - 2.0d0;
		u(1:n,1) =   5.0d0;
		
		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end
	
% ----------------------------------------------------------------------

	% 	Toi8 (TRIDIA)
	
	% Toint, P.L.: Test problems for partially separable optimization and results for the routine PSPMIN.
	% Technical Report, The University of Namur, Department of Mathematics, Belgium (1983)
	% See also: Mita,K. ,Fukuda,E.H. ,Yamashita,N.: Nonmonotone linesearches for unconstrained multiobjective 
	% optimization problems. J. Glob. Optim. 75(1), 63???90 (2019)
	
	if ( strcmp(problem,'Toi8' )  ) 
	
		% Number of variables

		n = 3;
		
		% Number of objectives
		
		m = n;

		% Box constraints

		l(1:n,1) = - 1.0d0;
		u(1:n,1) =   1.0d0;
		
		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end
		
% ----------------------------------------------------------------------

	% 	Toi9 (Shifted TRIDIA)
	
	% Toint, P.L.: Test problems for partially separable optimization and results for the routine PSPMIN.
	% Technical Report, The University of Namur, Department of Mathematics, Belgium (1983)
	% See also: Mita,K. ,Fukuda,E.H. ,Yamashita,N.: Nonmonotone linesearches for unconstrained multiobjective 
	% optimization problems. J. Glob. Optim. 75(1), 63???90 (2019)
	
	if ( strcmp(problem,'Toi9' )  ) 
	
		% Number of variables

		n = 4;
												
		% Number of objectives
		
		m = n;
		
		% Box constraints

		l(1:n,1) = - 1.0d0;
		u(1:n,1) =   1.0d0;
		
		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
	
		return
    end
	
% ----------------------------------------------------------------------

	% 	Toi10 (Rosenbrock)
	
	% Toint, P.L.: Test problems for partially separable optimization and results for the routine PSPMIN.
	% Technical Report, The University of Namur, Department of Mathematics, Belgium (1983)
	% See also: Mita,K. ,Fukuda,E.H. ,Yamashita,N.: Nonmonotone linesearches for unconstrained multiobjective 
	% optimization problems. J. Glob. Optim. 75(1), 63???90 (2019)
	
	if ( strcmp(problem,'Toi10' )  ) 
	
		% Number of variables

		n = 4;
										
		% Number of objectives
		
		m = n - 1;
		
		% Box constraints

		l(1:n,1) = - 2.0d0;
		u(1:n,1) =   2.0d0;
		
		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end			
	
% ----------------------------------------------------------------------

	%   VU1
	%   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'VU1' )  ) 
	
		% Number of variables

		n = 2;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = - 3.0d0;
		u(1:n,1) =   3.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return
    end
	
% ----------------------------------------------------------------------

	%   VU2
	%   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'VU2' )  ) 
	
		% Number of variables

		n = 2;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = - 3.0d0;
		u(1:n,1) =   3.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end
	 
	
% ----------------------------------------------------------------------

	%   ZDT1
	%   Comparison of Multiobjective Evolutionary Algorithms: Empirical Results
	
	if ( strcmp(problem,'ZDT1' )  ) 
	
		% Number of variables

		n = 30;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = 1.0d-2;
		u(1:n,1) = 1.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 	
		
		return
    end
	
% ----------------------------------------------------------------------

	%   ZDT2
	%   Comparison of Multiobjective Evolutionary Algorithms: Empirical Results
	
	if ( strcmp(problem,'ZDT2' )  ) 
	
		% Number of variables

		n = 30;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = 0.0d0;
		u(1:n,1) = 1.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end
	
% ----------------------------------------------------------------------

	%   ZDT3
	%   Comparison of Multiobjective Evolutionary Algorithms: Empirical Results
	
	if ( strcmp(problem,'ZDT3' )  ) 
	
		% Number of variables

		n = 30;
		
		% Number of objectives
		
		m = 2;
		
		% Box constraints

		l(1:n,1) = 1.0d-2;
		u(1:n,1) = 1.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end
	
% ----------------------------------------------------------------------

	%   ZDT4
	%   Comparison of Multiobjective Evolutionary Algorithms: Empirical Results
	
	if ( strcmp(problem,'ZDT4' )  ) 
	
		% Number of variables

		n = 30;
		
		% Number of objectives
		
		m = 2;

		% Box constraints

		l(1) = 1.0d-2;
		u(1) = 1.0d0;
		
		l(2:n) = - 5.0d0;
		u(2:n) =   5.0d0;
        
        l = l';
        u = u';

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return
    end
	
% ----------------------------------------------------------------------

	%   ZDT6
	%   Comparison of Multiobjective Evolutionary Algorithms: Empirical Results
	
	if ( strcmp(problem,'ZDT6' )  ) 
	
		% Number of variables

		n = 10;
		
		% Number of objectives
		
		m = 2;

		% Box constraints

		l(1:n,1) = 0.01d0;
		u(1:n,1) = 1.0d0;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 

		return
    end
	
% ----------------------------------------------------------------------

	%   ZLT1
	%   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( strcmp(problem,'ZLT1' )  ) 
	
		% Number of variables

		n = 10;
		
		% Number of objectives
		
		m = 5;
		
		% Box constraints

		l(1:n,1) = - 1.0d3;
		u(1:n,1) =   1.0d3;

		% Initial point

		 for i = 1:n
			x(i,1) = l(i) + ( u(i) - l(i) ) *  rand;
		 end 
		
		return

    end
	 			



	