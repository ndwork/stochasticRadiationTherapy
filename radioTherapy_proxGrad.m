
function [xStar,objectiveValues] = radioTherapy_proxGrad( A, d, nVoxsPerStructure, w, w0, varargin )

  p = inputParser;
  p.addParameter( 't', 0.001, @isnumeric );
  p.addParameter( 'N', 100, @isnumeric );  % Number of proxGrad iterations
  p.addParameter( 'verbose', true, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  t = p.Results.t;
  N = p.Results.N;
  verbose = p.Results.verbose;


  y = zeros( size( A, 2 ), 1 );

  calculateObjectiveValues = 0;
  if nargout > 1
    objectiveValues = zeros(N,1);
    calculateObjectiveValues = 1;
  end

  AT = A';
  A0 = A( 1 : nVoxsPerStructure(1), : );
  d0 = d( 1 : nVoxsPerStructure(1) );

  function out = g( Ax )
    tmp0 = max( d0 - Ax(1:nVoxsPerStructure(1)), 0 );
    out = 0.5 * sum( w0 .* tmp0 .* tmp0 );
    tmp = max( Ax - d, 0 );
    out = out + 0.5 * sum( w .* tmp .* tmp );
  end

  for k=0:N-1

    Ay = AT' * y;

    if calculateObjectiveValues > 0
      % Calculate objective on y
      if min( y ) < 0, hOfy=Inf; else, hOfy=0; end
      objectiveValues(k+1) = g(Ay) + hOfy;
    end

    gPrime_y = AT * ( w .* max( Ay - d, 0 ) ) ...
               - A0' * ( w0 .* max( d0 - Ay(1:nVoxsPerStructure(1)), 0 ) );

    x = y - t * gPrime_y;
    y = max( x, 0 );  % proximal operator of h

    if verbose
      formatString = ['%', num2str(ceil(log10(N))), '.', num2str(ceil(log10(N))), 'i' ];
      verboseString = [ 'proxGrad Iteration: ', num2str(k,formatString) ];
      if calculateObjectiveValues > 0
        verboseString = [ verboseString, ',  objective: ', num2str( objectiveValues(k+1) ) ];   %#ok<AGROW>
      end
      disp( verboseString );
    end

  end

  xStar = y;
end
