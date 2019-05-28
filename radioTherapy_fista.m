
function [xStar,objectiveValues] = radioTherapy_fista( A, d, nVoxsPerStructure, w, w0, varargin )

  p = inputParser;
  p.addParameter( 't', 0.001, @isnumeric );
  p.addParameter( 'N', 100, @isnumeric );  % Number of proxGrad iterations
  p.addParameter( 'verbose', true, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  t = p.Results.t;
  N = p.Results.N;
  verbose = p.Results.verbose;

  AT = A';
  A0 = A( 1 : nVoxsPerStructure(1), : );
  d0 = d( 1 : nVoxsPerStructure(1) );

  function out = g( x )
    Ax = AT' * x;
    tmp0 = max( d0 - Ax(1:nVoxsPerStructure(1)), 0 );
    out = 0.5 * sum( w0 .* tmp0 .* tmp0 );
    tmp = max( Ax - d, 0 );
    out = out + 0.5 * sum( w .* tmp .* tmp );
  end

  function out = gPrime( x )
    Ax = AT' * x;
    out = -A0' * ( w0 .* max( d0 - Ax(1:nVoxsPerStructure(1)), 0 ) );
    out = out + A' * ( w .* max( Ax - d, 0 ) );
  end

  function out = h( x )
    out=0;
    if min( x ) < 0, out = Inf; end
  end

  proxth = @(x,t) max( x, 0 );


  x0 = zeros( size( A, 2 ), 1 );
  if nargout > 1
    [xStar,objectiveValues] = fista( x0, @g, @gPrime, proxth, 't', t, 'h', @h, ...
      'N', N, 'verbose', verbose );
  else
    xStar = fista( x0, @g, @gPrime, proxth, 't', t, 'h', @h, ...
      'N', N, 'verbose', verbose );
  end
end
