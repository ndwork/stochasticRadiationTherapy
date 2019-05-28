
function [xStar,objectiveValues] = radioTherapy_stochProxGrad( ...
  A, d, nVoxsPerStructure, w, w0, varargin )

  p = inputParser;
  p.addParameter( 't', 0.001, @isnumeric );  % step size
  p.addParameter( 'batchSize', 16, @isnumeric );
  p.addParameter( 'N', 100, @isnumeric );  % Number of iterations
  p.addParameter( 'verbose', true, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  t = p.Results.t;
  batchSize = p.Results.batchSize;
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

  nRowsA = size( A, 1 );
  flagsA0 = zeros( nRowsA, 1 );
  flagsA0( 1 : nVoxsPerStructure(1) ) = 1;

  nBatches = ceil( nRowsA / batchSize );

  for k=0:N-1

    if calculateObjectiveValues > 0
      % Calculate objective on y
      Ay = AT' * y;
      if min( y ) < 0, hOfy=Inf; else, hOfy=0; end
      objectiveValues(k+1) = g(Ay) + hOfy;
    end

    permIndxs = randperm( nRowsA )';
    permA0Flags = flagsA0( permIndxs );

    for batchIndx = 1 : nBatches

      currentIndx = (batchIndx-1)*batchSize + 1;
      endIndx = min( currentIndx+batchSize-1, nRowsA );
      theseIndxs = permIndxs( currentIndx : endIndx );
      theseA0Flags = permA0Flags( currentIndx : endIndx );

      subAT = AT( :, theseIndxs );
      subA = subAT';
      subd = d( theseIndxs );
      subw = w( theseIndxs );

      subAy = subAT' * y;
      gPrime_y = subA' * ( subw .* max( subAy - subd, 0 ) );

      if max( theseA0Flags ) > 0
        theseA0Indxs = theseIndxs( theseA0Flags == 1 );
        subA0 = A0( theseA0Indxs, : );
        ptvAy = subAy( theseA0Flags == 1 );
        ptvd0 = d0( theseA0Indxs );
        ptvW0 = w0( theseA0Indxs );
        gPrime_y = gPrime_y - subA0' * ( ptvW0 .* max( ptvd0 - ptvAy, 0 ) );
      end

      x = y - t * gPrime_y;
      y = max( x, 0 );  % proximal operator of h

    end

    if verbose
      formatString = ['%', num2str(ceil(log10(N))), '.', num2str(ceil(log10(N))), 'i' ];
      verboseString = [ 'stochProxGrad Iteration: ', num2str(k,formatString) ];
      if calculateObjectiveValues > 0
        verboseString = [ verboseString, ',  objective: ', num2str( objectiveValues(k+1) ) ];   %#ok<AGROW>
      end
      disp( verboseString );
    end
  end

  xStar = y;
end
