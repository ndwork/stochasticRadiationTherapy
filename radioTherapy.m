
function [plan,objectiveValues] = radioTherapy( A, d, nVoxsPerStructure, w, w0, varargin )

  p = inputParser;
  p.addOptional( 'alg', 'stochProxGrad', @(x) true );
  p.addParameter( 'batchSize', 32, @isnumeric );
  p.addParameter( 't', [], @isnumeric );
  p.addParameter( 'N', 100, @isnumeric );
  p.parse( varargin{:} );
  alg = p.Results.alg;
  batchSize = p.Results.batchSize;
  t = p.Results.t;
  N = p.Results.N;

  switch alg
    case 'fista'
      if nargout > 1
        [plan,objectiveValues] = radioTherapy_fista( A, d, nVoxsPerStructure, w, w0, ...
          't', t, 'N', N );
      else
        plan = radioTherapy_fista( A, d, nVoxsPerStructure, w, w0, 't', t, 'N', N );
      end

    case 'fista_wLS'
      if nargout > 1
        [plan,objectiveValues] = radioTherapy_fista_wLS( A, d, nVoxsPerStructure, w, w0, ...
          'N', N );
      else
        plan = radioTherapy_fista_wLS( A, d, nVoxsPerStructure, w, w0, 'N', N );
      end

    case 'proxGrad'
      if nargout > 1
        [plan,objectiveValues] = radioTherapy_proxGrad( A, d, nVoxsPerStructure, w, w0, ...
          't', t, 'N', N );
      else
        plan = radioTherapy_proxGrad( A, d, nVoxsPerStructure, w, w0, 't', t, 'N', N );
      end

    case 'stochProxGrad'
      if nargout > 1
        [plan,objectiveValues] = radioTherapy_stochProxGrad( A, d, nVoxsPerStructure, w, w0, ...
          't', t, 'N', N, 'batchSize', batchSize );
      else
        plan = radioTherapy_stochProxGrad( A, d, nVoxsPerStructure, w, w0, ...
          't', t, 'N', N, 'batchSize', batchSize );
      end
      
    case 'proxSVRG'
      if nargout > 1
        [plan,objectiveValues] = radioTherapy_proxSVRG( A, d, nVoxsPerStructure, w, w0, ...
          't', t, 'N', N, 'batchSize', batchSize );
      else
        plan = radioTherapy_proxSVRG( A, d, nVoxsPerStructure, w, w0, ...
          't', t, 'N', N, 'batchSize', batchSize );
      end
      
    otherwise
      error( 'Unrecognized algorithm' );
  end

end
