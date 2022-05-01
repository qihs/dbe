function BVS = BooleanVectorSearch(A,PRECISION)
  
  if nargin <2
    PRECISION = 1.0e-14; % set the default precision 
  end
  
  BVS = [];
  [m,n] = size(A);  %% compute the dimension of the input matrix
  y = A(1:m,1);  %% the first column
  B = A(1:m,2:n)-y;  %% generators for the subspace
  BT = B';
  CT = rref(BT); %% rref is supposed to compute the canonical form
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Compute the rank profile %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  RankProfile = [];  % Initialization of the rank profile
  CRankP = []; % Initialization of the complementary of the rank profile
  j = 0; 
  temp = eye(n-1);
  for i=1:n-1
    k=j+1;
    for j=k:m
      if CT(i,j) ~= 0 
	if compareV(CT(:,j) - temp(:,i),PRECISION) %% 'if' condition verify the assumption for 'rref'
	  RankProfile = [RankProfile, j];
	  break;
	else
	  error('error in rankprofile');
	end
      else
        CRankP = [CRankP, j];
      end
    end
    
    if i==n-1 | j == m
      if CT(i,j)~=0
        i = i + 1;
      end % if this is satisfied, the following rows are zero rows
      break
    end
   
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  nc = i-1; % record the columns of the basis
  if nc == m
    BVS = eye(m);
  else
    C = CT(1:nc,:)'; % columns of C constitute a basis of the subspace
    D = C(CRankP,1:nc); % the submatrix of C by removing the identity matrix
    z = y(CRankP) - D*y(RankProfile); 
    temp = eye(m);
    
    %%% Check the rank profile
    for i=1:nc
      if compareV(z+D(:,i),PRECISION)
        BVS = [BVS;temp(:,RankProfile(i))'];
      end
    end

    %%% Check the complementary of rank profile
    for i = CRankP
%       z- temp(CRankP,i)
%       compareV(z-temp(CRankP,i),PRECISION)
      if compareV(z-temp(CRankP,i),PRECISION)
        BVS = [BVS;temp(:,i)'];
      end
    end
    
    BVS = BVS';
  end
end

function Bool = compareV(v,PRECISION)
  m = size(v,1);
  Bool = 1;
  for i = 1:m
    if abs(v(i)) > PRECISION
      Bool = 0;
      break;
    end
  end
end


