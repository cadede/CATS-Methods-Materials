function    v=norm2(X)

%     v=norm2(X)
%     Returns the row-wise 2-norm of matrix X, i.e., the square-root
%		of the sum of the squares for each row. If X is a vector, norm2()
%	   is equivalent to the built-in function norm(). But if X is a matrix
%		e.g., a triaxial accelerometer or magnetometer matrix, norm() gives
%		the overall norm of the matrix whereas norm2() gives the vector norm
%		of each row (i.e., the field strength in the case of a magnetometer
%		matrix.
%
%		Input:
%		X is a vector or matrix.
%
%		Returns:
%		v is the row-wise vector norm of X if X is a matrix. If X is a vector
%		 (row or column), v is the vector norm.
%
%		Example:
%		 v = norm2([0.2 0.4 -0.7;-0.3 1.1 0.1])
%		 returns: v=[0.83066;1.14455]
%
%     Valid: Matlab, Octave
%     markjohnson@st-andrews.ac.uk
%     Last modified: 10 May 2017

[m,n] = size(X) ;
if m==1 | n==1,
   v = norm(X) ;
else
   v = sqrt(sum(abs(X).^2,2)) ;
end

