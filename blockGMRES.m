function [x, error, iter, flag] = myblockgmres1( A, x, b, max_it, tol )

% input   A        complex nonsymmetric positive definite matrix
%         x        complex initial guess vector block
%         b        complex right hand side vector block
%         max_it   INTEGER maximum number of iterations
%         tol      complex error tolerance
%
% output  x        complex solution vector block
%         error    real error norm
%         iter     INTEGER number of iterations performed
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it

   iter = 0;                                         
   flag = 0;
   numRHS = min(size(b));
   
   bnrm2 = norm( b );
   if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end

   r =  b-A*x ;
   error = norm( r ) / bnrm2;
   if ( error < tol ) return, end

   [n,n] = size(A);                                  
   m = max_it;
   p = numRHS;
   V(1:n, 1:numRHS) = zeros(n,numRHS);
   H = [];
   QH = [];
   cs = [];
   sn = [];
   e1 = zeros(n, numRHS);
   for i=1:numRHS
    e1(i,i) = 1 + 0i;    
   end

   r = b-A*x;
   %calculate QR of initial residue r0
   [r0,R] = qr(r, 0);
   V(:,1:numRHS) = r0(:,1:numRHS);
   g = e1*R(1:numRHS, 1:numRHS);
    
     for j = 1:max_it,
         Avj = A*V(:,(j-1)*p+1:j*p);
         for i = 1:j
             H((i-1)*p+1:i*p, (j-1)*p+1:j*p) = V(:,(i-1)*p+1:i*p)'*Avj;
             Avj = Avj - V(:,(i-1)*p+1:i*p)*H((i-1)*p+1:i*p, (j-1)*p+1:j*p);
         end
         [Q2, R2] = qr(Avj, 0); 
         H(j*p+1:(j+1)*p, (j-1)*p+1:j*p) = R2(1:numRHS, 1:numRHS);
         V(:, j*p+1:(j+1)*p) = Q2(:, 1:numRHS);
         %QH is H which gets reduced through Given's
         QH = [QH; zeros(p, (j-1)*p)];
         QH = [QH zeros((j+1)*p, p)];
         QH(:, (j-1)*p+1:j*p) = H(:, (j-1)*p+1:j*p);
         
         %apply Givens rotation in all previous iterations
         %to new block of columns
         for k = (j-1)*p+1:j*p
             %columns-wise
             for l = 1:k-1
                %reverse order of rows in that block
                for i = p:-1:1
                        % apply Givens rotation
                        col = k;
                        row1 = l + i - 1;
                        row2 = l + i;
                        rr = QH(row1, col);
                        hh = QH(row2, col);
                        QH(row1, col)  = cs(l, i)*rr - sn(l, i)*hh;
                        QH(row2, col) = conj(sn(l, i))*rr + cs(l, i)*hh;
                end
             end
             %new column reduction
             for i = p:-1:1
                col = k;
                l = k;
                row1 = l + i - 1;
                row2 = l + i;
                rr = QH(row1, col);
                hh = QH(row2, col);
                rnorm = norm(rr);
                hnorm = norm(hh);
                rhnorm = sqrt(rnorm*rnorm + hnorm*hnorm);
                sn(l, i) =  -1*(hh*rr) / (rnorm*rhnorm);
                cs(l, i) =  real(rnorm/rhnorm);
                temp = conj(sn(l, i))*g(row1, :) + cs(l, i)*g(row2, :);
                g(row1, :) = g(row1, :)*cs(l, i) + sn(l, i)*g(row2, :);
                g(row2, :) = temp;
                QH(row1, col) = cs(l, i)*rr - sn(l, i)*hh;
                QH(row2, col) = 0.0;
             end
         end
         'iter=', j
         noexit = 0;
         %F-norm
         for i=1:numRHS
             error(i)  = norm(g(j*p+1:(j+1)*p, i));
         end
         norm(error)
         if(norm(error) > tol)
           noexit = 1;
         end
         %solve triangular systtem
         if ( noexit == 0 || j == max_it),  
            y = zeros(size(b));
            y = QH(1:j*p,1:j*p) \ g(1:j*p,:);
            %y = y*norm(b);
            y1(1:j*p, :) = y(1:j*p, :);
            x = x + V(:,1:j*p)*y1;
            break;
         end
     end

   if ( error > tol ) flag = 1; end;
