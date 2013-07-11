package Jama;
import Jama.util.*;

/** Eigenvalues and eigenvectors of a real matrix. 
<P>
    If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is
    diagonal and the eigenvector matrix V is orthogonal.
    I.e. A = V.times(D.times(V.transpose())) and 
    V.times(V.transpose()) equals the identity matrix.
<P>
    If A is not symmetric, then the eigenvalue matrix D is block diagonal
    with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
    lambda + i*mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda].  The
    columns of V represent the eigenvectors in the sense that A*V = V*D,
    i.e. A.times(V) equals V.times(D).  The matrix V may be badly
    conditioned, or even singular, so the validity of the equation
    A = V*D*inverse(V) depends upon V.cond().
**/

public class EigenvalueDecomposition implements java.io.Serializable {

/* ------------------------
   Class variables
 * ------------------------ */

   /** Row and column dimension (square matrix).
   @serial matrix dimension.
   */
   private int n;

   /** Symmetry flag.
   @serial internal symmetry flag.
   */
   private boolean issymmetric;

   /** Arrays for internal storage of eigenvalues.
   @serial internal storage of eigenvalues.
   */
   private float[] d, e;

   /** Array for internal storage of eigenvectors.
   @serial internal storage of eigenvectors.
   */
   private float[][] V;

   /** Array for internal storage of nonsymmetric Hessenberg form.
   @serial internal storage of nonsymmetric Hessenberg form.
   */
   private float[][] H;

   /** Working storage for nonsymmetric algorithm.
   @serial working storage for nonsymmetric algorithm.
   */
   private float[] ort;

/* ------------------------
   Private Methods
 * ------------------------ */

   // Symmetric Householder reduction to tridiagonal form.

   private void tred2 () {
        System.arraycopy(V[n-1], 0, d, 0, n);

      // Householder reduction to tridiagonal form.
   
      for (int i = n-1; i > 0; i--) {
   
         // Scale to avoid under/overflow.
   
         float scale = 0.0f;
         float h = 0.0f;
         for (int k = 0; k < i; k++) {
            scale = scale + Math.abs(d[k]);
         }
         if (scale == 0.0) {
            e[i] = d[i-1];
            for (int j = 0; j < i; j++) {
               d[j] = V[i-1][j];
               V[i][j] = 0.0f;
               V[j][i] = 0.0f;
            }
         } else {
   
            // Generate Householder vector.
   
            for (int k = 0; k < i; k++) {
               d[k] /= scale;
               h += d[k] * d[k];
            }
            float f = d[i-1];
            float g = (float) Math.sqrt(h);
            if (f > 0) {
               g = -g;
            }
            e[i] = scale * g;
            h = h - f * g;
            d[i-1] = f - g;
            for (int j = 0; j < i; j++) {
               e[j] = 0.0f;
            }
   
            // Apply similarity transformation to remaining columns.
   
            for (int j = 0; j < i; j++) {
               f = d[j];
               V[j][i] = f;
               g = e[j] + V[j][j] * f;
               for (int k = j+1; k <= i-1; k++) {
                  g += V[k][j] * d[k];
                  e[k] += V[k][j] * f;
               }
               e[j] = g;
            }
            f = 0.0f;
            for (int j = 0; j < i; j++) {
               e[j] /= h;
               f += e[j] * d[j];
            }
            float hh = f / (h + h);
            for (int j = 0; j < i; j++) {
               e[j] -= hh * d[j];
            }
            for (int j = 0; j < i; j++) {
               f = d[j];
               g = e[j];
               for (int k = j; k <= i-1; k++) {
                  V[k][j] -= (f * e[k] + g * d[k]);
               }
               d[j] = V[i-1][j];
               V[i][j] = 0.0f;
            }
         }
         d[i] = h;
      }
   
      // Accumulate transformations.
   
      for (int i = 0; i < n-1; i++) {
         V[n-1][i] = V[i][i];
         V[i][i] = 1.0f;
         float h = d[i+1];
         if (h != 0.0) {
            for (int k = 0; k <= i; k++) {
               d[k] = V[k][i+1] / h;
            }
            for (int j = 0; j <= i; j++) {
               float g = 0.0f;
               for (int k = 0; k <= i; k++) {
                  g += V[k][i+1] * V[k][j];
               }
               for (int k = 0; k <= i; k++) {
                  V[k][j] -= g * d[k];
               }
            }
         }
         for (int k = 0; k <= i; k++) {
            V[k][i+1] = 0.0f;
         }
      }
      for (int j = 0; j < n; j++) {
         d[j] = V[n-1][j];
         V[n-1][j] = 0.0f;
      }
      V[n-1][n-1] = 1.0f;
      e[0] = 0.0f;
   } 

   // Symmetric tridiagonal QL algorithm.
   
   private void tql2 () {

   //  This is derived from the Algol procedures tql2, by
   //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
   //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
   //  Fortran subroutine in EISPACK.
   
      for (int i = 1; i < n; i++) {
         e[i-1] = e[i];
      }
      e[n-1] = 0.0f;
   
      float f = 0.0f;
      float tst1 = 0.0f;
      float eps = (float) Math.pow(2.0,-52.0);
      for (int l = 0; l < n; l++) {

         // Find small subdiagonal element
   
         tst1 = Math.max(tst1,Math.abs(d[l]) + Math.abs(e[l]));
         int m = l;
         while (m < n) {
            if (Math.abs(e[m]) <= eps*tst1) {
               break;
            }
            m++;
         }
   
         // If m == l, d[l] is an eigenvalue,
         // otherwise, iterate.
   
         if (m > l) {
            int iter = 0;
            do {
               iter = iter + 1;  // (Could check iteration count here.)
   
               // Compute implicit shift
   
               float g = d[l];
               float p = (d[l+1] - g) / (2.0f * e[l]);
               float r = (float) Maths.hypot(p,1.0f);
               if (p < 0) {
                  r = -r;
               }
               d[l] = e[l] / (p + r);
               d[l+1] = e[l] * (p + r);
               float dl1 = d[l+1];
               float h = g - d[l];
               for (int i = l+2; i < n; i++) {
                  d[i] -= h;
               }
               f = f + h;
   
               // Implicit QL transformation.
   
               p = d[m];
               float c = 1.0f;
               float c2 = c;
               float c3 = c;
               float el1 = e[l+1];
               float s = 0.0f;
               float s2 = 0.0f;
               for (int i = m-1; i >= l; i--) {
                  c3 = c2;
                  c2 = c;
                  s2 = s;
                  g = c * e[i];
                  h = c * p;
                  r = (float) Maths.hypot(p,e[i]);
                  e[i+1] = s * r;
                  s = e[i] / r;
                  c = p / r;
                  p = c * d[i] - s * g;
                  d[i+1] = h + s * (c * g + s * d[i]);
   
                  // Accumulate transformation.
   
                  for (int k = 0; k < n; k++) {
                     h = V[k][i+1];
                     V[k][i+1] = s * V[k][i] + c * h;
                     V[k][i] = c * V[k][i] - s * h;
                  }
               }
               p = -s * s2 * c3 * el1 * e[l] / dl1;
               e[l] = s * p;
               d[l] = c * p;
   
               // Check for convergence.
   
            } while (Math.abs(e[l]) > eps*tst1);
         }
         d[l] = d[l] + f;
         e[l] = 0.0f;
      }
     
      // Sort eigenvalues and corresponding vectors.
   
      for (int i = 0; i < n-1; i++) {
         int k = i;
         float p = d[i];
         for (int j = i+1; j < n; j++) {
            if (d[j] < p) {
               k = j;
               p = d[j];
            }
         }
         if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (int j = 0; j < n; j++) {
               p = V[j][i];
               V[j][i] = V[j][k];
               V[j][k] = p;
            }
         }
      }
   }

   // Nonsymmetric reduction to Hessenberg form.

   private void orthes () {
   
      //  This is derived from the Algol procedures orthes and ortran,
      //  by Martin and Wilkinson, Handbook for Auto. Comp.,
      //  Vol.ii-Linear Algebra, and the corresponding
      //  Fortran subroutines in EISPACK.
   
      int low = 0;
      int high = n-1;
   
      for (int m = low+1; m <= high-1; m++) {
   
         // Scale column.
   
         float scale = 0.0f;
         for (int i = m; i <= high; i++) {
            scale = scale + Math.abs(H[i][m-1]);
         }
         if (scale != 0.0) {
   
            // Compute Householder transformation.
   
            float h = 0.0f;
            for (int i = high; i >= m; i--) {
               ort[i] = H[i][m-1]/scale;
               h += ort[i] * ort[i];
            }
            float g = (float) Math.sqrt(h);
            if (ort[m] > 0) {
               g = -g;
            }
            h = h - ort[m] * g;
            ort[m] = ort[m] - g;
   
            // Apply Householder similarity transformation
            // H = (I-u*u'/h)*H*(I-u*u')/h)
   
            for (int j = m; j < n; j++) {
               float f = 0.0f;
               for (int i = high; i >= m; i--) {
                  f += ort[i]*H[i][j];
               }
               f = f/h;
               for (int i = m; i <= high; i++) {
                  H[i][j] -= f*ort[i];
               }
           }
   
           for (int i = 0; i <= high; i++) {
               float f = 0.0f;
               for (int j = high; j >= m; j--) {
                  f += ort[j]*H[i][j];
               }
               f = f/h;
               for (int j = m; j <= high; j++) {
                  H[i][j] -= f*ort[j];
               }
            }
            ort[m] = scale*ort[m];
            H[m][m-1] = scale*g;
         }
      }
   
      // Accumulate transformations (Algol's ortran).

      for (int i = 0; i < n; i++) {
         for (int j = 0; j < n; j++) {
            V[i][j] = (i == j ? 1.0f : 0.0f);
         }
      }

      for (int m = high-1; m >= low+1; m--) {
         if (H[m][m-1] != 0.0) {
            for (int i = m+1; i <= high; i++) {
               ort[i] = H[i][m-1];
            }
            for (int j = m; j <= high; j++) {
               float g = 0.0f;
               for (int i = m; i <= high; i++) {
                  g += ort[i] * V[i][j];
               }
               // float division avoids possible underflow
               g = (g / ort[m]) / H[m][m-1];
               for (int i = m; i <= high; i++) {
                  V[i][j] += g * ort[i];
               }
            }
         }
      }
   }


   // Complex scalar division.

   private transient float cdivr, cdivi;
   private void cdiv(float xr, float xi, float yr, float yi) {
      float r,dd;
      if (Math.abs(yr) > Math.abs(yi)) {
         r = yi/yr;
         dd = yr + r*yi;
         cdivr = (xr + r*xi)/dd;
         cdivi = (xi - r*xr)/dd;
      } else {
         r = yr/yi;
         dd = yi + r*yr;
         cdivr = (r*xr + xi)/dd;
         cdivi = (r*xi - xr)/dd;
      }
   }


   // Nonsymmetric reduction from Hessenberg to real Schur form.

   private void hqr2 () {
   
      //  This is derived from the Algol procedure hqr2,
      //  by Martin and Wilkinson, Handbook for Auto. Comp.,
      //  Vol.ii-Linear Algebra, and the corresponding
      //  Fortran subroutine in EISPACK.
   
      // Initialize
   
      int nn = this.n;
      int nnn = nn-1;
      int low = 0;
      int high = nn-1;
      float eps = (float) Math.pow(2.0,-52.0);
      float exshift = 0.0f;
      float p=0,q=0,r=0,s=0,z=0,t,w,x,y;
   
      // Store roots isolated by balanc and compute matrix norm
   
      float norm = 0.0f;
      for (int i = 0; i < nn; i++) {
         if (i < low | i > high) {
            d[i] = H[i][i];
            e[i] = 0.0f;
         }
         for (int j = Math.max(i-1,0); j < nn; j++) {
            norm = norm + Math.abs(H[i][j]);
         }
      }
   
      // Outer loop over eigenvalue index
   
      int iter = 0;
      while (nnn >= low) {
   
         // Look for single small sub-diagonal element
   
         int l = nnn;
         while (l > low) {
            s = Math.abs(H[l-1][l-1]) + Math.abs(H[l][l]);
            if (s == 0.0) {
               s = norm;
            }
            if (Math.abs(H[l][l-1]) < eps * s) {
               break;
            }
            l--;
         }
       
         // Check for convergence
         // One root found
   
         if (l == nnn) {
            H[nnn][nnn] = H[nnn][nnn] + exshift;
            d[nnn] = H[nnn][nnn];
            e[nnn] = 0.0f;
            nnn--;
            iter = 0;
   
         // Two roots found
   
         } else if (l == nnn-1) {
            w = H[nnn][nnn-1] * H[nnn-1][nnn];
            p = (H[nnn-1][nnn-1] - H[nnn][nnn]) / 2.0f;
            q = p * p + w;
            z = (float) Math.sqrt(Math.abs(q));
            H[nnn][nnn] = H[nnn][nnn] + exshift;
            H[nnn-1][nnn-1] = H[nnn-1][nnn-1] + exshift;
            x = H[nnn][nnn];
   
            // Real pair
   
            if (q >= 0) {
               if (p >= 0) {
                  z = p + z;
               } else {
                  z = p - z;
               }
               d[nnn-1] = x + z;
               d[nnn] = d[nnn-1];
               if (z != 0.0) {
                  d[nnn] = x - w / z;
               }
               e[nnn-1] = 0.0f;
               e[nnn] = 0.0f;
               x = H[nnn][nnn-1];
               s = Math.abs(x) + Math.abs(z);
               p = x / s;
               q = z / s;
               r = (float) Math.sqrt(p * p+q * q);
               p = p / r;
               q = q / r;
   
               // Row modification
   
               for (int j = nnn-1; j < nn; j++) {
                  z = H[nnn-1][j];
                  H[nnn-1][j] = q * z + p * H[nnn][j];
                  H[nnn][j] = q * H[nnn][j] - p * z;
               }
   
               // Column modification
   
               for (int i = 0; i <= nnn; i++) {
                  z = H[i][nnn-1];
                  H[i][nnn-1] = q * z + p * H[i][nnn];
                  H[i][nnn] = q * H[i][nnn] - p * z;
               }
   
               // Accumulate transformations
   
               for (int i = low; i <= high; i++) {
                  z = V[i][nnn-1];
                  V[i][nnn-1] = q * z + p * V[i][nnn];
                  V[i][nnn] = q * V[i][nnn] - p * z;
               }
   
            // Complex pair
   
            } else {
               d[nnn-1] = x + p;
               d[nnn] = x + p;
               e[nnn-1] = z;
               e[nnn] = -z;
            }
            nnn = nnn - 2;
            iter = 0;
   
         // No convergence yet
   
         } else {
   
            // Form shift
   
            x = H[nnn][nnn];
            y = 0.0f;
            w = 0.0f;
            if (l < nnn) {
               y = H[nnn-1][nnn-1];
               w = H[nnn][nnn-1] * H[nnn-1][nnn];
            }
   
            // Wilkinson's original ad hoc shift
   
            if (iter == 10) {
               exshift += x;
               for (int i = low; i <= nnn; i++) {
                  H[i][i] -= x;
               }
               s = Math.abs(H[nnn][nnn-1]) + Math.abs(H[nnn-1][nnn-2]);
               x = y = 0.75f * s;
               w = -0.4375f * s * s;
            }

            // MATLAB's new ad hoc shift

            if (iter == 30) {
                s = (y - x) / 2.0f;
                s = s * s + w;
                if (s > 0) {
                    s = (float) Math.sqrt(s);
                    if (y < x) {
                       s = -s;
                    }
                    s = x - w / ((y - x) / 2.0f + s);
                    for (int i = low; i <= nnn; i++) {
                       H[i][i] -= s;
                    }
                    exshift += s;
                    x = y = w = 0.964f;
                }
            }
   
            iter = iter + 1;   // (Could check iteration count here.)
   
            // Look for two consecutive small sub-diagonal elements
   
            int m = nnn-2;
            while (m >= l) {
               z = H[m][m];
               r = x - z;
               s = y - z;
               p = (r * s - w) / H[m+1][m] + H[m][m+1];
               q = H[m+1][m+1] - z - r - s;
               r = H[m+2][m+1];
               s = Math.abs(p) + Math.abs(q) + Math.abs(r);
               p = p / s;
               q = q / s;
               r = r / s;
               if (m == l) {
                  break;
               }
               if (Math.abs(H[m][m-1]) * (Math.abs(q) + Math.abs(r)) <
                  eps * (Math.abs(p) * (Math.abs(H[m-1][m-1]) + Math.abs(z) +
                  Math.abs(H[m+1][m+1])))) {
                     break;
               }
               m--;
            }
   
            for (int i = m+2; i <= nnn; i++) {
               H[i][i-2] = 0.0f;
               if (i > m+2) {
                  H[i][i-3] = 0.0f;
               }
            }
   
            // float QR step involving rows l:n and columns m:n
   

            for (int k = m; k <= nnn-1; k++) {
               boolean notlast = (k != nnn-1);
               if (k != m) {
                  p = H[k][k-1];
                  q = H[k+1][k-1];
                  r = (notlast ? H[k+2][k-1] : 0.0f);
                  x = Math.abs(p) + Math.abs(q) + Math.abs(r);
                  if (x == 0.0) {
                      continue;
                  }
                  p = p / x;
                  q = q / x;
                  r = r / x;
               }

               s = (float) Math.sqrt(p * p + q * q + r * r);
               if (p < 0) {
                  s = -s;
               }
               if (s != 0) {
                  if (k != m) {
                     H[k][k-1] = -s * x;
                  } else if (l != m) {
                     H[k][k-1] = -H[k][k-1];
                  }
                  p = p + s;
                  x = p / s;
                  y = q / s;
                  z = r / s;
                  q = q / p;
                  r = r / p;
   
                  // Row modification
   
                  for (int j = k; j < nn; j++) {
                     p = H[k][j] + q * H[k+1][j];
                     if (notlast) {
                        p = p + r * H[k+2][j];
                        H[k+2][j] = H[k+2][j] - p * z;
                     }
                     H[k][j] = H[k][j] - p * x;
                     H[k+1][j] = H[k+1][j] - p * y;
                  }
   
                  // Column modification
   
                  for (int i = 0; i <= Math.min(nnn,k+3); i++) {
                     p = x * H[i][k] + y * H[i][k+1];
                     if (notlast) {
                        p = p + z * H[i][k+2];
                        H[i][k+2] = H[i][k+2] - p * r;
                     }
                     H[i][k] = H[i][k] - p;
                     H[i][k+1] = H[i][k+1] - p * q;
                  }
   
                  // Accumulate transformations
   
                  for (int i = low; i <= high; i++) {
                     p = x * V[i][k] + y * V[i][k+1];
                     if (notlast) {
                        p = p + z * V[i][k+2];
                        V[i][k+2] = V[i][k+2] - p * r;
                     }
                     V[i][k] = V[i][k] - p;
                     V[i][k+1] = V[i][k+1] - p * q;
                  }
               }  // (s != 0)
            }  // k loop
         }  // check convergence
      }  // while (n >= low)
      
      // Backsubstitute to find vectors of upper triangular form

      if (norm == 0.0) {
         return;
      }
   
      for (nnn = nn-1; nnn >= 0; nnn--) {
         p = d[nnn];
         q = e[nnn];
   
         // Real vector
   
         if (q == 0) {
            int l = nnn;
            H[nnn][nnn] = 1.0f;
            for (int i = nnn-1; i >= 0; i--) {
               w = H[i][i] - p;
               r = 0.0f;
               for (int j = l; j <= nnn; j++) {
                  r = r + H[i][j] * H[j][nnn];
               }
               if (e[i] < 0.0) {
                  z = w;
                  s = r;
               } else {
                  l = i;
                  if (e[i] == 0.0) {
                     if (w != 0.0) {
                        H[i][nnn] = -r / w;
                     } else {
                        H[i][nnn] = -r / (eps * norm);
                     }
   
                  // Solve real equations
   
                  } else {
                     x = H[i][i+1];
                     y = H[i+1][i];
                     q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
                     t = (x * s - z * r) / q;
                     H[i][nnn] = t;
                     if (Math.abs(x) > Math.abs(z)) {
                        H[i+1][nnn] = (-r - w * t) / x;
                     } else {
                        H[i+1][nnn] = (-s - y * t) / z;
                     }
                  }
   
                  // Overflow control
   
                  t = Math.abs(H[i][nnn]);
                  if ((eps * t) * t > 1) {
                     for (int j = i; j <= nnn; j++) {
                        H[j][nnn] = H[j][nnn] / t;
                     }
                  }
               }
            }
   
         // Complex vector
   
         } else if (q < 0) {
            int l = nnn-1;

            // Last vector component imaginary so matrix is triangular
   
            if (Math.abs(H[nnn][nnn-1]) > Math.abs(H[nnn-1][nnn])) {
               H[nnn-1][nnn-1] = q / H[nnn][nnn-1];
               H[nnn-1][nnn] = -(H[nnn][nnn] - p) / H[nnn][nnn-1];
            } else {
               cdiv(0.0f,-H[nnn-1][nnn],H[nnn-1][nnn-1]-p,q);
               H[nnn-1][nnn-1] = cdivr;
               H[nnn-1][nnn] = cdivi;
            }
            H[nnn][nnn-1] = 0.0f;
            H[nnn][nnn] = 1.0f;
            for (int i = nnn-2; i >= 0; i--) {
               float ra,sa,vr,vi;
               ra = 0.0f;
               sa = 0.0f;
               for (int j = l; j <= nnn; j++) {
                  ra = ra + H[i][j] * H[j][nnn-1];
                  sa = sa + H[i][j] * H[j][nnn];
               }
               w = H[i][i] - p;
   
               if (e[i] < 0.0) {
                  z = w;
                  r = ra;
                  s = sa;
               } else {
                  l = i;
                  if (e[i] == 0) {
                     cdiv(-ra,-sa,w,q);
                     H[i][nnn-1] = cdivr;
                     H[i][nnn] = cdivi;
                  } else {
   
                     // Solve complex equations
   
                     x = H[i][i+1];
                     y = H[i+1][i];
                     vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
                     vi = (d[i] - p) * 2.0f * q;
                     if (vr == 0.0 & vi == 0.0) {
                        vr = eps * norm * (Math.abs(w) + Math.abs(q) +
                        Math.abs(x) + Math.abs(y) + Math.abs(z));
                     }
                     cdiv(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi);
                     H[i][nnn-1] = cdivr;
                     H[i][nnn] = cdivi;
                     if (Math.abs(x) > (Math.abs(z) + Math.abs(q))) {
                        H[i+1][nnn-1] = (-ra - w * H[i][nnn-1] + q * H[i][nnn]) / x;
                        H[i+1][nnn] = (-sa - w * H[i][nnn] - q * H[i][nnn-1]) / x;
                     } else {
                        cdiv(-r-y*H[i][nnn-1],-s-y*H[i][nnn],z,q);
                        H[i+1][nnn-1] = cdivr;
                        H[i+1][nnn] = cdivi;
                     }
                  }
   
                  // Overflow control

                  t = Math.max(Math.abs(H[i][nnn-1]),Math.abs(H[i][nnn]));
                  if ((eps * t) * t > 1) {
                     for (int j = i; j <= nnn; j++) {
                        H[j][nnn-1] = H[j][nnn-1] / t;
                        H[j][nnn] = H[j][nnn] / t;
                     }
                  }
               }
            }
         }
      }
   
      // Vectors of isolated roots
   
      for (int i = 0; i < nn; i++) {
         if (i < low | i > high) {
              System.arraycopy(H[i], i, V[i], i, nn - i);
         }
      }
   
      // Back transformation to get eigenvectors of original matrix
   
      for (int j = nn-1; j >= low; j--) {
         for (int i = low; i <= high; i++) {
            z = 0.0f;
            for (int k = low; k <= Math.min(j,high); k++) {
               z = z + V[i][k] * H[k][j];
            }
            V[i][j] = z;
         }
      }
   }


/* ------------------------
   Constructor
 * ------------------------ */

   /** Check for symmetry, then construct the eigenvalue decomposition
       Structure to access D and V.
   @param Arg    Square matrix
   */

   public EigenvalueDecomposition (Matrix Arg) {
      float[][] A = Arg.getArray();
      n = Arg.getColumnDimension();
      V = new float[n][n];
      d = new float[n];
      e = new float[n];

      issymmetric = true;
      for (int j = 0; (j < n) & issymmetric; j++) {
         for (int i = 0; (i < n) & issymmetric; i++) {
            issymmetric = (A[i][j] == A[j][i]);
         }
      }

      if (issymmetric) {
         for (int i = 0; i < n; i++) {
              System.arraycopy(A[i], 0, V[i], 0, n);
         }
   
         // Tridiagonalize.
         tred2();
   
         // Diagonalize.
         tql2();

      } else {
         H = new float[n][n];
         ort = new float[n];
         
         for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
               H[i][j] = A[i][j];
            }
         }
   
         // Reduce to Hessenberg form.
         orthes();
   
         // Reduce Hessenberg to real Schur form.
         hqr2();
      }
   }

/* ------------------------
   Public Methods
 * ------------------------ */

   /** Return the eigenvector matrix
   @return     V
   */

   public Matrix getV () {
      return new Matrix(V,n,n);
   }

   /** Return the real parts of the eigenvalues
   @return     real(diag(D))
   */

   public float[] getRealEigenvalues () {
      return d;
   }

   /** Return the imaginary parts of the eigenvalues
   @return     imag(diag(D))
   */

   public float[] getImagEigenvalues () {
      return e;
   }

   /** Return the block diagonal eigenvalue matrix
   @return     D
   */

   public Matrix getD () {
      Matrix X = new Matrix(n,n);
      float[][] D = X.getArray();
      for (int i = 0; i < n; i++) {
         for (int j = 0; j < n; j++) {
            D[i][j] = 0.0f;
         }
         D[i][i] = d[i];
         if (e[i] > 0) {
            D[i][i+1] = e[i];
         } else if (e[i] < 0) {
            D[i][i-1] = e[i];
         }
      }
      return X;
   }
  private static final long serialVersionUID = 1;
}
