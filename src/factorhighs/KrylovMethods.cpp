#include "KrylovMethods.h"

#include "auxiliary/VectorOperations.h"

void IpmMatrix::reset(const HighsSparseMatrix& A,
                      const std::vector<double>& scaling, bool use_as) {
  A_ = &A;
  scaling_ = &scaling;
  use_as_ = use_as;
}

void IpmMatrix::apply(std::vector<double>& x) const {
  const int n = A_->num_col_;
  const int m = A_->num_row_;

  if (use_as_) {
    // compute [y_1; y_2] = [-Theta^-1 * x_1 + A^T * x_2; A * x_1]

    // split x
    std::vector<double> x_1(x.begin(), x.begin() + n);
    std::vector<double> x_2(x.begin() + n, x.end());

    std::vector<double> y_1(n);
    std::vector<double> y_2(m);

    // y_1 = A^T * x_2
    A_->alphaProductPlusY(1.0, x_2, y_1, true);

    // y_1 -= Theta^-1 * x_1
    for (int i = 0; i < n; ++i) y_1[i] -= (*scaling_)[i] * x_1[i];

    // y_2 = A * x_1
    A_->alphaProductPlusY(1.0, x_1, y_2);

    // put result back into x
    x.clear();
    x.insert(x.begin(), y_1.begin(), y_1.end());
    x.insert(x.begin() + n, y_2.begin(), y_2.end());

  } else {
    // compute A * Theta * A^T * x

    std::vector<double> w(n);

    // w = A^T * x
    A_->alphaProductPlusY(1.0, x, w, true);

    // w = Theta * w
    for (int i = 0; i < n; ++i) w[i] /= (*scaling_)[i];

    // x = A * w
    A_->alphaProductPlusY(1.0, w, x);
  }
}

IpmFactor::IpmFactor(const Numeric& N) : N_{N} {}
void IpmFactor::apply(std::vector<double>& x) const { N_.solve(x); }

void NeDiagPrec::reset(const HighsSparseMatrix& A,
                       const std::vector<double>& scaling) {
  // prepare diagonal preconditioner
  diag.assign(A.num_row_, 0.0);

  // build diagonal of normal equations
  for (int c = 0; c < A.num_col_; ++c) {
    for (int el = A.start_[c]; el < A.start_[c + 1]; ++el) {
      int r = A.index_[el];
      double v = A.value_[el];

      diag[r] += v * v / scaling[c];
    }
  }

  // compute inverse of diagonal entries
  for (int i = 0; i < diag.size(); ++i) diag[i] = 1.0 / diag[i];
}
void NeDiagPrec::apply(std::vector<double>& x) const {
  // apply diagonal preconditioner
  for (int i = 0; i < diag.size(); ++i) x[i] *= diag[i];
}

void applyRotation(double& x, double& y, double c, double s) {
  double t = c * x + s * y;
  y = -s * x + c * y;
  x = t;
}
void getRotation(double x, double y, double& c, double& s) {
  if (y == 0.0) {
    c = 1.0;
    s = 0.0;
  } else if (std::abs(y) > std::abs(x)) {
    double t = x / y;
    s = 1.0 / sqrt(1.0 + t * t);
    c = t * s;
  } else {
    double t = y / x;
    c = 1.0 / sqrt(1.0 + t * t);
    s = t * c;
  }
}
void update(std::vector<double>& x, int k, std::vector<double>& H, int ldh,
            std::vector<double>& s, std::vector<std::vector<double>>& V) {
  // Solve H * y = s
  std::vector<double> y = s;
  for (int i = k; i >= 0; --i) {
    y[i] /= H[i + ldh * i];
    for (int j = i - 1; j >= 0; --j) {
      y[j] -= H[j + ldh * i] * y[i];
    }
  }
  // x += V * y
  for (int j = 0; j <= k; ++j) vectorAdd(x, V[j], y[j]);
}

int Gmres(const AbstractMatrix* M, const AbstractMatrix* P,
          const std::vector<double>& b, std::vector<double>& x, double tol,
          int maxit) {
  // Attempt to solve M * x = b using GMRES, with the preconditioner P.
  // Return the number of iterations taken.
  //
  // Based on Netlib implementation and Kelley "Iterative Methods for Linear and
  // Nonlinear Equations":
  // https://www.netlib.org/templates/cpp/gmres.h
  // https://www.netlib.org/templates/templates.pdf
  // https://en.wikipedia.org/wiki/Generalized_minimal_residual_method#Example_code

  // sizes
  const int n = x.size();
  const int m = maxit;

  // vectors
  std::vector<std::vector<double>> V;

  // Hessenberg matrix
  std::vector<double> H((m + 1) * m);
  const int ldh = m + 1;

  // Givens rotations
  std::vector<double> sn(m + 1);
  std::vector<double> cs(m + 1);

  // workspace
  std::vector<double> w(n);

  // preconditioned residual P^-1 * (b - M * x)
  std::vector<double> r = x;
  M->apply(r);
  vectorAdd(r, b, 1.0, -1.0);
  if (P) P->apply(r);
  double beta = norm2(r);

  // preconditioned rhs P^-1 * b
  w = b;
  if (P) P->apply(w);
  double normb = norm2(w);

  if (normb == 0.0) normb = 1.0;
  if (beta <= tol) return 0;

  V.push_back(r);
  vectorScale(V[0], 1.0 / beta);

  // s = beta * e1
  std::vector<double> s(m + 1);
  s[0] = beta;

  // main loop
  int i;
  for (i = 0; i < maxit; ++i) {
    // Arnoldi
    w = V[i];
    M->apply(w);
    if (P) P->apply(w);
    V.push_back(w);
    for (int k = 0; k <= i; ++k) {
      H[k + ldh * i] = dotProd(V.back(), V[k]);
      vectorAdd(V.back(), V[k], -H[k + ldh * i]);
    }
    H[i + 1 + ldh * i] = norm2(V.back());

    // re-orthogonalize
    // if (norm2(V.back()) / norm2(w) < 1e-6) {
    // printf("Re-orthogonalize============\n");
    for (int k = 0; k <= i; ++k) {
      double temp = dotProd(V.back(), V[k]);
      H[k + ldh * i] += temp;
      vectorAdd(V.back(), V[k], -temp);
    }
    H[i + 1 + ldh * i] = norm2(V.back());
    // }

    // normalize
    vectorScale(V.back(), 1.0 / H[i + 1 + ldh * i]);

    // apply Givens rotations
    for (int k = 0; k < i; ++k)
      applyRotation(H[k + ldh * i], H[k + 1 + ldh * i], cs[k], sn[k]);

    // get latest rotation and apply it also to s
    getRotation(H[i + ldh * i], H[i + 1 + ldh * i], cs[i], sn[i]);
    applyRotation(H[i + ldh * i], H[i + 1 + ldh * i], cs[i], sn[i]);
    applyRotation(s[i], s[i + 1], cs[i], sn[i]);

    // printf("%d: %e\n", i, std::abs(s[i + 1]));

    // check termination
    if (std::abs(s[i + 1]) < tol) break;
  }

  // compute solution
  update(x, i, H, ldh, s, V);

  return i + 1;
}

int Cg(const AbstractMatrix* M, const AbstractMatrix* P,
       const std::vector<double>& b, std::vector<double>& x, double tol,
       int maxit) {
  // Attempt to solve M * x = b using CG, with the preconditioner P.
  // Return the number of iterations taken.

  int n = b.size();

  std::vector<double> w(n);

  std::vector<double> r = x;
  M->apply(r);
  vectorAdd(r, b, 1.0, -1.0);

  std::vector<double> z = r;
  if (P) P->apply(z);

  std::vector<double> p(z);

  double rho_old = dotProd(r, z);
  double rho_new;
  double norm_b = norm2(b);

  int iter = 0;
  while (iter < maxit) {
    w = p;
    M->apply(w);
    double alpha = rho_old / dotProd(p, w);
    vectorAdd(x, p, alpha);
    vectorAdd(r, w, -alpha);

    z = r;
    if (P) P->apply(z);
    rho_new = dotProd(r, z);

    if (norm2(r) < tol * norm_b) break;

    vectorAdd(p, z, 1.0, rho_new / rho_old);
    rho_old = rho_new;
    ++iter;

    if (isNanVector(x)) {
      printf("CG: x is nan at iter %d\n", iter);
      break;
    }
  }

  return iter;
}