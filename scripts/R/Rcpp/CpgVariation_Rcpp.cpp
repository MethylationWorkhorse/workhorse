#include <RcppArmadillo.h>
#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <math.h>
using namespace std;
//[[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// Function for calculating median 
double findMedian(double a[], int n) 
{ 
  // First we sort the array 
  sort(a, a + n);
  
  // check for even case 
  if (n % 2 != 0) 
    return (double)a[n / 2]; 
  
  return (double)(a[(n - 1) / 2] + a[n / 2]) / 2.0; 
}

double sum(const vector<double>& a)
{
  double s = 0;
  for (int i = 0; i < a.size(); i++)
  {
    s += a[i];
  }
  return s;
}

double mean(const vector<double>& a)
{
  return sum(a) / a.size();
}

double sqsum(const vector<double>& a)
{
  double s = 0;
  for (int i = 0; i < a.size(); i++)
  {
    s += pow(a[i], 2);
  }
  return s;
}

double stdev(const vector<double>& nums)
{
  double N = nums.size();
  return pow(sqsum(nums) / N - pow(sum(nums) / N, 2), 0.5);
}

vector<double> operator-(const vector<double>& a, double b)
{
  vector<double> retvect;
  for (int i = 0; i < a.size(); i++)
  {
    retvect.push_back(a[i] - b);
  }
  return retvect;
}

vector<double> operator*(const vector<double>& a, const vector<double>& b)
{
  vector<double> retvect;
  for (int i = 0; i < a.size() ; i++)
  {
    retvect.push_back(a[i] * b[i]);
  }
  return retvect;
}

// [[Rcpp::export]]
double pearsoncoeff(const vector<double>& X, const vector<double>& Y)
{
  return sum((X - mean(X))*(Y - mean(Y))) / (X.size()*stdev(X)* stdev(Y));
}

// [[Rcpp::export]]
int C_n_choose_m(int n, int r)
{
  if (r == 0) return 1;
  
  // Extra computation saving for large R,
  // using property:
  // N choose R = N choose (N-R)
  if (r > n / 2) return C_n_choose_m(n, n - r); 
  
  int res = 1; 
  // long res = 1; 
  
  for (int k = 1; k <= r; ++k)
  {
    res *= n - k + 1;
    res /= k;
  }
  return res;
}

// [[Rcpp::export]]
NumericMatrix na_matrix(int n, int m){
  NumericMatrix mat(n,m) ;
  std::fill( mat.begin(), mat.end(), NumericVector::get_na() ) ;
  return mat ;
}

/*
 * New Code Start::
 * 
 * Needs::
 *  1. Diff between two matricies
 *  2. 
 *     
 */

// [[Rcpp::export]]
NumericMatrix C_lociRSquared(NumericMatrix x, NumericVector idxVec, CharacterVector samVec, 
                             Nullable<CharacterVector> cpgVec=R_NilValue, 
                             double minDelta2=0.2, double minDelta1=0.1, double minDelta0=0.05,
                             int verbose=0, int vt=1) {
  int dat_nrow = x.nrow(), dat_ncol = x.ncol();
  int out_nrows = 12;
  int sam_ncols = idxVec.size();
  
  // Initialize Numeric Output Matrix
  // NumericMatrix out_mat(dat_nrow, (out_nrows*sam_ncols)+1);
  NumericMatrix out_mat = na_matrix(dat_nrow, (out_nrows*sam_ncols)+1);
  vector<string> out_cols_vec;
  out_cols_vec.push_back("RSquared");
  
  if (verbose>=vt) Rcout << "Matrix Dim: nrow=" << dat_nrow << ", ncol=" << dat_ncol << "\n";
  int sampleCnt = idxVec.size();
  
  for (int loci_idx = 0; loci_idx < dat_nrow; loci_idx++) {
    vector<double> x_vec;   // Pair wise first value (no NA's)
    vector<double> y_vec;   // Pair wise second value (no NA's)
    double r2 = 0;
    
    int col_int_idx = 0; // Column Intialization Index
    for (int sIdx = 0; sIdx < sampleCnt; sIdx++) {
      int col_beg_idx = col_int_idx;
      int col_end_idx = col_beg_idx + idxVec[sIdx];
      string sampleName = (string)samVec[sIdx];
      
      vector<double> val_vec;   // The actual values (no NA's)
      vector<double> dbv_vec;  // The delta betas (no NA's)
      
      int tot_cnt = 0;
      int nan_cnt = 0;
      int cmp_cnt = 0;
      
      int db2_cnt = 0;
      int db1_cnt = 0;
      int db0_cnt = 0;
      
      if (verbose>=vt) 
        Rcout << "loci=" << loci_idx << ", sIdx=" << sIdx << ", sample=" << sampleName << ", beg=" << col_beg_idx << ", end=" << col_end_idx << "\n";
      for (int col_idxA = col_beg_idx; col_idxA < col_end_idx; col_idxA++) {
        tot_cnt++;
        if (NumericMatrix::is_na(x(loci_idx,col_idxA))) {
          nan_cnt++;
        } else {
          val_vec.push_back(x(loci_idx,col_idxA));
        }
        for (int col_idxB = col_idxA+1; col_idxB < col_end_idx; col_idxB++) {
          
          if (NumericMatrix::is_na(x(loci_idx,col_idxA))) {
          } else if (NumericMatrix::is_na(x(loci_idx,col_idxB))) {
          } else {
            cmp_cnt++;
            
            double valA = x(loci_idx, col_idxA);
            double valB = x(loci_idx, col_idxB);
            
            double db = fabs( valA - valB );
            dbv_vec.push_back(db);
            if (db > minDelta2) db2_cnt++;
            if (db > minDelta1) db1_cnt++;
            if (db > minDelta0) db0_cnt++;
            
            x_vec.push_back(valA);
            y_vec.push_back(valB);
            
            if (verbose>=vt) 
              Rcout << "\tidxA=" << col_idxA << ", idxB=" << col_idxB << ", db=" << db << ", valA=" << valA << ", valB=" << valB << "\n";
          }
        }
      }
      
      // Calculate value(mean, median, and SD)
      double val_avg = NA_REAL;
      double val_med = NA_REAL;
      double val_std = NA_REAL;
      if (val_vec.size() != 0) {
        val_avg = mean(val_vec);
        val_med = findMedian(&val_vec[0], val_vec.size());
        val_std = stdev(val_vec);
      }
      
      // Add delta-value(Mean, Median and SD)
      double dbv_avg = NA_REAL;
      double dbv_med = NA_REAL;
      double dbv_std = NA_REAL;
      
      if (dbv_vec.size() != 0) {
        dbv_avg = mean(dbv_vec);
        dbv_med = findMedian(&dbv_vec[0], dbv_vec.size());
        dbv_std = stdev(dbv_vec);
      }
      
      // Total, NA, Comparison Counts::
      int out_idx = sIdx*out_nrows;
      out_mat(loci_idx,out_idx+1) = tot_cnt;
      out_mat(loci_idx,out_idx+2) = nan_cnt;
      out_mat(loci_idx,out_idx+3) = cmp_cnt;
      
      // Value(Beta) Average, Median, SD::
      out_mat(loci_idx,out_idx+4) = val_avg;
      out_mat(loci_idx,out_idx+5) = val_med;
      out_mat(loci_idx,out_idx+6) = val_std;
      
      // Delta-Beta Counts::
      out_mat(loci_idx,out_idx+7) = db2_cnt;
      out_mat(loci_idx,out_idx+8) = db1_cnt;
      out_mat(loci_idx,out_idx+9) = db0_cnt;
      
      // Delta Beta Average, Median, SD::
      out_mat(loci_idx,out_idx+10) = dbv_avg;
      out_mat(loci_idx,out_idx+11) = dbv_med;
      out_mat(loci_idx,out_idx+12) = dbv_std;
      
      if (loci_idx==0) {
        out_cols_vec.push_back(sampleName+"_Tot_Count");
        out_cols_vec.push_back(sampleName+"_Nan_Count");
        out_cols_vec.push_back(sampleName+"_Cmp_Count");
        
        out_cols_vec.push_back(sampleName+"_avg");
        out_cols_vec.push_back(sampleName+"_med");
        out_cols_vec.push_back(sampleName+"_sd");
        
        out_cols_vec.push_back(sampleName+"_db2_Count");
        out_cols_vec.push_back(sampleName+"_db1_Count");
        out_cols_vec.push_back(sampleName+"_db0_Count");
        
        out_cols_vec.push_back(sampleName+"_db_avg");
        out_cols_vec.push_back(sampleName+"_db_med");
        out_cols_vec.push_back(sampleName+"_db_sd");
      }
      
      col_int_idx = col_end_idx;
    }
    if (x_vec.size() > 1) r2 = pearsoncoeff(x_vec, y_vec);
    if (verbose>=vt) Rcout << "R2=" << r2 << ", out_cols_len" << out_cols_vec.size() << "\n\n";
    
    out_mat(loci_idx,0) = r2;
  }
  colnames(out_mat) = wrap(out_cols_vec);
  rownames(out_mat) = rownames(x);
  // if (cpgVec.isNotNull()) rownames(out_mat) = cpgVec;
  
  return out_mat;
}

// [[Rcpp::export]]
NumericMatrix C_lociRowsToCols(NumericMatrix x, NumericVector idxVec, CharacterVector samVec, 
                               Nullable<CharacterVector> cpgVec=R_NilValue, //CharacterVector cpgVec=R_NilValue,
                               double minDelta2=0.2, double minDelta1=0.1, double minDelta0=0.05,
                               int verbose=0, int vt=1) {
  int dat_nrow = x.nrow(), dat_ncol = x.ncol();
  if (verbose>=vt) Rcout << "Inp Matrix Dim: dat_nrow=" << dat_nrow << ", dat_ncol=" << dat_ncol << "\n";

  // Two columns for each row in output
  int out_ncols = dat_nrow * 2;
  if (verbose>=vt) Rcout << "Out Matrix Dim: " << "out_ncols=" << out_ncols << "\n";
  
  // To calculate out_nrows we need to calculate the number of sample combinations (choose function)
  //  for each sample. We'll create the output column name vector at the same time.
  int out_nrows = 0;
  vector<string> out_rows_vec;
  for (int ii=0; ii < idxVec.size(); ii++) {
    int choose_cnt = C_n_choose_m(idxVec[ii], 2);
    out_nrows += choose_cnt;
    for (int jj=0; jj < choose_cnt; jj++) {
      out_rows_vec.push_back((string)samVec[ii]);
      if (verbose>=vt) Rcout << "Out Matrix Cols: ii/jj=(" << ii << "/" << jj << "), SampleName=" << samVec[ii] << "\n";
    }
  }
  if (verbose>=vt) Rcout << "Out Matrix Dim: " << "out_nrows=" << out_nrows << "\n";
  
  // Intialize output matrix and set rownames(samples)::
  // NumericMatrix out_mat(dat_nrow, out_ncols, NA);
  // NumericMatrix out_mat(dat_nrow, out_ncols);
  NumericMatrix out_mat = na_matrix(out_nrows, out_ncols);
  rownames(out_mat) = wrap(out_rows_vec);
  
  // if (1) return(out_mat);
  
  CharacterVector inp_cpg_vev = rownames(x);
  vector<string> out_cols_vec;

  if (verbose>=vt) Rcout << "Matrix Dim: nrow=" << dat_nrow << ", ncol=" << dat_ncol << "\n";
  int sampleCnt = idxVec.size();
  
  for (int loci_idx = 0; loci_idx < dat_nrow; loci_idx++) {
    int row_out_idx = 0;
    int col_out_idx = loci_idx * 2;
    
    // vector<double> x_vec;   // Pair wise first value (no NA's)
    // vector<double> y_vec;   // Pair wise second value (no NA's)
    
    // Add cpgs to column names::
    string cpgA = (string)inp_cpg_vev(loci_idx) + "_A";
    string cpgB = (string)inp_cpg_vev(loci_idx) + "_B";
    out_cols_vec.push_back(cpgA);
    out_cols_vec.push_back(cpgB);
    
    int col_int_idx = 0; // Column Intialization Index
    for (int sIdx = 0; sIdx < sampleCnt; sIdx++) {
      int col_beg_idx = col_int_idx;
      int col_end_idx = col_beg_idx + idxVec[sIdx];
      string sampleName = (string)samVec[sIdx];
      
      if (verbose>=vt) 
        Rcout << "loci=" << loci_idx << ", sIdx=" << sIdx << ", sample=" << sampleName << ", beg=" << col_beg_idx << ", end=" << col_end_idx << "\n";
      for (int col_idxA = col_beg_idx; col_idxA < col_end_idx; col_idxA++) {
        for (int col_idxB = col_idxA+1; col_idxB < col_end_idx; col_idxB++) {
          
          double valA = x(loci_idx, col_idxA);
          double valB = x(loci_idx, col_idxB);
          
          out_mat(row_out_idx, col_out_idx+0) = valA; // x(loci_idx, col_idxA);
          out_mat(row_out_idx, col_out_idx+1) = valB; // x(loci_idx, col_idxB);
          
          if (verbose>=vt) 
            Rcout << "row/col_out_idx=(" << row_out_idx << "," << col_out_idx << "), vals(A/B)=(" << valA << "," << valB << "), " <<
              "\tidxA=" << col_idxA << ", idxB=" << col_idxB << ", valA=" << valA << ", valB=" << valB << "\n";
          
          row_out_idx++;
        }
      }
      col_int_idx = col_end_idx;
    }
  }
  colnames(out_mat) = wrap(out_cols_vec);

  return out_mat;
}

// [[Rcpp::export]]
NumericMatrix C_lociVariation(NumericMatrix x, double minDelta) {
  int dat_nrow = x.nrow(), dat_ncol = x.ncol();
  int out_nrows = 5;
  // int c = NCR_C(ncol, 2);
  // double R::choose(double n, double k);
  NumericMatrix out(dat_nrow, out_nrows);
  
  // NumericMatrix out(nrow, ncol);
  // Rcout << "Matrix Dim: nrow=" << nrow << ", ncol=" << ncol << "\n";
  NumericVector v =
    NumericVector::create(1, NA_REAL, R_NaN, R_PosInf, R_NegInf);
  
  // int na_idx = ncol;
  int idxA = 0;
  int idxB = idxA + 1;

  for (int i = 0; i < dat_nrow; i++) {
    // x(i, _ ) = na_omit(x(i, _ ));
    int ii = 0;
    
    // First count the number of na values
    double nan_cnt = 0;
    double mat_cnt = 0;
    double cmp_cnt = 0;
    double dbf_cnt = 0;
    double r2 = 0;
    
    for (int j = idxA; j < idxB; j++) {
      if (NumericMatrix::is_na(x(i,j))) {
        nan_cnt++;
      }
    }
    mat_cnt = dat_ncol - nan_cnt;
    cmp_cnt = C_n_choose_m(mat_cnt, 2);
    
    vector<double> x_vec;
    vector<double> y_vec;
    
    for (int j = 0; j < dat_ncol; j++) {
      for (int k = j+1; k < dat_ncol; k++) {
        if (NumericMatrix::is_na(x(i,j))) {
          // na_idx--;
          // Rcout << "\tna_idx=" << na_idx << ", val=" << x(i,j) << "\n";
          // out(i,na_idx) = R_NaN;
          
          // out(i,na_idx) = x(i,j)-0.00001;
          // out(i,ii) = x(i,j);
        } else if (NumericMatrix::is_na(x(i,k))) {
          // na_idx--;
          // out(i,na_idx) = R_NaN;
          
          // out(i,na_idx) = x(i,k)-0.00002;
          // out(i,ii) = x(i,k);
        } else {
          // out(i,ii) = x(i,j)-x(i,k);
          
          double db = fabs( x(i,j)-x(i,k) );
          if (db > minDelta) dbf_cnt++;
          
          x_vec.push_back(x(i,j));
          y_vec.push_back(x(i,k));
          
          ii++;
          Rcout << "\ti=" << i << ", j=" << j << "k=" << k << ", val[i,j]=" << x(i,j) << ", val[i,k]=" << x(i,k) << ", db=" << db << "\n";
        }
      }
    }
    Rcout << "\ti=" << i << ", x.size=" << x_vec.size() << ", y.size=" << y_vec.size() << "\n";
    for (int xi = 0; xi < x_vec.size(); xi++) {
      Rcout << "\txi=" << xi << ", [" << x_vec[xi] << ", " << y_vec[xi] << "]\n";
    }
    
    if (mat_cnt > 1) r2 = pearsoncoeff(x_vec, y_vec);
    out(i,0) = nan_cnt;
    out(i,1) = dbf_cnt;
    out(i,2) = cmp_cnt;
    out(i,3) = r2;
    out(i,4) = dat_ncol;
  }
  return out;
}



/*
 * Usage::
 * 
 * UHM(0.02, as.matrix(sam_tib$U), as.matrix(sam_tib$H), as.matrix(sam_tib$M))
 *     
 */

// [[Rcpp::export]]
NumericMatrix UHM(double pval_min_tit,
                  double pval_min_cel,
                  double db_cutoff,
                  NumericVector LenVec,
                  NumericMatrix DatMat) {
  
  int m_nrow = DatMat.nrow();
  int l_size = LenVec.size();
  int ret_data_size = 15;
  bool verbose = false;

  Rcout << "\t[UHM]: m_nrow=" << m_nrow << ", " <<
    "l_size=" << l_size <<
    "\n";
  
  int output_fields = (l_size*ret_data_size) + 13;
  NumericMatrix UHM(m_nrow, output_fields);
  
  for (int cgn_idx=0; cgn_idx < m_nrow; cgn_idx++) {
    int ii=0;
    
    double maxM = -1, minU = 1;
    
    
    if (verbose) {
      Rcout << "cgn_idx=" << cgn_idx << ", init(ii)=" << ii
            << ", l_size=" << l_size << "\n";
    }
    for (int ll=0; ll < l_size; ll++) {
      double pval_min = pval_min_cel;
      if (ll < 4) pval_min = pval_min_tit;

      vector<double> pass_beta_vec;
      // vector<double> pass_pval_vec;
      // vector<double> fail_beta_vec;
      // vector<double> fail_pval_vec;
      
      int data_len = LenVec[ll];
      int pass_pval_len = 0;

      double pass_beta_min = -1, pass_beta_max = -1, pass_beta_med = -1, pass_pval_max = -1;
      double fail_beta_min = -1, fail_beta_max = -1, fail_pval_min = -1, fail_pval_max = -1;
      // double fail_beta_med = -1, pass_beta_avg = -1
      // int fail_beta_len;
      
      int pp_pval_dbeta_cnt = 0;
      int pf_pval_dbeta_cnt = 0;
      int fp_pval_dbeta_cnt = 0;
      int ff_pval_dbeta_cnt = 0;
      for (int jj=0; jj < data_len; jj+=2) {
        double bval = DatMat(cgn_idx, ii);
        double pval = DatMat(cgn_idx, ii+1);
        // bval = printf("%.3f", DatMat(cgn_idx,ii));
        // pval = printf("%.6f", DatMat(cgn_idx,ii+1));
        if (verbose) {
          Rcout << "\t\tjj=" << jj << ", ii=" << ii
                << ", data_len=" << data_len
                << ", pval=" << pval
                << ", bval=" << bval
                << "\n";
        }
        
        bool pass_pval1 = false;
        if (pval <= pval_min) {
          pass_pval1 = true;
          pass_beta_vec.push_back(bval);
          // pass_pval_vec.push_back(pval);
          if (pval > pass_pval_max) pass_pval_max = pval;
          if (bval > maxM) maxM = bval;
          if (bval < minU) minU = bval;
        } else {
          if (fail_beta_min==-1 || fail_beta_min > bval) fail_beta_min = bval;
          if (fail_beta_max==-1 || fail_beta_max < bval) fail_beta_max = bval;
          if (fail_pval_min==-1 || fail_pval_min > pval) fail_pval_min = pval;
          if (fail_pval_max==-1 || fail_pval_max < pval) fail_pval_max = pval;

          // fail_beta_vec.push_back(bval);
          // fail_pval_vec.push_back(pval);
        }
        int dd = ii+2;
        for (int kk=jj+2; kk < data_len; kk+=2) {
          // double bval2 = DatMat(cgn_idx, kk+ii);
          // double pval2 = DatMat(cgn_idx, kk+ii+1);
          double bval2 = DatMat(cgn_idx, dd);
          double pval2 = DatMat(cgn_idx, dd+1);
          double delta = fabs(bval-bval2);
          
          bool pass_delta = false;
          bool pass_pval2 = false;
          if (delta <= db_cutoff) pass_delta = true;
          if (pval2 <= pval_min) pass_pval2 = true;
          
          char v1='F';
          char v2='F';
          if (pass_pval1 && pass_pval2) {
            v1='P';
            if (pass_delta) {
              v2='P';
              pp_pval_dbeta_cnt++;
            } else {
              pf_pval_dbeta_cnt++;
            }
          } else {
            if (pass_delta) {
              v2='P';
              fp_pval_dbeta_cnt++;
            } else {
              ff_pval_dbeta_cnt++;
            }
          }
          if (verbose) {
            Rcout << "\t\t" << v1 << v2 << " DB: delta=" << delta << " > " << db_cutoff
                  << ", dd=" << dd
                  << ", jj=" << jj
                  << ", bval1=" << bval
                  << ", pval1=" << pval
                  << ", kk=" << kk
                  << ", bval2=" << bval2
                  << ", pval2=" << pval2
                  << "\n";
          }
          dd += 2;
        }
        ii+=2;
      }
      pass_pval_len = pass_beta_vec.size();
      if (pass_pval_len>1) {
        int mid = (int)(pass_pval_len/2);
        sort(pass_beta_vec.begin(), pass_beta_vec.end());
        pass_beta_min = pass_beta_vec[0];
        pass_beta_max = pass_beta_vec[pass_pval_len-1];
        pass_beta_med = pass_beta_vec[mid];
        // pass_beta_avg = pass_beta_vec[0];
      } else if (pass_pval_len==1) {
        pass_beta_min = pass_beta_vec[0];
        pass_beta_max = pass_beta_vec[0];
        pass_beta_med = pass_beta_vec[0];
        // pass_beta_avg = pass_beta_vec[0];
      } else if (pass_pval_len==0) {
        pass_beta_min = -1;
        pass_beta_max = -1;
        pass_beta_med = -1;
        // pass_beta_avg = -1;
      } else {
        // HOW ARE WE HERE???
        Rcerr << "\n\n\t[UHM]: HOW ARE WE HERE\n\n";
      }
      if (verbose) {
        Rcout << "\tll=" << ll << ", ii=" << ii
              << ", data_len=" <<  data_len/2
              << "\n"
              << "\tll=" << ll << ", ii=" << ii
              << ", pass_pval_len=" <<  pass_pval_len
              << ", pass_beta_min=" <<  pass_beta_min
              << ", pass_beta_med=" <<  pass_beta_med
              << ", pass_beta_max=" <<  pass_beta_max
              << "\n"
              << "\tll=" << ll << ", ii=" << ii
              << ", fail_beta_min=" <<  fail_beta_min
              << ", fail_beta_max=" <<  fail_beta_max
              << ", fail_pval_min=" <<  fail_pval_min
              << ", pass_pval_max=" <<  pass_pval_max
              << "\n";
      }
      
      int idx_offset = (ret_data_size * ll) - 1;
      UHM(cgn_idx, ++idx_offset) = data_len/2;
      UHM(cgn_idx, ++idx_offset) = pass_pval_len;

      UHM(cgn_idx, ++idx_offset) = pass_beta_min;
      UHM(cgn_idx, ++idx_offset) = pass_beta_med;
      UHM(cgn_idx, ++idx_offset) = pass_beta_max;
      UHM(cgn_idx, ++idx_offset) = pass_pval_max;

      // Failure Stats::
      UHM(cgn_idx, ++idx_offset) = fail_beta_min;
      UHM(cgn_idx, ++idx_offset) = fail_beta_max;
      UHM(cgn_idx, ++idx_offset) = fail_pval_min;
      UHM(cgn_idx, ++idx_offset) = fail_pval_max;
      
      // Delta Beta Stats::
      UHM(cgn_idx, ++idx_offset) = pp_pval_dbeta_cnt;
      UHM(cgn_idx, ++idx_offset) = pf_pval_dbeta_cnt;
      UHM(cgn_idx, ++idx_offset) = fp_pval_dbeta_cnt;
      UHM(cgn_idx, ++idx_offset) = ff_pval_dbeta_cnt;
      
      // Delta Beta Passing Percent
      int pp_perc = 0;
      int pn_pval_dbeta_cnt = pp_pval_dbeta_cnt + pf_pval_dbeta_cnt;
      // fp_pval_dbeta_cnt + ff_pval_dbeta_cnt) );
      if (pn_pval_dbeta_cnt != 0)
        pp_perc = (100 * pp_pval_dbeta_cnt / pn_pval_dbeta_cnt);
      UHM(cgn_idx, ++idx_offset) = pp_perc;

      pass_beta_vec.clear();
      // fail_beta_vec.clear();
      // Rcout << "\n\n";
    }
    
    /*
     * Cross Sample Validation::
     * 
     * 
     */
    int idx_offSetL = (ret_data_size * l_size) - 1;
    int idx_offsetU = ret_data_size * 0;
    int idx_offsetH = ret_data_size * 1;
    int idx_offsetM = ret_data_size * 2;

    double HU_MAX = UHM(cgn_idx, idx_offsetH + 4) - UHM(cgn_idx, idx_offsetU + 2);
    double HU_MED = UHM(cgn_idx, idx_offsetH + 3) - UHM(cgn_idx, idx_offsetU + 3);
    double HU_MIN = UHM(cgn_idx, idx_offsetH + 2) - UHM(cgn_idx, idx_offsetU + 4);

    double MH_MAX = UHM(cgn_idx, idx_offsetM + 4) - UHM(cgn_idx, idx_offsetH + 2);
    double MH_MED = UHM(cgn_idx, idx_offsetM + 3) - UHM(cgn_idx, idx_offsetH + 3);
    double MH_MIN = UHM(cgn_idx, idx_offsetM + 2) - UHM(cgn_idx, idx_offsetH + 4);

    double MU_MAX = UHM(cgn_idx, idx_offsetM + 4) - UHM(cgn_idx, idx_offsetU + 2);
    double MU_MED = UHM(cgn_idx, idx_offsetM + 3) - UHM(cgn_idx, idx_offsetU + 3);
    double MU_MIN = UHM(cgn_idx, idx_offsetM + 2) - UHM(cgn_idx, idx_offsetU + 4);
    
    double MH_MU_MAX = min(MH_MAX, HU_MAX);
    double MH_MU_MED = min(MH_MED, HU_MED);
    double MH_MU_MIN = min(MH_MIN, HU_MIN);

    double MU_ALL = maxM - minU;

    if (verbose) {
      Rcout << "Indexes:\n" 
            << "\tidx_offSetL=" << idx_offSetL << "\n"
            << "\tidx_offsetU=" << idx_offsetU << "\n"
            << "\tidx_offsetH=" << idx_offsetH << "\n"
            << "\tidx_offsetM=" << idx_offsetM << "\n"
            << "\n";
      Rcout << "\tMU_ALL: maxM - minU: " << maxM << " - " << minU << " == " 
            << maxM - minU << " == " << MU_ALL << "\n\n";
      Rcout << "\tMU_MAX: " << MU_MAX << " == "
            << UHM(cgn_idx, idx_offsetM + 4) << " - " 
            << UHM(cgn_idx, idx_offsetU + 2) << "\n";
      Rcout << "\tMU_MED: " << MU_MED << " == "
            << UHM(cgn_idx, idx_offsetM + 3) << " - " 
            << UHM(cgn_idx, idx_offsetU + 3) << "\n";
      Rcout << "\tMU_MIN: " << MU_MIN << " == "
            << UHM(cgn_idx, idx_offsetM + 2) << " - " 
            << UHM(cgn_idx, idx_offsetU + 4) << "\n";
    }
    
    UHM(cgn_idx, ++idx_offSetL)  = HU_MAX;
    UHM(cgn_idx, ++idx_offSetL)  = HU_MED;
    UHM(cgn_idx, ++idx_offSetL)  = HU_MIN;
    
    UHM(cgn_idx, ++idx_offSetL)  = MH_MAX;
    UHM(cgn_idx, ++idx_offSetL)  = MH_MED;
    UHM(cgn_idx, ++idx_offSetL)  = MH_MIN;

    UHM(cgn_idx, ++idx_offSetL)  = MU_MAX;
    UHM(cgn_idx, ++idx_offSetL) = MU_MED;
    UHM(cgn_idx, ++idx_offSetL) = MU_MIN;
    
    UHM(cgn_idx, ++idx_offSetL)  = MH_MU_MAX;
    UHM(cgn_idx, ++idx_offSetL)  = MH_MU_MED;
    UHM(cgn_idx, ++idx_offSetL)  = MH_MU_MIN;

    UHM(cgn_idx, ++idx_offSetL) = MU_ALL;
    
    // Rcout << "\n\n";
  }
  return(UHM);
}

/*
 * Original Unused Functions::
 */

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

// [[Rcpp::export]]
NumericVector rowSums_C(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow);
  // Rcout << "Matrix Dim: nrow=" << nrow << ", ncol=" << ncol << "\n";
  
  for (int i = 0; i < nrow; i++) {
    double total = 0;
    for (int j = 0; j < ncol; j++) {
      if (!NumericMatrix::is_na(x(i,j))) {
        // total += x(i, j);
        total++;
      }
    }
    out[i] = total;
    // out[i] = 1;
  }
  return out;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)

set.seed(1014)
x <- matrix(sample(100), 10)
# rowSums(x)
# rowSumsC(x)
# rowStats_C(x)
# evalCpp("2+2")


*/
