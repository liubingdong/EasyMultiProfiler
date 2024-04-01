#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace std;
using namespace RcppParallel;
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::depends(RcppParallel)]]

/* 计算矩阵内部的特征 */
// pre-compute sum and stdev
struct cor_p1 : public Worker {
	const RMatrix<double> mat;
	const int rstart, rend, nperiod;

	RVector<double> rsum, rstdev;

	cor_p1(const NumericMatrix& mat, const int rstart, const int rend,
		NumericVector rsum, NumericVector rstdev)
		: mat(mat), rstart(rstart), rend(rend), nperiod(rend - rstart), rsum(rsum), rstdev(rstdev) { }

	void operator()(size_t begin, size_t end) {
		for (size_t c = begin; c < end; c++) {
			double sum, sum2;
			sum = sum2 = 0;

			for (int r = rstart; r < rend; r++) {
				double d = mat(r,c);
				sum += d;
				sum2 += pow(d,2);
			}

			rsum[c] = sum;
			rstdev[c] = sqrt(nperiod * sum2 - pow(sum,2));
		}
	}
};


/* 根据第一、第二个矩阵的特征计算相关性和P值：用于两个矩阵时 */
// compute correlation
struct cor_p2 : public Worker {
	const RMatrix<double> mat;
	const RMatrix<double> mat2;
	const int rstart, rend, nperiod;
	const size_t mat2ncol;
	const RVector<double> sum, stdev;
	const RVector<double> sum2, stdev2;
    
	RMatrix<double> rmat;
	RMatrix<double> pmat;
	cor_p2(const NumericMatrix& mat, const NumericMatrix& mat2, const int rstart, const int rend, 
		const NumericVector& sum, const NumericVector& stdev, const NumericVector& sum2, const NumericVector& stdev2,
		NumericMatrix rmat, NumericMatrix pmat)
		: mat(mat),mat2(mat2), rstart(rstart), rend(rend),  nperiod(rend - rstart),  mat2ncol(mat2.ncol()), sum(sum), stdev(stdev), sum2(sum2), stdev2(stdev2),rmat(rmat),pmat(pmat) {}

	void operator()(size_t begin, size_t end) {
		for (size_t c1 = begin; c1 < end; c1++) {
			for (size_t c2 = 0; c2 < mat2ncol; c2++) {/* 这里为矩阵2的列数 */
				double sXY = 0;
				for (int r = rstart; r < rend; r++)
					sXY += mat(r,c1) * mat2(r,c2);
				
				rmat(c1,c2) = (nperiod * sXY - sum[c1] * sum2[c2]) / (stdev[c1] * stdev2[c2]);
				double tmpdt = (abs(rmat(c1,c2))) * (sqrt(nperiod-2));
				double tmpdt2 = rmat(c1,c2)*rmat(c1,c2);
				tmpdt2 = sqrt(1-tmpdt2);
				tmpdt2 = tmpdt/tmpdt2;
				pmat(c1,c2) = (1-(R::pt(tmpdt2,nperiod-2,true,false)))*2.0;				
			}
		}
	}
};



/* 根据第一个矩阵的特征计算相关性：用于单个矩阵时 */
struct cor_p3 : public Worker {
	const RMatrix<double> mat;
	const int rstart, rend, nperiod;
	const RVector<double> sum, stdev;
    
	RMatrix<double> rmat;
	RMatrix<double> pmat;
	cor_p3(const NumericMatrix& mat, const int rstart, const int rend,
		const NumericVector& sum, const NumericVector& stdev, 
		NumericMatrix rmat,NumericMatrix pmat)
		: mat(mat), rstart(rstart), rend(rend), nperiod(rend - rstart), sum(sum), stdev(stdev), rmat(rmat), pmat(pmat) {}

	void operator()(size_t begin, size_t end) {
		for (size_t c1 = begin; c1 < end; c1++) {
			rmat(c1,c1) = 1;
			pmat(c1,c1) = 1;
			for (size_t c2 = 0; c2 < c1; c2++) {
				double sXY = 0;
				for (int r = rstart; r < rend; r++)
					sXY += mat(r,c1) * mat(r,c2);

				rmat(c1,c2) = (nperiod * sXY - sum[c1] * sum[c2]) / (stdev[c1] * stdev[c2]);
				rmat(c2,c1) = rmat(c1,c2);
				double tmpdt = (abs(rmat(c1,c2))) * (sqrt(nperiod-2));
				double tmpdt2 = rmat(c1,c2)*rmat(c1,c2);
				tmpdt2 = sqrt(1-tmpdt2);
				tmpdt2 = tmpdt/tmpdt2;
				pmat(c1,c2) = (1-(R::pt(tmpdt2,nperiod-2,true,false)))*2.0;
				pmat(c2,c1) = pmat(c1,c2);				
			}
		}
	}
};


/* 开始计算单个矩阵的相关性 */
inline List cp_cor_s_helper(const NumericMatrix& mat, const int rstart, const int rend) {
	int nc = mat.ncol();
	NumericVector rsum(nc), rstdev(nc);

	cor_p1 p1(mat, rstart, rend, rsum, rstdev);
	parallelFor(0, nc, p1);

	NumericMatrix rmat(nc, nc);
	NumericMatrix pmat(nc, nc);
	List R_P = List::create(Named("R_matrix") = rmat, Named("P_matrix") = pmat);
	cor_p3 p2(mat, rstart, rend, rsum, rstdev, rmat, pmat);
	parallelFor(0, nc, p2);

	return R_P;
}

/* 导出函数，R可调用 */
// [[Rcpp::export]]
List cp_cor_s(NumericMatrix mat) {
	return cp_cor_s_helper(mat, 0, mat.nrow());
}

/* 开始计算两个矩阵的相关性 */
inline List cp_cor_t_helper(const NumericMatrix& mat, const NumericMatrix& mat2, const int rstart, const int rend) {
	int nc = mat.ncol();
	int nc2 = mat2.ncol();
	NumericVector rsum(nc), rstdev(nc);
	NumericVector rsum2(nc2), rstdev2(nc2);
	
	cor_p1 p1(mat, rstart, rend, rsum, rstdev);
	parallelFor(0, nc, p1);
	cor_p1 p2(mat2, rstart, rend, rsum2, rstdev2);
	parallelFor(0, nc2, p2);
	NumericMatrix rmat(nc, nc2);
	NumericMatrix pmat(nc, nc2);
	List R_P = List::create(Named("R_matrix") = rmat, Named("P_matrix") = pmat);
	cor_p2 p3(mat,mat2, rstart, rend, rsum, rstdev,rsum2, rstdev2, rmat, pmat);
	parallelFor(0, nc, p3);

	return R_P;
}

// [[Rcpp::export]]
List cp_cor_t(NumericMatrix mat,NumericMatrix mat2) {
	return cp_cor_t_helper(mat,mat2, 0, mat.nrow());
}