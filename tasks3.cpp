#include<iostream>
#include<cmath>
#include<algorithm>
#include<fstream>
#include<string>

#ifndef MVECTOR_H // the 'include guard'
#define MVECTOR_H // see C++ Primer Sec. 2.9.2

#include <vector>

// Class that represents a mathematical vector
class MVector
{

//Operator overload for <<
friend std::ostream& operator<< (std::ostream& os, const MVector& output)
{
	os<<"(";
	for (int i=0;i<output.size()-1;i++)
	{
		os<<output[i]<<",";
	}
	os<<output[output.size()-1]<<")";
	return os;
}

public:
	// constructors
	MVector() {}
	explicit MVector(int n) : v(n) {}
	//the explicit keyword prevents implicit conversions from type integer to MVector. without the explicit keyword (and using the commented out section at the bottom) then a is a length 2 vector when added
	//to y because of the implicit conversion. Adding the two does not change the first 2 cpts, so z[0]=y[0], and z[1]=y[1] i.e. where there is an a[i] value, but  z[2] is very strange.
	MVector(int n, double x) : v(n, x) {}

	// access element (lvalue) (see example sheet 5, q5.6)
	double &operator[](int index) { return v[index]; }

	// access element (rvalue) (see example sheet 5, q5.7)
	double operator[](int index) const { return v[index]; }

	int size() const { return v.size(); } // number of elements

	double LInfNorm() const
	{
		double maxAbs = 0;
		std::size_t s = size();
		for (int i=0; i<s; i++)
		{
			maxAbs = std::max(std::abs(v[i]), maxAbs);
		}
		return maxAbs;
	}
	
	double L2Norm() const
	{
		double rootsumsquares=0;
		for (int i=0;i<size();i++)
		{
			rootsumsquares+= v[i]*v[i];
		}
		rootsumsquares=pow(rootsumsquares,0.5);
		return rootsumsquares;
	}

private:
	std::vector<double> v;
};

// Operator overload for "scalar * vector"
inline MVector operator*(const double& lhs,const MVector& rhs)
{
	MVector temp(rhs);
	for (int i=0;i<temp.size();i++) temp[i]*=lhs;
	return temp;
}
// Operator overload for "vector * scalar"
inline MVector operator*(const MVector& lhs,const double& rhs)
{
	MVector temp(lhs);
	for (int i=0;i<temp.size();i++) temp[i]*=rhs;
	return temp;
}
//Operator overload for "vector + vector"
inline MVector operator+ (const MVector& lhs, const MVector& rhs)
{
	if (lhs.size()==rhs.size())
	{
		MVector temp(lhs);
		for (int i=0; i<temp.size();i++) temp[i]+=rhs[i];
		return temp;
	}
	else 
	{
		std::cout<<"ERROR: Attempted addition of two vectors of different length."<<std::endl;
		exit(1);
	}
} 
//Operator overload for "vector - vector"
inline MVector operator- (const MVector& lhs, const MVector& rhs)
{
	if (lhs.size()==rhs.size())
	{
		MVector temp(lhs);
		for (int i=0; i<temp.size();i++) temp[i]-=rhs[i];
		return temp;
	}
	else
	{
		std::cout<<"ERROR: Attempted subtraction of two vectors of different length."<<std::endl;
		exit(1);
	}
}
//Operator overload for "vector / scalar"
inline MVector operator/ (const MVector& lhs, const double& rhs)
{
	MVector temp(lhs);
	for (int i=0;i<temp.size();i++) temp[i]/=rhs;
	return temp;
}

double dot(const MVector& lhs, const MVector& rhs)
{
	double dotproduct=0;
	if (lhs.size()==rhs.size())
	{
		for (int i=0; i<lhs.size(); i++)
		{
			dotproduct+=lhs[i]*rhs[i];
		}
		return dotproduct;
	}
	else 
	{
		std::cout<<"ERROR: Attempted dot product of two vectors of different length."<<std::endl;
		exit(1);
	}
}
#endif

#ifndef MMATRIX_H // the 'include guard'
#define MMATRIX_H

#include <vector>
#include <iostream>

// Class that represents a mathematical matrix
class MMatrix
{

friend std::ostream& operator<< (std::ostream& os, MMatrix A)
{
	for (int i=0;i<A.Rows();i++)
	{
		for (int j=0;j<A.Cols();j++)
		{
			os.precision(2);
			os.width(3);os<<A(i,j);
		}
		os<<std::endl;
	}
	return os;
}

public:
	// constructors
	MMatrix() : nRows(0), nCols(0) {}
	MMatrix(int n, int m, double x = 0) : nRows(n), nCols(m), A(n * m, x) {}

	// set all matrix entries equal to a double
	MMatrix &operator=(double x)
	{
		for (int i = 0; i < nRows * nCols; i++) A[i] = x;
		return *this;
	}

	// access element, indexed by (row, column) [rvalue]
	double operator()(int i, int j) const
	{
		return A[j + i * nCols];
	}

	// access element, indexed by (row, column) [lvalue]
	double &operator()(int i, int j)
	{
		return A[j + i * nCols];
	}

	// size of matrix
	int Rows() const { return nRows; }
	int Cols() const { return nCols; }


private:
	unsigned int nRows, nCols;
	std::vector<double> A;
};
MVector operator*(const MMatrix& A, const MVector& x)
{
	MVector b(A.Rows());
	if (A.Cols() == x.size())
	{
		double sum;
		for (int i=0;i<A.Rows();i++)
		{
			sum=0;
			for (int j=0;j<A.Cols();j++)
			{
				sum+=A(i,j)*x[j];
			}
			b[i]=sum;
		}
		return b;
	}
	else 
	{
		std::cout<<"ERROR: Attempted matrix * vector multiplication where number of matrix columns != vector length."<<std::endl;
		exit(1);
	}
}
#endif


MMatrix poissonMatrix(int n)
{
	MMatrix A(n,n);
	for (int i=0;i<n;i++)
	{
		for (int j=0;j<n;j++)
		{
			if (i==j) A(i,j)=2;
			else if (i-j==1 || j-i==1) A(i,j)=-1;
		}
	}
	return A;
}

MVector ConjGradMethod(MMatrix A, MVector x, MVector b, std::ofstream& fileName)
{
	int maxIterations = 1000;
	double tol = 1e-6;
	MVector r=b-A*x;
	MVector p=r;
	double alpha=0, beta=0;
	for (int iter=0; iter<maxIterations; iter++)
	{
		alpha=dot(r,r)/dot(p,A*p);
		x=x+alpha*p;
		MVector rkm1=r;
		r=r-alpha*(A*p);
		// ...calculate new values for x and r here...
		// check if solution is accurate enough
		if (r.L2Norm() < tol) break;
		{
			beta=dot(r,r)/dot(rkm1,rkm1);
			p=r+beta*p;
		// ...calculate new conjugate vector p here...
		}
	}
	return x;
}

int main()
{
	std::string extension(".txt");
	MVector N(3); N[0]=10; N[1]=25; N[2]=100;
	std::ofstream fileName;
	for(int i=0;i<3;i++)
	{
		int n=N[i];
		MMatrix A= poissonMatrix(n);
		MVector b(n), x0(n), r0(n);
		for (int i=0;i<n;i++) 
		{
			b[i]=1/pow(n+1,2);
			x0[i]=0;
		}
		r0=b-A*x0;
		std::string fileNametxt(std::to_string(n)+extension);
		fileName.open(fileNametxt);
		if (!fileName) return 1;
		MVector c=ConjGradMethod(A,x0,b,fileName);
		fileName<<c<<std::endl;
		fileName.close();
	}

	
	return 0;
}
