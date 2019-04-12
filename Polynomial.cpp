//============================================================================
// Name        : Polynomial.cpp
// Author      : Christine Trac
//============================================================================

#include <iostream>
#include "Polynomial.h"

#ifndef MARMOSET_TESTING
int main();
#endif

#ifndef MARMOSET_TESTING

int main(){
	double array1[3] = {2, 3, 6};
	poly_t poly1{nullptr, 0};
	init_poly(poly1,array1,2);

//	poly_multiply(poly1,poly2);

//	poly_divide(poly1, 5);

//	poly_diff(poly2);

	poly_approx_int(poly1, 1, 4, 2);
	std::cout << poly_approx_int(poly1, 1, 4, 2) << std::endl;

	return 0;
}

#endif

void init_poly(poly_t &p, double const init_coeffs[], unsigned int const init_degree){
	if(p.a_coeffs != nullptr){
		delete[] p.a_coeffs;
		p.a_coeffs = nullptr;
	}
	p.degree = init_degree;
	p.a_coeffs = new double[p.degree + 1];
	for(unsigned int i = 0; i<=p.degree; ++i){
		p.a_coeffs[i] = init_coeffs[i];
	}
}

void destroy_poly (poly_t &p){
	delete[] p.a_coeffs;
	p.a_coeffs = nullptr;
}

unsigned int poly_degree(poly_t const &p){
	if(p.a_coeffs == nullptr){
		throw 0;
	}
	else{
		return p.degree;
	}
}

double poly_coeff(poly_t const &p, unsigned int n){
	if(p.a_coeffs == nullptr){
		throw 0;
	}
	else{
		return p.a_coeffs[n];
	}
}

double poly_val(poly_t const &p, double const x){
	double value = p.a_coeffs[p.degree];
	if(p.a_coeffs == nullptr){
		throw 0;
	}
	else{
		for (unsigned int i = p.degree; i>0; --i){
			value = (value*x) + p.a_coeffs[i-1];
		}
		return value;
	}
}

void poly_add(poly_t &p, poly_t const &q){
	if(p.a_coeffs == nullptr || q.a_coeffs == nullptr){
		throw 0;
	}
	else{
		double temp[(p.degree + 1)];
		for(unsigned int i = 0; i<(p.degree+1); ++i){
			temp[i] = p.a_coeffs[i];
		}

        delete[] p.a_coeffs;
	    p.a_coeffs=nullptr;

	    if(p.degree >= q.degree){
	    	p.a_coeffs = new double[p.degree + 1];

	    	for(unsigned int i = 0; i<(p.degree+1); ++i){
	    		if(i <= q.degree){
	    			p.a_coeffs[i] = q.a_coeffs[i] + temp[i];
	    		}else{
	    			p.a_coeffs[i] = temp[i];
	    		}
//	    		std::cout << p.a_coeffs[i] << std::endl;
	    	}
	    }
	    else{
	    	p.a_coeffs = new double[q.degree+1];
	    	p.degree = q.degree;
	    	for(unsigned int i = 0; i<(q.degree+1); ++i){
	    		if(i <= p.degree){
	    			p.a_coeffs[i] = q.a_coeffs[i] + temp[i];
	    		}else{
	    			p.a_coeffs[i] = q.a_coeffs[i];
	    		}
//	    		std::cout << p.a_coeffs[i] << std::endl;
	    	}
	    }
	}
}

void poly_subtract( poly_t &p, poly_t const &q ){
	if(p.a_coeffs == nullptr || q.a_coeffs == nullptr){
		throw 0;
	}
	else{
		double temp[(p.degree + 1)];
		for(unsigned int i = 0; i<(p.degree+1); ++i){
			temp[i] = p.a_coeffs[i];
		}
	    delete[] p.a_coeffs;
		p.a_coeffs=nullptr;

		if(p.degree >= q.degree){
		    p.a_coeffs = new double[p.degree + 1];

		    for(unsigned int i = 0; i<(p.degree+1); ++i){
		    	if(i <= q.degree){
		    	    p.a_coeffs[i] = temp[i] - q.a_coeffs[i];
		    	 }else{
		    	    p.a_coeffs[i] = temp[i];
		    	 }
//		    	 std::cout << p.a_coeffs[i] << std::endl;
		    }
		 }
		 else{
			 p.a_coeffs = new double[q.degree+1];
			 p.degree = q.degree;

			 for(unsigned int i = 0; i<(q.degree+1); ++i){
		    	    if(i <= p.degree){
		    	    	p.a_coeffs[i] = temp[i] - q.a_coeffs[i];
		    	    }else{
		    	    	p.a_coeffs[i] = -q.a_coeffs[i];
		    	    }
//		    	    std::cout << p.a_coeffs[i] << std::endl;
		    	   }
		   	   }
		}
}

void poly_multiply( poly_t &p, poly_t const &q ){
	if(p.a_coeffs == nullptr || q.a_coeffs == nullptr){
		throw 0;
	}
	else{
		double temp[(p.degree + 1)];
		for(unsigned int i = 0; i<(p.degree+1); ++i){
			temp[i] = p.a_coeffs[i];
		}

		delete[] p.a_coeffs;
	    p.a_coeffs=nullptr;
		p.a_coeffs = new double[p.degree + q.degree + 1];

		for(unsigned int i = 0; i<(p.degree + q.degree + 1); ++i){
			p.a_coeffs[i] = 0;
		}

		for(unsigned int k = 0; k<(p.degree+1); ++k){
			for(unsigned int j = 0; j<(q.degree+1); ++j){
				p.a_coeffs[k+j] += temp[k]*q.a_coeffs[j];
			}
		}

		p.degree = p.degree + q.degree;

//		for(unsigned int m = 0; m<(p.degree+q.degree+1); ++m){
//			std::cout << p.a_coeffs[m] << std::endl;
//		}
	}
}

double poly_divide( poly_t &p, double r ){
	double value = 0.0;
	double remainder = 0.0;
	if(p.a_coeffs == nullptr){
		throw 0;
	}
	else{
		double temp[p.degree + 1];
		for(unsigned int i = 0; i<(p.degree+1); ++i){
			temp[i] = p.a_coeffs[i];
		}

		delete[] p.a_coeffs;
		p.a_coeffs = nullptr;
		p.a_coeffs = new double[p.degree];

		for (int i = (p.degree-1); i>=0; --i){
			if(i == (p.degree-1)){
				p.a_coeffs[i] = temp[i+1];
			}
			else{
			value = r*p.a_coeffs[i+1];
			p.a_coeffs[i] = value + temp[i+1];
			}
		}

		value = r*p.a_coeffs[0];
		remainder = value + temp[0];

		p.degree = p.degree - 1;

//		for(unsigned int m = 0; m<(p.degree+1); ++m){
//			std::cout << p.a_coeffs[m] << std::endl;
//		}

//		std::cout << std::endl;
//		std::cout << remainder << std::endl;

		return remainder;
	}
}

void poly_diff( poly_t &p ){
	if(p.a_coeffs == nullptr){
		throw 0;
	}
	else{
		double temp[p.degree + 1];
		for(unsigned int i = 0; i<(p.degree+1); ++i){
			temp[i] = p.a_coeffs[i];
		}

		delete[] p.a_coeffs;
		p.a_coeffs = nullptr;
		p.a_coeffs = new double[p.degree];

		for(unsigned int k = 0; k<(p.degree); ++k){
			double coefficient = (k+1)*temp[k+1];
			p.a_coeffs[k] = coefficient;
		}

		p.degree = p.degree - 1;
	}

//		for(unsigned int m = 0; m<(p.degree+1); ++m){
//			std::cout << p.a_coeffs[m] << std::endl;
//		}
}

double poly_approx_int( poly_t const &p, double a, double b, unsigned int n ){
	double half_h = ((b-a)/n)/2;
	double height = (b-a)/n;
	double sum = 0.0;
	double k = 0.0;
	double increment = height*k;
	if(p.a_coeffs == nullptr){
		throw 0;
	}
	else{
		for (double i = a; i<= b; i+=increment){
			double value = poly_val(p, i);
			if((i==b)||(i==a)){
				sum += value;
			}
			else{
				sum += (2*value);
			}
			++k;
		}
		double area = half_h*sum;
		std::cout << area << std::endl;
		return area;
	}
}

/*double poly_approx_int( poly_t const &p, double a, double b, unsigned int n ){
	double half_h = ((b-a)/n)/2;
	double sum = 0.0;

	if(p.a_coeffs == nullptr){
		throw 0;
	}
	else{
		double value = p.a_coeffs[p.degree];
		double x = b;
		for (unsigned int i = p.degree; i>0; --i){
			value = (value*x) + p.a_coeffs[i-1];
			if(x==b || x==a){
				sum += value;
			} else{
				sum += (2*value);
			}
			--x;
		}
		double area = half_h*sum;
//		std::cout << area << std::endl;
		return area;
	}
} */
