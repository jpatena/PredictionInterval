// Jaron Patena
// Estimate LOC using prediction interval
// Compile: c++ -o LOC LOCestimate.cpp
// Execute: ./LOC
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <boost/math/distributions/students_t.hpp>

using namespace std;

class StdDeviation
{
private:
	int max;
	double value[100];
	double mean;

public:
	double CalcMean();
	double CalcVariance();
	double CalcSampleVariance();
	int SetValues(double *p, int count);
	double CalcStdDeviation();
	double CalcSampleStdDeviation();
};


double StdDeviation::CalcMean()
{
	double sum = 0.0;
	for(int i = 0; i < max; i++)
		sum += value[i];
	return (sum / max);
}

double StdDeviation::CalcVariance()
{
	mean = CalcMean();
	double temp = 0;

	for(int i = 0; i < max; i++)
		temp += (value[i] - mean) * (value[i] - mean) ;
	return temp / max;
}

double StdDeviation::CalcSampleVariance()
{
	mean = CalcMean();
	double temp = 0;

	for(int i = 0; i < max; i++)
		temp += (value[i] - mean) * (value[i] - mean) ;
	return temp / (max - 1);
}

int StdDeviation::SetValues(double *p, int count)
{
  if(count > 100) return -1;
	max = count;

	for(int i = 0; i < count; i++)
		value[i] = p[i];
	return 0;
}

double StdDeviation::CalcStdDeviation()
{
  return sqrt(CalcVariance());
}

double StdDeviation::CalcSampleStdDeviation()
{
	return sqrt(CalcSampleVariance());
}

class StatsCalc
{
private:
	double XSeries[100];
	double YSeries[100];
	double xmean;
	double ymean;
	int max;
	double p;
	double BOne;
	double BZero;
	double conf;
	double t;
	double xk;
	double yk;
	double range;
	double UPI;
	double LPI;
	StdDeviation x;
	StdDeviation y;
public:
	void SetValues(double *xvalues, double *yvalues, int count);
	double CalcCovariance();
	double CalcCorrelation();
	double CalcBone();
	double CalcBzero();
	double LeastSquares(double x);
	double CalcVar();
	double CalcStdDev();
	double CalctDistribution(double confidence);
	double CalcRange();
	double CalcUPI();
	double CalcLPI();
	void DispResults();
};

void StatsCalc::SetValues(double *xvalues, double *yvalues, int count)
{
	for(int i = 0; i < count; i++)
	{
		XSeries[i] = xvalues[i];
		YSeries[i] = yvalues[i];
	}
	x.SetValues(xvalues, count);
	y.SetValues(yvalues, count);
	max = count; // n
}

double StatsCalc::CalcCovariance()
{
	xmean = x.CalcMean();
	ymean = y.CalcMean();

	double total = 0;
	for(int i = 0; i < max; i++)
		total += (XSeries[i] - xmean) * (YSeries[i] - ymean);
	return total / max;
}

double StatsCalc::CalcCorrelation() // r
{
	double cov = CalcCovariance();
	double correlation = cov / (x.CalcStdDeviation() * y.CalcStdDeviation());
	
	return correlation;
}

// line regression
double StatsCalc::CalcBone() // b1
{
	return CalcCorrelation()*(y.CalcStdDeviation() / x.CalcStdDeviation());
}

double StatsCalc::CalcBzero() // b0
{
	return (ymean-(CalcBone()*xmean));
}

// prediction interval
double StatsCalc::CalcVar()
{
	double step1, step2;
	int i=0;
	// 1/n-2
	step1 = 1.0/(max-2);
	// sigma notation(yi-b0-b1xi)^2
	for(; i < 8; i++)
	{
		step2 += pow(YSeries[i]-CalcBzero()-(CalcBone()*XSeries[i]),2);
	}
	return step1*step2;
}

double StatsCalc::CalcStdDev()
{
	return sqrt(CalcVar());
}

double StatsCalc::LeastSquares(double x) // y = b0 + b1x (prediction)
{
	xk = x; // estimated object LOC
	yk = CalcBzero()+(CalcBone()*xk);
	return yk;
}

double StatsCalc::CalctDistribution(double confidence) // t(95%) or t(85%)
{
	conf = confidence;
	double alpha = 1.0-conf;
	boost::math::students_t dist(max-2.0);
	t = boost::math::quantile(complement(dist, alpha/2.0));
	return t;
}

double StatsCalc::CalcRange()
{
	double step1, step2, step3, step4, step5, step6;
	// tdistribution*std dev
	step1 = t*CalcStdDev();
	
	// (xk-xmean)^2
	step2 = pow((xk-xmean),2.0);

	//	sigma notation(xi-xmean)^2
	for(int i=0; i < max; i++)
	{
		step3 += pow(XSeries[i]-xmean,2.0);
	}
	// (xk-xmean)^2 / sigma notation(xi-xmean)^2
	step4 = step2/step3;

	// sqrt[1 + 1/n +  (xk-xmean)^2 / sigma notation(xi-xmean)^2]
	step5 = 1.0+(1.0/max)+step4;
	step6 = sqrt(step5);
	range = step1*step6;

	UPI = yk + range;
	LPI = yk - range;

	return range;
}

double StatsCalc::CalcUPI()
{
	return UPI;
}
double StatsCalc::CalcLPI()
{
	return LPI;
}

void StatsCalc::DispResults()
{
	CalcCorrelation();
	cout << "b0: " << CalcBzero() << endl;
	cout << "b1: " << CalcBone() << endl;
	cout << "prediction (y): " << LeastSquares(400) << endl;
	cout << "std dev: " << CalcStdDev() << endl << endl;

	cout << "t(85%): " << CalctDistribution(0.85) << endl;
	cout << "Range(85%): " << CalcRange() << endl;
	cout << "UPI: " << CalcUPI() << " LPI: " << CalcLPI() << endl << endl;

	cout << "t(95%): " << CalctDistribution(0.95) << endl;
	cout << "Range(95%): " << CalcRange() << endl;
	cout << "UPI: " << CalcUPI() << " LPI: " << CalcLPI() << endl;
}

int main() // main program
{
	StatsCalc num;
	// Estimated object LOC
	double xarr[] = { 130, 650, 99, 150, 128, 302, 95, 945};
	// Actual New and Changed LOC
	double yarr[] = { 186, 699, 132, 272, 291, 331, 199, 1890};
	
	num.SetValues(xarr,yarr,sizeof(xarr) / sizeof(xarr[0]));
	num.DispResults();

	return 0;
}
