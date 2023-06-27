#ifndef __function__
#define __function__

#include <cmath>

class function {

public:
  virtual double Eval(double) const = 0;			//pure virtual method for evaluating the function
  virtual ~function() { ; };
};

class fun : public function {

public:
  fun() {}
  ~fun() {}
  double Eval(double x) const override { return M_PI / 2 * cos(M_PI / 2 * x); }
};

#endif // __function__