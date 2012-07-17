Defining Custom Weak Forms
--------------------------

The best study material is the library of default weak forms that is located in 
the directories

::

    hermes/hermes2d/include/weakform_library/

and

::

    hermes/hermes2d/src/weakform_library/

of the Hermes library repo. In the former you will find the following 
header files::

    weakforms_elasticity.h  
    weakforms_hcurl.h    
    weakforms_neutronics.h
    weakforms_h1.h          
    weakforms_maxwell.h

The latter contains the corresponding *.cpp files. Let us describe 
the definition of a few of them here.

DefaulMatrixFormVol
~~~~~~~~~~~~~~~~~~~

For starters, let us see how the default volumetric
weak form corresponding to the integral

.. math::

      \int_{area} coeff(x, y) \, u \, v \, \mbox{d}\bfx

is defined. Below is the header. Copy it into your definitions.h file located 
in your example's directory, rename, and adjust to your needs:

.. sourcecode::
    .

      template<typename Scalar>
      class HERMES_API DefaultMatrixFormVol : public MatrixFormVol<Scalar>
      {
      public:
        DefaultMatrixFormVol(int i, int j,
			     Hermes2DFunction<Scalar>* coeff = HERMES_ONE, std::string area = HERMES_ANY,
                             SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

        DefaultMatrixFormVol(int i, int j, 
          Hermes2DFunction<Scalar>* coeff, Hermes::vector<std::string> areas,
          SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

        ~DefaultMatrixFormVol();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
          Geom<double> *e, ExtData<Scalar> *ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
          Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

        virtual MatrixFormVol<Scalar>* clone();

      private:

        Hermes2DFunction<Scalar>* coeff;
        GeomType gt;
      };

.. latexcode::
    .

      template<typename Scalar>
      class HERMES_API DefaultMatrixFormVol : public MatrixFormVol<Scalar>
      {
      public:
        DefaultMatrixFormVol(int i, int j,
			     Hermes2DFunction<Scalar>* coeff = HERMES_ONE, std::string area = HERMES_ANY,
                             SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

        DefaultMatrixFormVol(int i, int j, 
          Hermes2DFunction<Scalar>* coeff, Hermes::vector<std::string> areas,
          SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

        ~DefaultMatrixFormVol();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
          Geom<double> *e, ExtData<Scalar> *ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
          Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

        virtual MatrixFormVol<Scalar>* clone();

      private:

        Hermes2DFunction<Scalar>* coeff;
        GeomType gt;
      };

This form is a descendant of MatrixFormVol, which will be the case also with any custom matrix form that 
you define. The form is rather general and thus it takes a number of parameters. 
These parameters are: 

* i, j (block index in the Jacobian matrix, needed only if the form is to be used in a systems of equations). If your custom form is written for only one equation, then you can just set i = 0 and j = 0 inside the form's body and drop these indices from the argument list. 
* Hermes2DFunction<Scalar>* coeff = HERMES_ONE is used to pass either a constant, cubic spline, or general space-dependent coefficient into the form. HERMES_ONE is a default parameter that represents a constant function whose value is 1.0.
* std::string area = HERMES_ANY is the material marker that the form should be assigned to. The default value, HERMES_ANY means all material markers. If you know that your form is an integral over the whole domain, you can drop this argument from the list.
* SymFlag sym = HERMES_NONSYM is a flag telling Hermes whether the form is symmetric with respect to the Jacobian matrix. For diagonal blocks (i = j) this is the same as symmetry with respect to the basis function 'u' and test function 'v'. However, if one creates a form at position i, j where i != j and uses HERMES_SYM, then Hermes will place the same form also at the position j, i. 
* GeomType gt = HERMES_PLANAR is a flag telling Hermes whether the form should be treated as planar, axisymmetric wrt. the x axis, or axisymmetric wrt. the y-axis. Default argument is HERMES_PLANAR which means planar.

The header contains the constructor::

        DefaultMatrixFormVol(int i, int j,
			     Hermes2DFunction<Scalar>* coeff = HERMES_ONE, std::string area = HERMES_ANY,
                             SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

Destructor::

        ~DefaultMatrixFormVol();

Mandatory method value() that returns value of the form::

    virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
      Geom<double> *e, ExtData<Scalar> *ext) const;

Mandatory method ord() that returns integration order for any pairs of basis function 'u' and text 
function 'v'::

    virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
      Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

Mandatory method clone() that is used if replication of the form is needed::

    virtual MatrixFormVol<Scalar>* clone();

You can also define any number of private variables in your form::

    private:

      Hermes2DFunction<Scalar>* coeff;
      GeomType gt;

When you define your custom weak form, the only work that you really need to do is to 
define the methods value() and ord().

The following shows the source of the DefaultMatrixFormVol. It is a bit lengthy because it 
covers the three geometrical cases (planar, axisymmetric in x and axisymmetric in y) as
well as the real and complex versions, as well as versions with one and multiple material 
markers. This is 12 cases in total. Up to this, its definition is straightforward.

Constructor (real case, single material marker)::

      template<>
      DefaultMatrixFormVol<double>::DefaultMatrixFormVol
      (int i, int j, Hermes2DFunction<double>* coeff, std::string area, SymFlag sym, GeomType gt)
        : MatrixFormVol<double>(i, j, area, sym), coeff(coeff), gt(gt)
      {
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if (coeff == HERMES_ONE)
          this->coeff = new Hermes2DFunction<double>(1.0);
      }

Constructor (complex case, single material marker)::

      template<>
      DefaultMatrixFormVol<std::complex<double> >::DefaultMatrixFormVol
        (int i, int j, Hermes2DFunction<std::complex<double> >* coeff, std::string area, SymFlag sym, GeomType gt)
        : MatrixFormVol<std::complex<double> >(i, j, area, sym), coeff(coeff), gt(gt)
      {
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if (coeff == HERMES_ONE)
          this->coeff = new Hermes2DFunction<std::complex<double> >(std::complex<double>(1.0, 0.0));
      }

Constructor (real case, multiple material markers)::

      template<>
      DefaultMatrixFormVol<double>::DefaultMatrixFormVol
        (int i, int j, 
        Hermes2DFunction<double>* coeff, Hermes::vector<std::string> areas, SymFlag sym, GeomType gt)
        : MatrixFormVol<double>(i, j, areas, sym), coeff(coeff), gt(gt)
      {
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if (coeff == HERMES_ONE)
          this->coeff = new Hermes2DFunction<double>(1.0);
      }

Constructor (complex case, multiple material markers)::

      template<>
      DefaultMatrixFormVol<std::complex<double> >::DefaultMatrixFormVol
        (int i, int j, 
        Hermes2DFunction<std::complex<double> >* coeff, Hermes::vector<std::string> areas, SymFlag sym, GeomType gt)
        : MatrixFormVol<std::complex<double> >(i, j, areas, sym), coeff(coeff), gt(gt)
      {
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if (coeff == HERMES_ONE)
          this->coeff = new Hermes2DFunction<std::complex<double> >(std::complex<double>(1.0, 0.0));
      }

Destructor::

      template<typename Scalar>
      DefaultMatrixFormVol<Scalar>::~DefaultMatrixFormVol()
      {
        // FIXME: Should be deleted here only if it was created here.
        //if (coeff != HERMES_ONE) delete coeff;
      };

Value (this is where your integral is really defined)::

      template<typename Scalar>
      Scalar DefaultMatrixFormVol<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext) const
      {
        Scalar result = 0;
        if (gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
          }
        }
        else {
          if (gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
            }
          }
        }

        return result;
      }

The arguments have the following meaning:

* n: number of integration points.
* wt[]: array of integration weights.
* u_ext[]: array of pointers to previous Newton iterations.
* u: pointer to basis function.
* v: pointer to test function.
* e: pointer to geometrical data such as physical coordinates of integration points, normal and tangential vectors, etc. 
* ext: pointer to user-defined external functions (must be registered with the weak form before they can be used).

In the form, one can access the value of the basis function at integration point i via u->val[i],
its x-derivative via u->dx[i] and its y->derivative via u->y[i]. Same holds for the test function 'v'.
Via e->x[i] and e->y[i] one can access the physics coordinates of integration point i.


Order::

      template<typename Scalar>
      Ord DefaultMatrixFormVol<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
      {
        Ord result = Ord(0);
        if (gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
          }
        }
        else {
          if (gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
            }
          }
        }

        return result;
      }

It is very important to define this method carefully since underestimating 
the integration order may cause unwanted numerical errors, and over-estimating it 
makes the computation time longer. Higher-order numerical quadrature is in general 
quite expensive. In most cases it is enough to copy and paste
the same code that is used to define the value of the form. Hermes has a parser that 
will determine the order automatically. 

However, one needs to be careful when 
using non-polynomial functions such as sin(), exp() log() as Hermes will automatically
set the highest available integration order. In such situations, the user needs to 
say what integration order he/she wants to use. Should it be a fixed order 10? Or should 
it depend on the polynomial degrees of the basis and test functions? Should 
it be their maximum or their sum? And how accurate should the FEM solution be anyway?
This really depends on the concrete application and Hermes cannot make this decision. 

Clone::

      template<typename Scalar>
      MatrixFormVol<Scalar>* DefaultMatrixFormVol<Scalar>::clone()
      {
        return new DefaultMatrixFormVol<Scalar>(*this);
      }

This method is a formality but it needs to be there.
