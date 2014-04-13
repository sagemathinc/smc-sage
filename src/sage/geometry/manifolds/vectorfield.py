r"""
Vector fields

The class :class:`VectorField` implements vector fields on differentiable 
manifolds over `\RR`. 


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013, 2014) : initial version

"""

#******************************************************************************
#       Copyright (C) 2013, 2014 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2013, 2014 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.free_module_tensor import FiniteFreeModuleElement
from tensorfield import TensorFieldParal

class VectorFieldParal(FiniteFreeModuleElement, TensorFieldParal):
    r"""
    Vector field on an open set of a differentiable manifold, 
    with values on parallelizable open subset of a differentiable manifold. 
    
    An instance of this class is a vector field along an open subset `U` 
    of some manifold `S` with values in a parallelizable open subset `V` 
    of a manifold `M`, via a differentiable mapping `\Phi: U \rightarrow V`. 
    The standard case of a vector field *on* a manifold corresponds to `S=M`, 
    `U=V` and `\Phi = \mathrm{Id}`. Another common case is `\Phi` being an
    immersion.

    INPUT:
    
    - ``vector_field_module`` -- free module `\mathcal{X}(U,\Phi)` of vector 
      fields along `U` with values on `\Phi(U)\subset V \subset M`
    - ``name`` -- (default: None) name given to the vector field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the vector field; 
      if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A vector field on a 3-dimensional manifold::
    
        sage: M = Manifold(3, 'M')
        sage: c_xyz.<x,y,z> = M.chart('x y z')
        sage: v = M.vector_field('V') ; v
        vector field 'V' on the 3-dimensional manifold 'M'
        sage: latex(v)
        V
    
    Vector fields are considered as elements of a module over the ring of
    scalar fields on `M`::
    
        sage: v.parent()
        free module X(M) of vector fields on the 3-dimensional manifold 'M'
        sage: v.parent().base_ring()
        ring of scalar fields on the 3-dimensional manifold 'M'
        sage: v.parent() is M.vector_field_module()
        True

    A vector field is a tensor field of rank 1 and of type (1,0)::
    
        sage: v.tensor_rank
        1
        sage: v.tensor_type
        (1, 0)

    Components of a vector field with respect to a given frame::
    
        sage: e = M.vector_frame('e') ; M.set_default_frame(e)
        sage: v[0], v[1], v[2] = (1, 4, 9)  # components on M's default frame (e)
        sage: v.comp()
        1-index components w.r.t. vector frame (M, (e_0,e_1,e_2))
    
    The totality of the components are accessed via the operator [:]::
    
        sage: v[:] = (1, 4, 9) # equivalent to v[0], v[1], v[2] = (1, 4, 9)
        sage: v[:]
        [1, 4, 9]
        
    The components are also read on the expansion on the frame 'e', as provided
    by the method :meth:`view`::
    
        sage: v.view()   # displays the expansion on the manifold's default frame (e)
        V = e_0 + 4 e_1 + 9 e_2
    
    A subset of the components can be accessed by means of Python's slice 
    notation::
        
        sage: v[1:] = (-2, -3)
        sage: v[:]
        [1, -2, -3]
        sage: v[:2]
        [1, -2]
        
    The components are instances of the class 
    :class:`~sage.tensor.modules.comp.Components`::
    
        sage: type(v.comp())
        <class 'sage.tensor.modules.comp.Components'>

    Components in another frame::
    
        sage: f = M.vector_frame('f')
        sage: for i in range(3):
        ...       v.set_comp(f)[i] = (i+1)**3
        ...
        sage: v.comp(f)[2]
        27
        sage: v.view(f)
        V = f_0 + 8 f_1 + 27 f_2

    The range of the indices depends on the convention set for the manifold::
        
        sage: M = Manifold(3, 'M', start_index=1)
        sage: c_xyz.<x,y,z> = M.chart('x y z')
        sage: e = M.vector_frame('e') ; M.set_default_frame(e)
        sage: v = M.vector_field('V')
        sage: (v[1], v[2], v[3]) = (1, 4, 9)
        sage: v[0]
        Traceback (most recent call last):
        ...
        IndexError: Index out of range: 0 not in [1,3]

    A vector field acts on scalar fields (derivation along the vector field)::
    
        sage: M = Manifold(2, 'M')            
        sage: c_cart.<x,y> = M.chart('x y')
        sage: f = M.scalar_field(x*y^2, name='f')  
        sage: v = M.vector_field('v')         
        sage: v[:] = (-y, x)
        sage: v.view()
        v = -y d/dx + x d/dy
        sage: v(f)
        scalar field 'v(f)' on the 2-dimensional manifold 'M'
        sage: v(f).expr()
        2*x^2*y - y^3
        sage: latex(v(f))
        v\left(f\right)

    """
    def __init__(self, vector_field_module, name=None, latex_name=None):
        FiniteFreeModuleElement.__init__(self, vector_field_module, name=name, 
                                         latex_name=latex_name)
        # TensorFieldParal attributes:
        self.domain = vector_field_module.domain
        self.ambient_domain = vector_field_module.ambient_domain
        # Initialization of derived quantities:
        TensorFieldParal._init_derived(self) 
        # Initialization of list of quantities depending on self:
        self._init_dependencies()
        
    def _repr_(self) :
        r"""
        String representation of the object.
        """
        description = "vector field "
        if self.name is not None:
            description += "'%s' " % self.name
        return self._final_repr(description)

    def _new_instance(self):
        r"""
        Create a :class:`VectorFieldParal` instance. 
        
        """
        return VectorFieldParal(self.fmodule)

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        TensorFieldParal._del_derived(self)
        self._del_dependencies()
        
    def _init_dependencies(self):
        r"""
        Initialize list of quantities that depend on ``self``
        """
        self._lie_der_along_self = {}

    def _del_dependencies(self):
        r"""
        Clear list of quantities that depend on ``self``
        """
        if self._lie_der_along_self != {}:
            for idtens, tens in self._lie_der_along_self.items():
                del tens._lie_derivatives[id(self)]
            self._lie_der_along_self.clear()

    def __call__(self, scalar):
        r"""
        Action on a scalar field.
            
        INPUT:
            
        - ``scalar`` -- scalar field `f`
            
        OUTPUT:
            
        - scalar field representing the derivative of `f` along the vector 
          field, i.e. `v^i \frac{\partial f}{\partial x^i}`
          
        EXAMPLES:
        
        Action of a vector field on a scalar field on a 2-dimensional manifold::
        
            sage: M = Manifold(2, 'M')            
            sage: c_cart.<x,y> = M.chart('x y')
            sage: f = M.scalar_field(x*y^2)  
            sage: v = M.vector_field()         
            sage: v[:] = (-y, x)
            sage: v(f)
            scalar field on the 2-dimensional manifold 'M'
            sage: v(f).expr()
            2*x^2*y - y^3
          
        """
        from diffform import OneFormParal
        from scalarfield import ScalarField, ZeroScalarField
        if isinstance(scalar, OneFormParal):  #!# it should be OneForm
            # This is actually the action of the vector field on a 1-form, 
            # as a tensor field of type (1,0):
            return scalar(self)
        if not isinstance(scalar, ScalarField):
            raise TypeError("The argument must be a scalar field")
        if not scalar.domain.is_subdomain(self.domain):
            raise ValueError("The scalar field and the vector are defined " +
                             "on different domains.")
        if isinstance(scalar, ZeroScalarField):
            return scalar
        # search for a commont chart: 
        chart = None
        def_chart = self.domain.def_chart
        if def_chart in scalar.express:
            if def_chart.frame in self.components:
                chart = def_chart
        else:
            for kchart in scalar.express:
                if kchart.frame in self.components: 
                    chart = kchart
                    break
        if chart is None:
            raise ValueError("No common chart found.")
        v = self.comp(chart.frame)
        f = scalar.function_chart(chart) 
        res = 0 
        for i in scalar.manifold.irange():
            res += v[i, chart] * f.diff(i)
        # Name of the output:
        res_name = None
        if self.name is not None and scalar.name is not None:
            res_name = self.name + "(" + scalar.name + ")"
        # LaTeX symbol for the output:
        res_latex = None
        if self.latex_name is not None and scalar.latex_name is not None:
            res_latex = self.latex_name + r"\left(" + scalar.latex_name + \
                        r"\right)"
        return res.scalar_field(name=res_name, latex_name=res_latex)
