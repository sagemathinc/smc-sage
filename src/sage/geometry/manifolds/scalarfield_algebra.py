r"""
Scalar field algebra

The class :class:`ScalarFieldAlgebra` implements the commutative algebra 
`C^\infty(U)` of differentiable scalar fields on some open subset `U` of a 
differentiable manifold `M` over `\RR`. By *scalar field*, it is meant a 
function `U\rightarrow \RR`. 
`C^\infty(U)` is an algebra over `\RR`, whose ring product is the pointwise 
multiplication of real-valued functions, which is clearly commutative. 
In the present implementation, the field `\RR`, over which the 
albegra `C^\infty(U)` is constructed, is represented by Sage's symbolic ring
SR. 

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014): initial version

"""

#******************************************************************************
#       Copyright (C) 2014 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2014 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.symbolic.ring import SR
from sage.geometry.manifolds.scalarfield import ScalarField, ZeroScalarField

class ScalarFieldAlgebra(UniqueRepresentation, Parent):
    r"""
    Commutative algebra `C^\infty(U)` of differentiable functions 
    `U \rightarrow \RR`, where `U` is an open subset of some differentiable 
    manifold over `\RR`. `C^\infty(U)` is a commutative algebra over `\RR`,
    the latter being represented by Sage's symbolic ring SR. 
    
    The class :class:`ScalarFieldAlgebra` inherits from 
    :class:`~sage.structure.parent.Parent`, with the category set to
    :class:`~sage.categories.commutative_algebras.CommutativeAlgebras`. 
    The corresponding Element class is 
    :class:`~sage.geometry.manifolds.scalarfield.ScalarField`. 

    INPUT:
    
    - ``domain`` -- the manifold open subset `U` on which the scalar fields are 
      defined (must be an instance of class 
      :class:`~sage.geometry.manifolds.domain.OpenDomain`)
    
    EXAMPLES:
    
    Ring of scalar fields on some open subset of a 2-dimensional manifold::
    
        sage: M = Manifold(2, 'M')
        sage: U = M.open_domain('U')
        sage: X.<x,y> = U.chart()
        sage: CU = U.scalar_field_algebra() ; CU
        algebra of scalar fields on the open domain 'U' on the 2-dimensional manifold 'M'
        sage: latex(CU)
        C^\infty(U)

    `C^\infty(U)` belongs to the category of commutative algebras over `\RR` 
    (represented here by Sage's symbolic ring)::
    
        sage: CU.category()
        Category of commutative algebras over Symbolic Ring

    The elements of `C^\infty(U)` are scalar fields::
    
        sage: CU.an_element()
        scalar field on the open domain 'U' on the 2-dimensional manifold 'M'
        sage: CU.an_element().view()
        (x, y) |--> 2
        sage: f = U.scalar_field(x + y^2, name='f')
        sage: f.parent()
        algebra of scalar fields on the open domain 'U' on the 2-dimensional manifold 'M'
        sage: f in CU
        True

    The zero element::
    
        sage: CU.zero()
        zero scalar field on the open domain 'U' on the 2-dimensional manifold 'M'
        sage: CU.zero().view()
        (x, y) |--> 0
        sage: CU.zero() == CU(0)
        True

    The unit element::
    
        sage: CU.one()
        scalar field on the open domain 'U' on the 2-dimensional manifold 'M'
        sage: CU.one().view()
        (x, y) |--> 1        
        sage: CU.one() == CU(1)
        True

    As for any parent in Sage, elements can be constructed by means of the 
    __call__ operator, by providing some coordinate expression::
    
        sage: g = CU(sin(x)) ; g  # coordinate expression in the default chart
        scalar field on the open domain 'U' on the 2-dimensional manifold 'M'
        sage: g.view()
        (x, y) |--> sin(x)
        sage: Y.<u,v> = U.chart()
        sage: g = CU({Y: exp(u+v)}) ; g  # coordinate expression in a non-default chart
        scalar field on the open domain 'U' on the 2-dimensional manifold 'M'
        sage: g.view(Y)
        (u, v) |--> e^(u + v)
     
    The algebra `C^\infty(U)` coerces to `C^\infty(V)` where `V` is any open 
    subset of `U`::
    
        sage: V = U.open_domain('V')
        sage: XV = X.restrict(V, x^2+y^2<1)
        sage: CV = V.scalar_field_algebra()
        sage: CV.has_coerce_map_from(CU)
        True

    The coercion map is nothing but the restriction to `V` of scalar fields
    on `U`::
    
        sage: f_V = CV(f) ; f_V
        scalar field on the open domain 'V' on the 2-dimensional manifold 'M'
        sage: f_V.view()
        (x, y) |--> y^2 + x

    The coercion map allows for the addition of elements of `C^\infty(U)` 
    with elements of `C^\infty(V)`, the result being an element of 
    `C^\infty(V)`::
    
        sage: h = V.scalar_field(x*y)
        sage: f.parent() , h.parent()
        (algebra of scalar fields on the open domain 'U' on the 2-dimensional manifold 'M',
         algebra of scalar fields on the open domain 'V' on the 2-dimensional manifold 'M')
        sage: s = f + h
        sage: s.parent()
        algebra of scalar fields on the open domain 'V' on the 2-dimensional manifold 'M'
        sage: s.view()
        (x, y) |--> x*y + y^2 + x
        sage: s == f_V + h
        True
        sage: s == h + f
        True

    """

    Element = ScalarField

    def __init__(self, domain):
        Parent.__init__(self, base=SR, category=CommutativeAlgebras(SR))
        self.domain = domain
        self._populate_coercion_lists_()
        
    #### Methods required for any Parent 
    def _element_constructor_(self, coord_expression=None, name=None, 
                              latex_name=None):
        r"""
        Construct a scalarfield
        """
        if coord_expression == 0:
            return ZeroScalarField(self.domain)
        if isinstance(coord_expression, ScalarField):
            if self.domain.is_subdomain(coord_expression.domain):
                # restriction of the scalar field to self.domain:
                sexpress = {}
                for chart, funct in coord_expression.express.items():
                    for schart in self.domain._atlas:
                        if schart in chart.subcharts:
                            sexpress[schart] = funct.expr()
                resu = self.element_class(self.domain, 
                                          coord_expression=sexpress, name=name, 
                                          latex_name=latex_name)
            else:
                raise TypeError("Cannot coerce the " + str(coord_expression) +
                                "to a scalar field on the " + str(self.domain))
        else:
            resu = self.element_class(self.domain, 
                                      coord_expression=coord_expression, 
                                      name=name, latex_name=latex_name)
        return resu

    def _an_element_(self):
        r"""
        Construct some (unamed) element of the module
        """
        return self.element_class(self.domain, coord_expression=2)
            
            
    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to self exists from other parent
        """
        if other is SR:
            return True  # coercion from the base ring (multiplication by the 
                         # algebra unit, i.e. self.one())
        elif other is ZZ:
           return True   # important to define self(1) (for self.one())
        elif other is QQ:
            return True
        elif isinstance(other, ScalarFieldAlgebra):
            return self.domain.is_subdomain(other.domain)
        else:
            return False

    #### End of methods required for any Parent 

    def _repr_(self):
        r"""
        String representation of the object.
        """
        return "algebra of scalar fields on the " + str(self.domain)
        
    def _latex_(self):
        r"""
        LaTeX representation of the object.
        """
        return r"C^\infty("  + self.domain.latex_name + ")"
