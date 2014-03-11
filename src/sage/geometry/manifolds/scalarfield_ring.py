r"""
Scalar field ring

The class :class:`ScalarFieldRing` implements the ring of differentiable scalar
fields on some open set of a  differentiable manifolds over `\RR`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014): initial version

EXAMPLES:

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
from sage.categories.commutative_rings import CommutativeRings
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.symbolic.ring import SR
from sage.geometry.manifolds.scalarfield import ScalarField, ZeroScalarField
#from sage.structure.unique_representation import UniqueRepresentation

class ScalarFieldRing(Parent):
    r"""
    Ring of differentiable scalar fields on an open set of a differentiable 
    manifolds over `\RR`

    INPUT:
    
    - ``domain`` -- the manifold open domain `U` on which the scalar fields are 
      defined (must be an instance of class :class:`OpenDomain`)
    
    """

    Element = ScalarField

    def __init__(self, domain):
        Parent.__init__(self, category=CommutativeRings())
        self.domain = domain
        
    #### Methods required for any Parent 
    def _element_constructor_(self, coord_expression=None, chart=None, 
                              name=None, latex_name=None):
        r"""
        Construct a scalarfield
        """
        if coord_expression == 0:
            return ZeroScalarField(self.domain)
        if isinstance(coord_expression, ScalarField):
            if self.domain.is_subdomain(coord_expression.domain):
                resu = self.element_class(self.domain, 
                                coord_expression=coord_expression.expr(chart), 
                                chart=chart, name=name, 
                                latex_name=latex_name)
            else:
                raise TypeError("Cannot coerce this scalar field to " 
                                                            + str(self.domain))
        else:
            resu = self.element_class(self.domain, 
                                      coord_expression=coord_expression, 
                                      chart=chart, name=name, 
                                      latex_name=latex_name)
        return resu

    def _an_element_(self):
        r"""
        Construct some (unamed) element of the module
        """
        resu = self.element_class(self.domain, coord_expression=1)
        return resu
            
            
    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to self exists from other parent
        """
        if other is SR:
            return True
        elif other is ZZ or QQ:
            return True
        else:
            return False

    #### End of methods required for any Parent 

    def _repr_(self):
        r"""
        String representation of the object.
        """
        return "Ring of the scalar fields on the " + str(self.domain)
        
    
