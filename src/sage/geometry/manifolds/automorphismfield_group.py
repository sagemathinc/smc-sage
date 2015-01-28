r"""
Group of tangent-space automorphism fields


AUTHORS:

- Eric Gourgoulhon (2015): initial version

"""

#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

#from sage.structure.unique_representation import UniqueRepresentation
#from sage.structure.parent import Parent
from sage.tensor.modules.free_module_linear_group import FreeModuleLinearGroup
from sage.geometry.manifolds.vectorfield_module import VectorFieldFreeModule
from sage.geometry.manifolds.automorphismfield import AutomorphismFieldParal

class AutomorphismFieldParalGroup(FreeModuleLinearGroup):
    r"""
    General linear group of the module of vector fields along an open subset
    `U` of some manifold `S` with values in a parallelizable open subset `V`
    of a manifold `M`.

    Given an open subset `U` of a manifold `S` and a differentiable mapping
    `\Phi: U \rightarrow V = \Phi(U) \subset M`, where `M` is a manifold
    and `V` is parallelizable, this class implements the general linear
    group `\mathrm{GL}(\mathcal{X}(U,\Phi))` of the free module
    `\mathcal{X}(U,\Phi)` of vector fields along `U` with values in
    `V=\Phi(U)`. Since `V` is parallelizable, `\mathcal{X}(U,\Phi)` is indeed
    a free module over `C^\infty(U)`, the ring (algebra) of smooth scalar
    fields on `U`. Elements of this group are fields along `U` of automorphisms
    of the tangent spaces of `V=\Phi(U)`. 

    If `V=\Phi(U)` is not parallelizable, the class
    :class:`AutomorphismFieldGroup` must be used instead.

    This is a Sage *parent* class, the corresponding *element* class being
    :class:`~sage.geometry.manifolds.automorphismfield.AutomorphismFieldParal`.

    INPUT:

    - ``vector_field_module`` -- free module `\mathcal{X}(U,\Phi)` of vector
      fields along `U` with values on `V`


    EXAMPLES:

    Group of tangent-space automorphism fields of a 2-dimensional
    parallelizable manifold::

        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: XM = M.vector_field_module() ; XM
        free module X(M) of vector fields on the 2-dimensional manifold 'M'
        sage: from sage.geometry.manifolds.automorphismfield_group import \
        ....:                                       AutomorphismFieldParalGroup
        sage: G = AutomorphismFieldParalGroup(XM)  ; G
        General linear group of the free module X(M) of vector fields on the
         2-dimensional manifold 'M'
        sage: latex(G)
        \mathrm{GL}\left( \mathcal{X}\left(M\right) \right)
        sage: type(G)
        <class 'sage.geometry.manifolds.automorphismfield_group.AutomorphismFieldParalGroup_with_category'>
    
    Instead of importing ``AutomorphismFieldParalGroup`` in the global name
    space, it is recommended to use the method
    :meth:`~sage.geometry.manifolds.domain.ManifoldOpenSubset.automorphism_field_group`::
    
        sage: G = M.automorphism_field_group() ; G
        General linear group of the free module X(M) of vector fields on the
         2-dimensional manifold 'M'

    There is a unique instance of this group::
    
        sage: G is M.automorphism_field_group()
        True

    ``G`` is nothing but the general linear group of the module
    `\mathcal{X}(M)`::

        sage: G is XM.general_linear_group()
        True

    ``G`` is a group::

        sage: G.category()
        Category of groups
        sage: G in Groups()
        True

    It is not an abelian group::

        sage: G in CommutativeAdditiveGroups()
        False

    ``G`` is a *parent* object, whose elements are tangent-space automorphisms::

        sage: G.Element
        <class 'sage.geometry.manifolds.automorphismfield.AutomorphismFieldParal'>
        sage: a = G.an_element() ; a
        field of tangent-space automorphisms on the 2-dimensional manifold 'M'
        sage: a.parent() is G
        True

    As automorphisms of `\mathcal{X}(M)`, the elements of ``G`` map a vector
    field to a vector field::
    
        sage: v = XM.an_element() ; v
        vector field on the 2-dimensional manifold 'M'
        sage: v.display()
        2 d/dx + 2 d/dy
        sage: a(v)
        vector field on the 2-dimensional manifold 'M'
        sage: a(v).display()
        2 d/dx - 2 d/dy

    Indeed the matrix of ``a`` with respect to the frame
    `(\partial_x,\partial_y)` is::

        sage: a[X.frame(),:]
        [ 1  0]
        [ 0 -1]    

    The elements of ``G`` can also be considered as tensor fields of
    type (1,1)::

        sage: a.tensor_type()
        (1, 1)
        sage: a.tensor_rank()
        2
        sage: a.domain()
        2-dimensional manifold 'M'
        sage: a.display()
        d/dx*dx - d/dy*dy
        
    The identity element of the group ``G`` is::

        sage: id = G.one() ; id
        field of tangent-space identity maps on the 2-dimensional manifold 'M'
        sage: id*a == a
        True
        sage: a*id == a
        True
        sage: a*a^(-1) == id
        True
        sage: a^(-1)*a == id
        True
    
    """

    Element = AutomorphismFieldParal

    def __init__(self, vector_field_module):
        r"""
        See :class:`AutomorphismfieldParalGroup` for documentation and
        examples.

        TEST::

            sage: M = Manifold(2, 'M') ; M
            2-dimensional manifold 'M'
            sage: X.<x,y> = M.chart()
            sage: from sage.geometry.manifolds.automorphismfield_group import \
            ....:                                   AutomorphismFieldParalGroup
            sage: G = AutomorphismFieldParalGroup(M.vector_field_module()) ; G
            General linear group of the free module X(M) of vector fields on
             the 2-dimensional manifold 'M'
            sage: TestSuite(G).run()
        
        """
        if not isinstance(vector_field_module, VectorFieldFreeModule):
            raise TypeError("{} is not a free module of vector fields".format(
                            vector_field_module))
        FreeModuleLinearGroup.__init__(self, vector_field_module)
