r"""
Tangent spaces

The class :class:`TangentSpace` implements vector spaces tangent to a
differentiable manifold and the class :class:`TangentVector` implements the
corresponding tangent vectors. 
 
AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014): initial version

REFERENCES:

- Chap. 3 of J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer
  (New York) (2013)

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

from sage.symbolic.ring import SR
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from sage.tensor.modules.free_module_tensor import FiniteRankFreeModuleElement

class TangentVector(FiniteRankFreeModuleElement):
    r"""
    Tangent vector to a differentiable manifold at a given point.
    
    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.geometry.manifolds.tangentspace.TangentSpace`.
    It inherits from
    :class:`~sage.tensor.modules.free_module_tensor.FiniteRankFreeModuleElement`
    since :class:`~sage.geometry.manifolds.tangentspace.TangentSpace`
    inherits from
    :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`. 

    INPUT:

    - ``parent`` -- the tangent space to which the vector belongs (must be an
      instance of
      :class:`~sage.geometry.manifolds.tangentspace.TangentSpace`)
    - ``name`` -- (default: ``None``) string; symbol given to the vector
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      the vector; if none is provided, ``name`` will be used

    EXAMPLES:

    Tangent vector on a 2-dimensional manifold::

        sage: M = Manifold(2, 'M')
        sage: c_xy.<x,y> = M.chart()
        sage: p = M.point((2,3), name='p')
        sage: Tp = p.tangent_space()
        sage: v = Tp((-2,1), name='v') ; v
        tangent vector v at point 'p' on 2-dimensional manifold 'M'
        sage: v.display()
        v = -2 d/dx + d/dy
        sage: v.parent()
        tangent space at point 'p' on 2-dimensional manifold 'M'
        sage: v in Tp
        True

    """
    def __init__(self, parent, name=None, latex_name=None):
        FiniteRankFreeModuleElement.__init__(self, parent, name=name, latex_name=latex_name)
        # Extra data (with respect to FiniteRankFreeModuleElement):
        self._point = parent._point

    def _repr_(self):
        r"""
        String representation of the object.
        """
        desc = "tangent vector"
        if self._name:
            desc += " " + str(self._name)
        desc += " at " + str(self._point)
        return desc

#******************************************************************************

class TangentSpace(FiniteRankFreeModule):
    r"""
    Tangent space to a differentiable manifold at a given point.

    Since the tangent space at a given point `p` of a differentiable
    manifold `M` of dimension `n` is a `n`-dimensional vector space over `\RR`
    without any distinguished basis, the class :class:`TangentSpace` inherits
    from
    :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`,
    which implements free modules of finite rank (hence vector spaces of finite
    dimension) without any distinguished basis. 

    This is a Sage *parent* class, the corresponding *element* class being
    :class:`~sage.geometry.manifolds.tangentspace.TangentVector`.

    INPUT:

    - ``point`` -- (instance of
      :class:`~sage.geometry.manifolds.point.ManifoldPoint`) point `p` at which the
      tangent space is defined.

    EXAMPLES:

    Tangent space on a 2-dimensional manifold::

        sage: M = Manifold(2, 'M')
        sage: c_xy.<x,y> = M.chart()
        sage: p = M.point((-1,2), name='p')
        sage: Tp = p.tangent_space() ; Tp
        tangent space at point 'p' on 2-dimensional manifold 'M'

    Tangent spaces actually belong to a dynamically generated subclass of
    :class:`TangentSpace`::

        sage: type(Tp)
        <class 'sage.geometry.manifolds.tangentspace.TangentSpace_with_category'>

    They are free modules of finite rank over Sage's Symbolic Ring (actually
    vector spaces of finite dimension over `\RR`)::

        sage: isinstance(Tp, FiniteRankFreeModule)
        True
        sage: Tp.base_ring()
        Symbolic Ring
        sage: Tp.category()
        Category of vector spaces over Symbolic Ring
        sage: Tp.rank()
        2
        sage: dim(Tp)
        2

    The tangent space is automatically endowed with bases deduced from the
    vector frames around the point::

        sage: Tp.bases()
        [Basis (d/dx,d/dy) on the tangent space at point 'p' on 2-dimensional manifold 'M']
        sage: M.frames()
        [coordinate frame (M, (d/dx,d/dy))]

    At this stage, only one basis has been defined in the tangent space, but
    new bases can be added from vector frames on the manifolds by means of the
    method :meth:`~sage.geometry.manifolds.vectorframe.VectorFrame.at`, for
    instance, from the frame associated with new coordinates on the manifold::

        sage: c_uv.<u,v> = M.chart()
        sage: c_uv.frame().at(p)
        Basis (d/du,d/dv) on the tangent space at point 'p' on 2-dimensional manifold 'M'
        sage: Tp.bases()
        [Basis (d/dx,d/dy) on the tangent space at point 'p' on 2-dimensional manifold 'M',
         Basis (d/du,d/dv) on the tangent space at point 'p' on 2-dimensional manifold 'M']

    All the bases defined on ``Tp`` are on the same footing. Accordingly the
    tangent space is not in the category of modules with a distinguished
    basis::

        sage: Tp in ModulesWithBasis(SR)
        False

    It is simply in the category of modules, as any instance of
    :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`::

        sage: Tp in Modules(SR)
        True

    Since the base ring is a field, it is actually in the category of
    vector spaces::

        sage: Tp in VectorSpaces(SR)
        True

    A typical element::

        sage: v = Tp.an_element() ; v
        tangent vector at point 'p' on 2-dimensional manifold 'M'
        sage: v.display()
        d/dx + 2 d/dy
        sage: v.parent()
        tangent space at point 'p' on 2-dimensional manifold 'M'

    The zero vector::

        sage: Tp.zero()
        tangent vector zero at point 'p' on 2-dimensional manifold 'M'
        sage: Tp.zero().display()
        zero = 0
        sage: Tp.zero().parent()
        tangent space at point 'p' on 2-dimensional manifold 'M'

    """

    Element = TangentVector

    def __init__(self, point):
        manif = point._manifold
        name = "T_" + str(point._name) + " " + str(manif._name)
        latex_name = r"T_{" + str(point._name) + "}\," + str(manif._latex_name)
        self._point = point
        self._manif = manif
        FiniteRankFreeModule.__init__(self, SR, manif._dim, name=name,
                                      latex_name=latex_name,
                                      start_index=manif._sindex)
        # Initialization of bases of the tangent space from existing vector
        # frames around the point:
        for frame in point._subset._top_frames:
            if point in frame._domain:
                basis = self.basis(symbol=frame._symbol,
                                   latex_symbol=frame._latex_symbol)
                # Names of basis vectors set to those of the frame vector
                # fields:
                n = manif._dim
                for i in range(n):
                    basis._vec[i]._name = frame._vec[i]._name
                    basis._vec[i]._latex_name = frame._vec[i]._latex_name
                basis._name = "(" + \
                        ",".join([basis._vec[i]._name for i in range(n)]) + ")"
                basis._latex_name = r"\left(" + \
                     ",".join([basis._vec[i]._latex_name for i in range(n)])+ \
                     r"\right)"
                basis._symbol = basis._name
                basis._latex_symbol = basis._latex_name
                # Names of cobasis linear forms set to those of the coframe
                # 1-forms:
                coframe = frame.coframe()
                cobasis = basis.dual_basis()
                for i in range(n):
                    cobasis._form[i]._name = coframe._form[i]._name
                    cobasis._form[i]._latex_name = coframe._form[i]._latex_name
                cobasis._name = "(" + \
                     ",".join([cobasis._form[i]._name for i in range(n)]) + ")"
                cobasis._latex_name = r"\left(" + \
                  ",".join([cobasis._form[i]._latex_name for i in range(n)])+ \
                  r"\right)"
                point._frame_bases[frame] = basis
        # The basis induced by the default frame of the manifold subset
        # in which the point has been created is declared the default
        # basis of self:
        def_frame = point._subset._def_frame
        if def_frame in point._frame_bases:
            self._def_basis = point._frame_bases[def_frame]
        # Initialization of the changes of bases from the existing changes of
        # frames around the point:
        for frame_pair, automorph in point._subset._frame_changes.iteritems():
            frame1 = frame_pair[0] ; frame2 = frame_pair[1]
            fr1, fr2 = None, None
            for frame in point._frame_bases:
                if frame1 in frame._subframes:
                    fr1 = frame
                    break
            for frame in point._frame_bases:
                if frame2 in frame._subframes:
                    fr2 = frame
                    break
            if fr1 is not None and fr2 is not None:
                basis1 = point._frame_bases[fr1]
                basis2 = point._frame_bases[fr2]
                auto = self.automorphism()
                for frame, comp in automorph._components.iteritems():
                    basis = None
                    if frame is frame1:
                        basis = basis1
                    if frame is frame2:
                        basis = basis2
                    if basis is not None:
                        cauto = auto.add_comp(basis)
                        for ind, val in comp._comp.iteritems():
                            cauto._comp[ind] = val(point)
                self._basis_changes[(basis1, basis2)] = auto

    def _repr_(self):
        r"""
        String representation of the object.
        """

        description = "tangent space at " + str(self._point)
        return description

    def _an_element_(self):
        r"""
        Construct some (unamed) vector in the tangent space
        """
        resu = self.element_class(self)
        if self._def_basis is not None:
            resu.set_comp()[:] = range(1, self._rank+1)
        return resu

    def dimension(self):
        r"""
        Return the vector space dimension of ``self``.
        """
        return self._rank

    dim = dimension
