r"""
Differential forms

Let `S` and `M` be two differentiable manifolds over `\RR`.
Given a positive integer `p`, an open subset `U` of `S`,  an open subset
`V` of `M` and differentiable mapping `\Phi: U \rightarrow V \subset M`,
a *differential form of degree* `p`, or *p-form*,
*along* `U` *with values in* `V` is a field along `U` of alternating
multilinear forms of degree `p` in the tangent spaces to `V`.
The standard case of a differential form *on* a manifold corresponds to `S=M`,
`U=V` and `\Phi = \mathrm{Id}_U`. Other common cases are `\Phi` being an
immersion and `\Phi` being a curve in `V` (`U` is then an open interval
of `\RR`).

Two classes implement differential forms, depending whether the open
set `V` is parallelizable:

* :class:`DiffFormParal` for the case where `V` is parallelizable
* :class:`DiffForm` for the generic case, i.e. `V` not assumed
  parallelizable.

.. NOTE::

    A difference with the preceding Sage class
    :class:`~sage.tensor.differential_form_element.DifferentialForm`
    is that the present class lies at the tensor field level. Accordingly, an
    instance of :class:`DiffForm` can have various sets of components, each in
    a different coordinate system or coframe, while the class
    :class:`~sage.tensor.differential_form_element.DifferentialForm` considers
    differential forms at the component level in a fixed chart. In this
    respect, the class
    :class:`~sage.tensor.differential_form_element.DifferentialForm` is closer
    to the class :class:`~sage.tensor.modules.comp.CompFullyAntiSym` than
    to :class:`DiffForm`

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013, 2014): initial version
- Joris Vankerschaver (2010): developed a previous class,
  :class:`~sage.tensor.differential_form_element.DifferentialForm` (cf. the
  above note), which inspired the storage of the non-zero components as a
  dictionary whose keys are the indices.

REFERENCES:

- S. Kobayashi & K. Nomizu : *Foundations of Differential Geometry*, vol. 1,
  Interscience Publishers (New York) (1963)
- J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer (New York)
  (2013)

"""

#******************************************************************************
#       Copyright (C) 2013, 2014 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2013, 2014 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2010 Joris Vankerschaver <joris.vankerschaver@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.free_module_alt_form import FreeModuleAltForm
from tensorfield import TensorField, TensorFieldParal

class DiffForm(TensorField):
    r"""
    Differential form with values in an open subset of a differentiable
    manifold.

    Given an positive integer `p`, an open set `U` of a manifold `S`,
    an open set `V` of a manifold `M`  and a differentiable mapping
    `\Phi: U \rightarrow V`, an instance of this class is a field along `U`
    of alternating multilinear forms of degree `p` in the tangent spaces to
    `V`. The standard case of a differential form *on* a manifold corresponds
    to `S=M`, `U=V` and `\Phi = \mathrm{Id}_U`. Other common cases are `\Phi`
    being an immersion and `\Phi` being a curve in `V` (`U` is then an open
    interval of `\RR`).

    If `V` is parallelizable, the class :class:`DiffFormParal` must be
    used instead.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.geometry.manifolds.diffform_module.DiffFormModule`.

    INPUT:

    - ``vector_field_module`` -- module `\mathcal{X}(U,\Phi)` of vector
      fields along `U` associated with the mapping `\Phi: U \rightarrow V`.
    - ``degree`` -- the degree of the differential form (i.e. its tensor rank)
    - ``name`` -- (default: ``None``) name given to the differential form
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the differential
      form; if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    Differential form of degree 2 on a non-parallelizable 2-dimensional manifold::

        sage: M = Manifold(2, 'M')
        sage: U = M.open_subset('U') ; V = M.open_subset('V')
        sage: M.declare_union(U,V)   # M is the union of U and V
        sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
        sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W', restrictions1= x>0, restrictions2= u+v>0)
        sage: inv = transf.inverse()
        sage: W = U.intersection(V)
        sage: eU = c_xy.frame() ; eV = c_uv.frame()
        sage: a = M.diff_form(2, name='a') ; a
        2-form 'a' on the 2-dimensional manifold 'M'
        sage: a.parent()
        Module /\^2(M) of 2-forms on the 2-dimensional manifold 'M'
        sage: a.degree()
        2

    Setting the components of a::

        sage: a[eU,0,1] = x*y^2 + 2*x
        sage: a.add_comp_by_continuation(eV, W, c_uv)
        sage: a.display(eU)
        a = (x*y^2 + 2*x) dx/\dy
        sage: a.display(eV)
        a = (-1/16*u^3 + 1/16*u*v^2 - 1/16*v^3 + 1/16*(u^2 - 8)*v - 1/2*u) du/\dv

    A 1-form on M::

        sage: a = M.one_form('a') ; a
        1-form 'a' on the 2-dimensional manifold 'M'
        sage: a.parent()
        Module /\^1(M) of 1-forms on the 2-dimensional manifold 'M'
        sage: a.degree()
        1

    Setting the components of the 1-form in a consistent way::

        sage: a[eU,:] = [-y, x]
        sage: a.add_comp_by_continuation(eV, W, c_uv)
        sage: a.display(eU)
        a = -y dx + x dy
        sage: a.display(eV)
        a = 1/2*v du - 1/2*u dv

    The exterior derivative of the 1-form is a 2-form::

        sage: da = a.exterior_der() ; da
        2-form 'da' on the 2-dimensional manifold 'M'
        sage: da.display(eU)
        da = 2 dx/\dy
        sage: da.display(eV)
        da = -du/\dv

    Another 1-form::

        sage: b = M.one_form('b')
        sage: b[eU,:] = [1+x*y, x^2]
        sage: b.add_comp_by_continuation(eV, W, c_uv)

    Adding two 1-forms results in another 1-form::

        sage: s = a + b ; s
        1-form 'a+b' on the 2-dimensional manifold 'M'
        sage: s.display(eU)
        a+b = ((x - 1)*y + 1) dx + (x^2 + x) dy
        sage: s.display(eV)
        a+b = (1/4*u^2 + 1/4*(u + 2)*v + 1/2) du + (-1/4*u*v - 1/4*v^2 - 1/2*u + 1/2) dv

    The exterior product of two 1-forms is a 2-form::

        sage: s = a.wedge(b) ; s
        2-form 'a/\b' on the 2-dimensional manifold 'M'
        sage: s.display(eU)
        a/\b = (-2*x^2*y - x) dx/\dy
        sage: s.display(eV)
        a/\b = (1/8*u^3 - 1/8*u*v^2 - 1/8*v^3 + 1/8*(u^2 + 2)*v + 1/4*u) du/\dv

    Multiplying a 1-form by a scalar field results in another 1-form::

        sage: f = M.scalar_field({c_xy: (x+y)^2, c_uv: u^2}, name='f')
        sage: s = f*a ; s
        1-form on the 2-dimensional manifold 'M'
        sage: s.display(eU)
        (-x^2*y - 2*x*y^2 - y^3) dx + (x^3 + 2*x^2*y + x*y^2) dy
        sage: s.display(eV)
        1/2*u^2*v du - 1/2*u^3 dv

    """
    def __init__(self, vector_field_module, degree, name=None, latex_name=None):
        r"""
        Construct a differential form.

        TEST::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:  intersection_name='W', restrictions1= x>0,
            ....:  restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: a = M.diff_form(2, name='a') ; a
            2-form 'a' on the 2-dimensional manifold 'M'
            sage: a[eU,0,1] = x*y^2 + 2*x
            sage: a.add_comp_by_continuation(eV, W, c_uv)

        """
        TensorField.__init__(self, vector_field_module, (0,degree), name=name,
                             latex_name=latex_name, antisym=range(degree),
                        parent=vector_field_module.dual_exterior_power(degree))
        self._init_derived() # initialization of derived quantities

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = str(self._tensor_rank) + "-form "
        if self._name is not None:
            description += "'%s' " % self._name
        return self._final_repr(description)

    def _new_instance(self):
        r"""
        Create an instance of the same class, of the same degree and on the
        same domain.
        """
        return self.__class__(self._vmodule, self._tensor_rank)

    def _init_derived(self):
        r"""
        Initialize the derived quantities.
        """
        TensorField._init_derived(self)
        self._exterior_derivative = None

    def _del_derived(self):
        r"""
        Delete the derived quantities.
        """
        TensorField._del_derived(self)
        self._exterior_derivative = None

    def exterior_der(self):
        r"""
        Compute the exterior derivative of the differential form.

        OUTPUT:

        - the exterior derivative of ``self``.

        EXAMPLE:

        """
        from sage.tensor.modules.format_utilities import format_unop_txt, \
                                                         format_unop_latex
        if self._exterior_derivative is None:
            vmodule = self._vmodule # shortcut
            rname = format_unop_txt('d', self._name)
            rlname = format_unop_latex(r'\mathrm{d}', self._latex_name)
            resu = vmodule.alternating_form(self._tensor_rank+1, name=rname,
                                            latex_name=rlname)
            for dom, rst in self._restrictions.iteritems():
                resu._restrictions[dom] = rst.exterior_der()
            self._exterior_derivative = resu
        return self._exterior_derivative

    def wedge(self, other):
        r"""
        Exterior product with another differential form.

        INPUT:

        - ``other``: another differential form

        OUTPUT:

        - instance of :class:`DiffForm` representing the exterior
          product self/\\other.

        """
        from sage.tensor.modules.free_module_alt_form import FreeModuleAltForm
        from sage.tensor.modules.format_utilities import is_atomic
        if self._domain.is_subset(other._domain):
            if not self._ambient_domain.is_subset(other._ambient_domain):
                raise TypeError("Incompatible ambient domains for exterior " +
                                "product.")
        elif other._domain.is_subset(self._domain):
            if not other._ambient_domain.is_subset(self._ambient_domain):
                raise TypeError("Incompatible ambient domains for exterior " +
                                "product.")
        dom_resu = self._domain.intersection(other._domain)
        ambient_dom_resu = self._ambient_domain.intersection(
                                                         other._ambient_domain)
        self_r = self.restrict(dom_resu)
        other_r = other.restrict(dom_resu)
        if ambient_dom_resu.is_manifestly_parallelizable():
            # call of the FreeModuleAltForm version:
            return FreeModuleAltForm.wedge(self_r, other_r)
        # otherwise, the result is created here:
        if self._name is not None and other._name is not None:
            sname = self._name
            oname = other._name
            if not is_atomic(sname):
                sname = '(' + sname + ')'
            if not is_atomic(oname):
                oname = '(' + oname + ')'
            resu_name = sname + '/\\' + oname
        if self._latex_name is not None and other._latex_name is not None:
            slname = self._latex_name
            olname = other._latex_name
            if not is_atomic(slname):
                slname = '(' + slname + ')'
            if not is_atomic(olname):
                olname = '(' + olname + ')'
            resu_latex_name = slname + r'\wedge ' + olname
        dest_map = self._vmodule._dest_map
        dest_map_resu = dest_map.restrict(dom_resu,
                                          subcodomain=ambient_dom_resu)
        vmodule = dom_resu.vector_field_module(dest_map=dest_map_resu)
        resu_degree = self._tensor_rank + other._tensor_rank
        resu = vmodule.alternating_form(resu_degree, name=resu_name,
                                        latex_name=resu_latex_name)
        for dom in self_r._restrictions:
            if dom in other_r._restrictions:
                resu._restrictions[dom] = self_r._restrictions[dom].wedge(
                                          other_r._restrictions[dom])
        return resu

    def hodge_star(self, metric):
        r"""
        Compute the Hodge dual of the differential form.

        If ``self`` is a `p`-form `A`, its Hodge dual is the `(n-p)`-form
        `*A` defined by (`n` being the manifold's dimension)

        .. MATH::

            *A_{i_1\ldots i_{n-p}} = \frac{1}{p!} A_{k_1\ldots k_p}
                \epsilon^{k_1\ldots k_p}_{\qquad\ i_1\ldots i_{n-p}}

        where `\epsilon` is the volume form associated with some
        pseudo-Riemannian metric `g` on the manifold, and the indices
        `k_1,\ldots, k_p` are raised with `g`.

        INPUT:

        - ``metric``: the pseudo-Riemannian metric `g` defining the Hodge dual,
          via the volume form `\epsilon`; must be an instance of
          :class:`~sage.geometry.manifolds.metric.Metric`

        OUTPUT:

        - the `(n-p)`-form `*A`

        EXAMPLES:

        """
        from sage.functions.other import factorial
        from sage.tensor.modules.format_utilities import format_unop_txt, \
                                                         format_unop_latex
        p = self._tensor_rank
        eps = metric.volume_form(p)
        args = range(p) + [eps] + range(p)
        resu = self.contract(*args)
        if p > 1:
            resu = resu / factorial(p)
        resu.set_name(name=format_unop_txt('*', self._name),
                     latex_name=format_unop_latex(r'\star ', self._latex_name))
        return resu

    def degree(self):
        r"""
        Return the degree of ``self``.

        """
        return self._tensor_rank


#******************************************************************************

class DiffFormParal(FreeModuleAltForm, TensorFieldParal):
    r"""
    Differential form with values in a parallelizable open subset of a
    differentiable manifold.

    Given an positive integer `p`, an open set `U` of a manifold `S`,
    a parallelizable open set `V` of a manifold `M`  and a differentiable
    mapping
    `\Phi: U \rightarrow V`, an instance of this class is a field along `U`
    of alternating multilinear forms of degree `p` in the tangent spaces to
    `V`. The standard case of a differential form *on* a manifold corresponds
    to `S=M`, `U=V` and `\Phi = \mathrm{Id}_U`. Other common cases are `\Phi`
    being an immersion and `\Phi` being a curve in `V` (`U` is then an open
    interval of `\RR`).

    If `V` is not parallelizable, the class :class:`DiffForm` must
    be used instead.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.geometry.manifolds.diffform_module.DiffFormFreeModule`.

    INPUT:

    - ``vector_field_module`` -- module `\mathcal{X}(U,\Phi)` of vector
      fields along `U` associated with the mapping `\Phi: U \rightarrow V`.
    - ``degree`` -- the degree of the differential form (i.e. its tensor rank)
    - ``name`` -- (default: ``None``) name given to the differential form
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the differential
      form; if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A 2-form on a 4-dimensional manifold::

        sage: M = Manifold(4, 'M')
        sage: c_txyz.<t,x,y,z> = M.chart()
        sage: a = M.diff_form(2, 'a') ; a
        2-form 'a' on the 4-dimensional manifold 'M'
        sage: a.parent()
        Free module /\^2(M) of 2-forms on the 4-dimensional manifold 'M'

    A differential form is a tensor field of purely covariant type::

        sage: a.tensor_type()
        (0, 2)

    It is antisymmetric, its components being instances of class
    :class:`~sage.tensor.modules.comp.CompFullyAntiSym`::

        sage: a.symmetries()
        no symmetry;  antisymmetry: (0, 1)
        sage: a[0,1] = 2
        sage: a[1,0]
        -2
        sage: a.comp()
        Fully antisymmetric 2-indices components w.r.t. coordinate frame (M, (d/dt,d/dx,d/dy,d/dz))
        sage: type(a.comp())
        <class 'sage.tensor.modules.comp.CompFullyAntiSym'>

    Setting a component with repeated indices to a non-zero value results in an
    error::

        sage: a[1,1] = 3
        Traceback (most recent call last):
        ...
        ValueError: by antisymmetry, the component cannot have a nonzero value for the indices (1, 1)
        sage: a[1,1] = 0  # OK, albeit useless
        sage: a[1,2] = 3  # OK

    The expansion of a differential form with respect to a given coframe is
    displayed via the method
    :meth:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm.display`::

        sage: a.display() # expansion with respect to the default coframe (dt, dx, dy, dz)
        a = 2 dt/\dx + 3 dx/\dy
        sage: latex(a.display()) # output for the notebook
        a = 2 \mathrm{d} t\wedge \mathrm{d} x + 3 \mathrm{d} x\wedge \mathrm{d} y

    Differential forms can be added or subtracted::

        sage: b = M.diff_form(2)
        sage: b[0,1], b[0,2], b[0,3] = (1,2,3)
        sage: s = a + b ; s
        2-form on the 4-dimensional manifold 'M'
        sage: a[:], b[:], s[:]
        (
        [ 0  2  0  0]  [ 0  1  2  3]  [ 0  3  2  3]
        [-2  0  3  0]  [-1  0  0  0]  [-3  0  3  0]
        [ 0 -3  0  0]  [-2  0  0  0]  [-2 -3  0  0]
        [ 0  0  0  0], [-3  0  0  0], [-3  0  0  0]
        )
        sage: s = a - b ; s
        2-form on the 4-dimensional manifold 'M'
        sage: s[:]
        [ 0  1 -2 -3]
        [-1  0  3  0]
        [ 2 -3  0  0]
        [ 3  0  0  0]

    An example of 3-form is the volume element on `\RR^3` in Cartesian
    coordinates::

        sage: M = Manifold(3, 'R3', '\RR^3', start_index=1)
        sage: c_cart.<x,y,z> = M.chart()
        sage: eps = M.diff_form(3, 'epsilon', r'\epsilon')
        sage: eps[1,2,3] = 1  # the only independent component
        sage: eps[:] # all the components are set from the previous line:
        [[[0, 0, 0], [0, 0, 1], [0, -1, 0]], [[0, 0, -1], [0, 0, 0], [1, 0, 0]], [[0, 1, 0], [-1, 0, 0], [0, 0, 0]]]
        sage: eps.display()
        epsilon = dx/\dy/\dz

    Spherical components of the volume element from the tensorial
    change-of-frame formula::

        sage: c_spher.<r,th,ph> = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi')
        sage: spher_to_cart = c_spher.coord_change(c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
        sage: cart_to_spher = spher_to_cart.set_inverse(sqrt(x^2+y^2+z^2), atan2(sqrt(x^2+y^2),z), atan2(y, x))
        Check of the inverse coordinate transformation:
          r == r
          th == arctan2(r*sin(th), r*cos(th))
          ph == arctan2(r*sin(ph)*sin(th), r*cos(ph)*sin(th))
          x == x
          y == y
          z == z
        sage: eps.comp(c_spher.frame()) # computation of the components in the spherical frame
        Fully antisymmetric 3-indices components w.r.t. coordinate frame (R3, (d/dr,d/dth,d/dph))
        sage: eps.comp(c_spher.frame())[1,2,3, c_spher]
        r^2*sin(th)
        sage: eps.display(c_spher.frame())
        epsilon = sqrt(x^2 + y^2 + z^2)*sqrt(x^2 + y^2) dr/\dth/\dph
        sage: eps.display(c_spher.frame(), c_spher)
        epsilon = r^2*sin(th) dr/\dth/\dph

    The exterior product of two differential forms is performed via the method
    :meth:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm.wedge`::

        sage: a = M.one_form('A')
        sage: a[:] = (x*y*z, -z*x, y*z)
        sage: b = M.one_form('B')
        sage: b[:] = (cos(z), sin(x), cos(y))
        sage: ab = a.wedge(b) ; ab
        2-form 'A/\B' on the 3-dimensional manifold 'R3'
        sage: ab[:]
        [                         0  x*y*z*sin(x) + x*z*cos(z)  x*y*z*cos(y) - y*z*cos(z)]
        [-x*y*z*sin(x) - x*z*cos(z)                          0   -(x*cos(y) + y*sin(x))*z]
        [-x*y*z*cos(y) + y*z*cos(z)    (x*cos(y) + y*sin(x))*z                          0]
        sage: ab.display()
        A/\B = (x*y*z*sin(x) + x*z*cos(z)) dx/\dy + (x*y*z*cos(y) - y*z*cos(z)) dx/\dz - (x*cos(y) + y*sin(x))*z dy/\dz

    The tensor product of a 1-form and a 2-form is not a 3-form but a tensor
    field of type (0,3) with less symmetries::

        sage: c = a*ab ; c
        tensor field 'A*(A/\B)' of type (0,3) on the 3-dimensional manifold 'R3'
        sage: c.symmetries()  #  the antisymmetry is only w.r.t. the last two arguments:
        no symmetry;  antisymmetry: (1, 2)
        sage: d = ab*a ; d
        tensor field '(A/\B)*A' of type (0,3) on the 3-dimensional manifold 'R3'
        sage: d.symmetries()  #  the antisymmetry is only w.r.t. the first two arguments:
        no symmetry;  antisymmetry: (0, 1)

    The exterior derivative of a differential form is obtained by means of the
    method :meth:`exterior_der`::

        sage: da = a.exterior_der() ; da
        2-form 'dA' on the 3-dimensional manifold 'R3'
        sage: da.display()
        dA = -(x + 1)*z dx/\dy - x*y dx/\dz + (x + z) dy/\dz
        sage: db = b.exterior_der() ; db
        2-form 'dB' on the 3-dimensional manifold 'R3'
        sage: db.display()
        dB = cos(x) dx/\dy + sin(z) dx/\dz - sin(y) dy/\dz
        sage: dab = ab.exterior_der() ; dab
        3-form 'd(A/\B)' on the 3-dimensional manifold 'R3'

    As a 3-form over a 3-dimensional manifold, d(A/\\B) is necessarily
    proportional to the volume 3-form::

        sage: dab == dab[[1,2,3]]/eps[[1,2,3]]*eps
        True

    We may also check that the classical anti-derivation formula is fulfilled::

        sage: dab == da.wedge(b) - a.wedge(db)
        True

    The Lie derivative of a 2-form is a 2-form::

        sage: v = M.vector_field('v')
        sage: v[:] = (y*z, -x*z, x*y)
        sage: ab.lie_der(v)
        2-form on the 3-dimensional manifold 'R3'

    Let us check Cartan formula, which expresses the Lie derivative in terms
    of exterior derivatives::

        sage: ab.lie_der(v) == v.contract(ab.exterior_der()) + (v.contract(ab)).exterior_der()
        True

    A 1-form on a `\RR^3`::

        sage: om = M.one_form('omega', r'\omega') ; om
        1-form 'omega' on the 3-dimensional manifold 'R3'

    A 1-form is of course a differential form::

        sage: isinstance(om, sage.geometry.manifolds.diffform.DiffFormParal)
        True
        sage: om.parent()
        Free module /\^1(R3) of 1-forms on the 3-dimensional manifold 'R3'
        sage: om.tensor_type()
        (0, 1)

    Setting the components w.r.t. the manifold's default frame::

        sage: om[:] = (2*z, x, x-y)
        sage: om[:]
        [2*z, x, x - y]
        sage: om.display()
        omega = 2*z dx + x dy + (x - y) dz

    A 1-form acts on vector fields::

        sage: v = M.vector_field('V')
        sage: v[:] = (x, 2*y, 3*z)
        sage: om(v)
        scalar field 'omega(V)' on the 3-dimensional manifold 'R3'
        sage: om(v).display()
        omega(V): R3 --> R
           (x, y, z) |--> 2*x*y + (5*x - 3*y)*z
           (r, th, ph) |--> 2*r^2*cos(ph)*sin(ph)*sin(th)^2 + r^2*(5*cos(ph)
                            - 3*sin(ph))*cos(th)*sin(th)
        sage: latex(om(v))
        \omega\left(V\right)

    The tensor product of two 1-forms is a tensor field of type (0,2)::

        sage: a = M.one_form('A')
        sage: a[:] = (1, 2, 3)
        sage: b = M.one_form('B')
        sage: b[:] = (6, 5, 4)
        sage: c = a*b ; c
        tensor field 'A*B' of type (0,2) on the 3-dimensional manifold 'R3'
        sage: c[:]
        [ 6  5  4]
        [12 10  8]
        [18 15 12]
        sage: c.symmetries()    # c has no symmetries:
        no symmetry;  no antisymmetry

    The exterior product of two 1-forms is a 2-form::

        sage: d = a.wedge(b) ; d
        2-form 'A/\B' on the 3-dimensional manifold 'R3'
        sage: d[:]
        [  0  -7 -14]
        [  7   0  -7]
        [ 14   7   0]
        sage: d.symmetries()
        no symmetry;  antisymmetry: (0, 1)

    We can check the standard formula relating the exterior product to the
    tensor product::

        sage: a.wedge(b) == a*b - b*a
        True

    """
    def __init__(self, vector_field_module, degree, name=None, latex_name=None):
        r"""
        Construct a differential form.

        TEST::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: a = M.diff_form(2, name='a') ; a
            2-form 'a' on the 2-dimensional manifold 'M'
            sage: a.parent()
            Free module /\^2(M) of 2-forms on the 2-dimensional manifold 'M'
            sage: a[0,1] = x*y
            sage: TestSuite(a).run()

        """
        FreeModuleAltForm.__init__(self, vector_field_module, degree,
                                   name=name, latex_name=latex_name)
        # TensorFieldParal attributes:
        self._vmodule = vector_field_module
        self._domain = vector_field_module._domain
        self._ambient_domain = vector_field_module._ambient_domain
        # initialization of derived quantities:
        self._init_derived()

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = str(self._tensor_rank) + "-form "
        if self._name is not None:
            description += "'%s' " % self._name
        return self._final_repr(description)

    def _new_instance(self):
        r"""
        Create an instance of the same class, of the same degree and on the
        same domain.
        """
        return self.__class__(self._fmodule, self._tensor_rank)

    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        TensorFieldParal._init_derived(self)
        self._exterior_derivative = None

    def _del_derived(self, del_restrictions=True):
        r"""
        Delete the derived quantities

        INPUT:

        - ``del_restrictions`` -- (default: True) determines whether the
          restrictions of ``self`` to subdomains are deleted.

        """
        TensorFieldParal._del_derived(self, del_restrictions=del_restrictions)
        self._exterior_derivative = None

    def __call__(self, *args):
        r"""
        Redefinition of
        :meth:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor.__call__`
        to allow for domain treatment
        """
        return TensorFieldParal.__call__(self, *args)

    def exterior_der(self):
        r"""
        Compute the exterior derivative of the differential form.

        OUTPUT:

        - the exterior derivative of ``self``.

        EXAMPLE:

        Exterior derivative of a 1-form on a 4-dimensional manifold::

            sage: M = Manifold(4, 'M')
            sage: c_txyz.<t,x,y,z> = M.chart()
            sage: a = M.one_form('A')
            sage: a[:] = (t*x*y*z, z*y**2, x*z**2, x**2 + y**2)
            sage: da = a.exterior_der() ; da
            2-form 'dA' on the 4-dimensional manifold 'M'
            sage: da.display()
            dA = -t*y*z dt/\dx - t*x*z dt/\dy - t*x*y dt/\dz + (-2*y*z + z^2) dx/\dy + (-y^2 + 2*x) dx/\dz + (-2*x*z + 2*y) dy/\dz
            sage: latex(da)
            \mathrm{d}A

        The exterior derivative is nilpotent::

            sage: dda = da.exterior_der() ; dda
            3-form 'ddA' on the 4-dimensional manifold 'M'
            sage: dda.display()
            ddA = 0
            sage: dda == 0
            True

        """
        from sage.calculus.functional import diff
        from sage.tensor.modules.format_utilities import format_unop_txt, \
                                                         format_unop_latex
        from sage.tensor.modules.comp import CompFullyAntiSym
        from vectorframe import CoordFrame
        if self._exterior_derivative is None:
            # A new computation is necessary:
            fmodule = self._fmodule # shortcut
            rname = format_unop_txt('d', self._name)
            rlname = format_unop_latex(r'\mathrm{d}', self._latex_name)
            self._exterior_derivative = fmodule.alternating_form(
                                                           self._tensor_rank+1,
                                                           name=rname,
                                                           latex_name=rlname)
            # 1/ List of all coordinate frames in which the components of self
            # are known
            coord_frames = []
            for frame in self._components:
                if isinstance(frame, CoordFrame):
                    coord_frames.append(frame)
            if coord_frames == []:
                # A coordinate frame is searched, at the price of a change of
                # frame, priveleging the frame of the domain's default chart
                dom = self._domain
                def_coordf = dom._def_chart._frame
                for frame in self._components:
                    if (frame, def_coordf) in dom._frame_changes:
                        self.comp(def_coordf, from_basis=frame)
                        coord_frames = [def_coordf]
                        break
                if coord_frames == []:
                    for chart in dom._atlas:
                        if chart != dom._def_chart: # the case def_chart is treated above
                            coordf = chart._frame
                            for frame in self._components:
                                if (frame, coordf) in dom._frame_changes:
                                    self.comp(coordf, from_basis=frame)
                                    coord_frames[coordf]
                                    break
                            if coord_frames != []:
                                break
            # 2/ The computation:
            for frame in coord_frames:
                chart = frame._chart
                sc = self._components[frame]
                dc = CompFullyAntiSym(fmodule._ring, frame,
                                      self._tensor_rank+1,
                                      start_index=fmodule._sindex,
                                     output_formatter=fmodule._output_formatter)
                for ind, val in sc._comp.iteritems():
                    for i in fmodule.irange():
                        ind_d = (i,) + ind
                        if len(ind_d) == len(set(ind_d)):
                            # all indices are different
                            dc[[ind_d]] += \
                               val.function_chart(chart).diff(i).scalar_field()
                self._exterior_derivative._components[frame] = dc
        return self._exterior_derivative

    def wedge(self, other):
        r"""
        Exterior product with another differential form.

        This is a redefinition of
        :meth:`sage.tensor.modules.free_module_alt_form.FreeModuleAltForm.wedge`
        to treat properly the domains.

        INPUT:

        - ``other``: another differential form

        OUTPUT:

        - instance of :class:`DiffFormParal` representing the exterior
          product self/\\other.

        """
        from sage.tensor.modules.free_module_alt_form import FreeModuleAltForm
        if self._domain.is_subset(other._domain):
            if not self._ambient_domain.is_subset(other._ambient_domain):
                raise TypeError("Incompatible ambient domains for exterior " +
                                "product.")
        elif other._domain.is_subset(self._domain):
            if not other._ambient_domain.is_subset(self._ambient_domain):
                raise TypeError("Incompatible ambient domains for exterior " +
                                "product.")
        dom_resu = self._domain.intersection(other._domain)
        self_r = self.restrict(dom_resu)
        other_r = other.restrict(dom_resu)
        return FreeModuleAltForm.wedge(self_r, other_r)

    def hodge_star(self, metric):
        r"""
        Compute the Hodge dual of the differential form.

        If ``self`` is a `p`-form `A`, its Hodge dual is the `(n-p)`-form
        `*A` defined by (`n` being the manifold's dimension)

        .. MATH::

            *A_{i_1\ldots i_{n-p}} = \frac{1}{p!} A_{k_1\ldots k_p}
                \epsilon^{k_1\ldots k_p}_{\qquad\ i_1\ldots i_{n-p}}

        where `\epsilon` is the volume form associated with some
        pseudo-Riemannian metric `g` on the manifold, and the indices
        `k_1,\ldots, k_p` are raised with `g`.

        INPUT:

        - ``metric``: the pseudo-Riemannian metric `g` defining the Hodge dual,
          via the volume form `\epsilon`; must be an instance of
          :class:`~sage.geometry.manifolds.metric.Metric`

        OUTPUT:

        - the `(n-p)`-form `*A`

        EXAMPLES:

        Hodge star of a 1-form in the Euclidean space `R^3`::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: g = M.metric('g')
            sage: g[1,1], g[2,2], g[3,3] = 1, 1, 1
            sage: a = M.one_form('A')
            sage: var('Ax Ay Az')
            (Ax, Ay, Az)
            sage: a[:] = (Ax, Ay, Az)
            sage: sa = a.hodge_star(g) ; sa
            2-form '*A' on the 3-dimensional manifold 'M'
            sage: sa.display()
            *A = Az dx/\dy - Ay dx/\dz + Ax dy/\dz
            sage: ssa = sa.hodge_star(g) ; ssa
            1-form '**A' on the 3-dimensional manifold 'M'
            sage: ssa.display()
            **A = Ax dx + Ay dy + Az dz
            sage: ssa == a  # must hold for a Riemannian metric in dimension 3
            True

        Hodge star of a 0-form (scalar field) in `R^3`::

            sage: f = M.scalar_field(function('F',x,y,z), name='f')
            sage: sf = f.hodge_star(g) ; sf
            3-form '*f' on the 3-dimensional manifold 'M'
            sage: sf.display()
            *f = F(x, y, z) dx/\dy/\dz
            sage: ssf = sf.hodge_star(g) ; ssf
            scalar field '**f' on the 3-dimensional manifold 'M'
            sage: ssf.display()
            **f: M --> R
               (x, y, z) |--> F(x, y, z)
            sage: ssf == f # must hold for a Riemannian metric
            True

        Hodge star of a 0-form in Minkowksi spacetime::

            sage: M = Manifold(4, 'M')
            sage: X.<t,x,y,z> = M.chart()
            sage: g = M.metric('g', signature=2)
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1, 1, 1, 1
            sage: g.display()  # Minkowski metric
            g = -dt*dt + dx*dx + dy*dy + dz*dz
            sage: var('f0')
            f0
            sage: f = M.scalar_field(f0, name='f')
            sage: sf = f.hodge_star(g) ; sf
            4-form '*f' on the 4-dimensional manifold 'M'
            sage: sf.display()
            *f = f0 dt/\dx/\dy/\dz
            sage: ssf = sf.hodge_star(g) ; ssf
            scalar field '**f' on the 4-dimensional manifold 'M'
            sage: ssf.display()
            **f: M --> R
               (t, x, y, z) |--> -f0
            sage: ssf == -f  # must hold for a Lorentzian metric
            True

        Hodge star of a 1-form in Minkowksi spacetime::

            sage: a = M.one_form('A')
            sage: var('At Ax Ay Az')
            (At, Ax, Ay, Az)
            sage: a[:] = (At, Ax, Ay, Az)
            sage: a.display()
            A = At dt + Ax dx + Ay dy + Az dz
            sage: sa = a.hodge_star(g) ; sa
            3-form '*A' on the 4-dimensional manifold 'M'
            sage: sa.display()
            *A = -Az dt/\dx/\dy + Ay dt/\dx/\dz - Ax dt/\dy/\dz - At dx/\dy/\dz
            sage: ssa = sa.hodge_star(g) ; ssa
            1-form '**A' on the 4-dimensional manifold 'M'
            sage: ssa.display()
            **A = At dt + Ax dx + Ay dy + Az dz
            sage: ssa == a  # must hold for a Lorentzian metric in dimension 4
            True

        Hodge star of a 2-form in Minkowksi spacetime::

            sage: F = M.diff_form(2, 'F')
            sage: var('Ex Ey Ez Bx By Bz')
            (Ex, Ey, Ez, Bx, By, Bz)
            sage: F[0,1], F[0,2], F[0,3] = -Ex, -Ey, -Ez
            sage: F[1,2], F[1,3], F[2,3] = Bz, -By, Bx
            sage: F[:]
            [  0 -Ex -Ey -Ez]
            [ Ex   0  Bz -By]
            [ Ey -Bz   0  Bx]
            [ Ez  By -Bx   0]
            sage: sF = F.hodge_star(g) ; sF
            2-form '*F' on the 4-dimensional manifold 'M'
            sage: sF[:]
            [  0  Bx  By  Bz]
            [-Bx   0  Ez -Ey]
            [-By -Ez   0  Ex]
            [-Bz  Ey -Ex   0]
            sage: ssF = sF.hodge_star(g) ; ssF
            2-form '**F' on the 4-dimensional manifold 'M'
            sage: ssF[:]
            [  0  Ex  Ey  Ez]
            [-Ex   0 -Bz  By]
            [-Ey  Bz   0 -Bx]
            [-Ez -By  Bx   0]
            sage: ssF.display()
            **F = Ex dt/\dx + Ey dt/\dy + Ez dt/\dz - Bz dx/\dy + By dx/\dz - Bx dy/\dz
            sage: F.display()
            F = -Ex dt/\dx - Ey dt/\dy - Ez dt/\dz + Bz dx/\dy - By dx/\dz + Bx dy/\dz
            sage: ssF == -F  # must hold for a Lorentzian metric in dimension 4
            True

        Test of the standard identity

        .. MATH::

            *(A\wedge B) = \epsilon(A^\sharp, B^\sharp, ., .)

        where `A` and `B` are any 1-forms and `A^\sharp` and `B^\sharp` the
        vectors associated to them by the metric `g` (index raising)::

            sage: b = M.one_form('B')
            sage: var('Bt Bx By Bz')
            (Bt, Bx, By, Bz)
            sage: b[:] = (Bt, Bx, By, Bz) ; b.display()
            B = Bt dt + Bx dx + By dy + Bz dz
            sage: epsilon = g.volume_form()
            sage: (a.wedge(b)).hodge_star(g) == epsilon.contract(0,a.up(g)).contract(0,b.up(g))
            True

        """
        from sage.functions.other import factorial
        from sage.tensor.modules.format_utilities import format_unop_txt, \
                                                         format_unop_latex
        p = self._tensor_rank
        eps = metric.volume_form(p)
        resu = self.contract(0, eps, 0)
        for j in range(1, p):
            resu = resu.trace(0, p-j)
        if p > 1:
            resu = resu / factorial(p)
        resu.set_name(name=format_unop_txt('*', self._name),
                     latex_name=format_unop_latex(r'\star ', self._latex_name))
        return resu
