r"""
Tensor fields with values on a parallelizable manifold

The class :class:`TensorFieldParal` implements tensor fields along a
differentiable manifolds with values on a parallelizable differentiable
manifold. For non-parallelizable manifolds, see the class
:class:`~sage.manifolds.differentiable.tensorfield.TensorField`.

Various derived classes of :class:`TensorFieldParal` are devoted to specific
tensor fields:

* :class:`~sage.manifolds.differentiable.vectorfield.VectorFieldParal` for
  vector fields (rank-1 contravariant tensor fields)

* :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismFieldParal`
  for fields of tangent-space automorphisms

* :class:`~sage.manifolds.differentiable.diff_form.DiffFormParal` for
  differential forms (fully antisymmetric covariant tensor fields)


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015) : initial version

REFERENCES:

- S. Kobayashi & K. Nomizu : *Foundations of Differential Geometry*, vol. 1,
  Interscience Publishers (New York) (1963)
- J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer (New York)
  (2013)
- B O'Neill : *Semi-Riemannian Geometry*, Academic Press (San Diego) (1983)

EXAMPLES:

A tensor field of type (1,1) on a 2-dimensional differentiable manifold::

    sage: M = Manifold(2, 'M', start_index=1)
    sage: c_xy.<x,y> = M.chart()
    sage: t = M.tensor_field(1, 1, 'T') ; t
    Tensor field T of type (1,1) on the 2-dimensional differentiable manifold M
    sage: t.tensor_type()
    (1, 1)
    sage: t.tensor_rank()
    2

Components w.r.t. the manifold's default frame are created by providing the
relevant indices inside square brackets::

    sage: t[1,1] = x^2

Unset components are initialized to zero::

    sage: t[:]  # list of components w.r.t. the manifold's default vector frame
    [x^2   0]
    [  0   0]

The full set of components w.r.t. a given vector frame is returned by the
method
:meth:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal.comp`; it is
an instance of the class :class:`~sage.tensor.modules.comp.Components`::

    sage: t.comp(c_xy.frame())
    2-indices components w.r.t. Coordinate frame (M, (d/dx,d/dy))
    sage: type(t.comp(c_xy.frame()))
    <class 'sage.tensor.modules.comp.Components'>

If no vector frame is mentionned in the argument of
:meth:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal.comp`, it is
assumed to be the manifold's default frame::

    sage: M.default_frame()
    Coordinate frame (M, (d/dx,d/dy))
    sage: t.comp() is t.comp(c_xy.frame())
    True

Individual components w.r.t. the manifold's default frame are accessed by
listing their indices inside double square brackets; they are scalar
fields on the manifold, and therefore instances of the class
:class:`~sage.manifolds.differentiable.scalarfield.DiffScalarField`::

    sage: t[[1,1]]
    Scalar field on the 2-dimensional differentiable manifold M
    sage: t[[1,1]].display()
    M --> R
    (x, y) |--> x^2
    sage: t[[1,2]]
    Scalar field zero on the 2-dimensional differentiable manifold M
    sage: t[[1,2]].display()
    zero: M --> R
       (x, y) |--> 0

A direct access to the coordinate expression of some component is obtained
via the single square brackets::

    sage: t[1,1]
    x^2
    sage: t[1,1] is t[[1,1]].coord_function() # the coordinate function
    True
    sage: t[1,1] is t[[1,1]].coord_function(c_xy)
    True
    sage: t[1,1].expr() is t[[1,1]].expr() # the symbolic expression
    True

In other words, the single square brackets return an instance of
:class:`~sage.manifolds.coord_func.CoordFunction` that is the
coordinate function representing the component in some chart (by default,
the manifold's default chart)::

    sage: type(t[1,1])    # single bracket --> coordinate function
    <class 'sage.manifolds.coord_func_symb.CoordFunctionSymb'>
    sage: type(t[[1,1]])  # double bracket --> scalar field
    <class 'sage.manifolds.differentiable.scalarfield.DiffScalarFieldAlgebra_with_category.element_class'>

Expressions in a chart different from the manifold's default one are
obtained by specifying the chart as the last argument inside the
single square brackets::

    sage: c_uv.<u,v> = M.chart()
    sage: xy_to_uv = c_xy.transition_map(c_uv, [x+y, x-y])
    sage: uv_to_xy = xy_to_uv.inverse()
    sage: t[1,1, c_uv]
    1/4*u^2 + 1/2*u*v + 1/4*v^2

Note that ``t[1,1, c_uv]`` is the component of the tensor t w.r.t. to
the coordinate frame associated to the chart (x,y) expressed in terms of
the coordinates (u,v). Indeed, ``t[1,1, c_uv]`` is a shortcut for
``t.comp(c_xy.frame())[[1,1]].coord_function(c_uv)``::

    sage: t[1,1, c_uv] is t.comp(c_xy.frame())[[1,1]].coord_function(c_uv)
    True

Similarly, ``t[1,1]`` is a shortcut for
``t.comp(c_xy.frame())[[1,1]].coord_function(c_xy)``::

    sage: t[1,1] is t.comp(c_xy.frame())[[1,1]].coord_function(c_xy)
    True
    sage: t[1,1] is t.comp()[[1,1]].coord_function()  # since c_xy.frame() and c_xy are the manifold's default values
    True

Internally, the components are stored as a dictionary (attribute
:attr:`_comp` of the class
:class:`~sage.tensor.modules.comp.Components`) whose
keys are the indices. Only the non-zero components and non-redundant
components (in case of symmetries) are stored::

    sage: t.comp()._comp
    {(1, 1): Scalar field on the 2-dimensional differentiable manifold M}

All the components can be set at once via [:]::

    sage: t[:] = [[1, -x], [x*y, 2]]
    sage: t[:]
    [  1  -x]
    [x*y   2]

To set the components in a vector frame different from the manifold's
default one, the method
:meth:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal.set_comp`
can be employed::

    sage: e = M.vector_frame('e')
    sage: t.set_comp(e)[1,1] = x+y
    sage: t.set_comp(e)[2,1], t.set_comp(e)[2,2] = y, -3*x

but, as a shortcut, one may simply specify the frame as the first argument
of the square brackets::

    sage: t[e,1,1] = x+y
    sage: t[e,2,1], t[e,2,2] = y, -3*x
    sage: t.comp(e)
    2-indices components w.r.t. Vector frame (M, (e_1,e_2))
    sage: t.comp(e)[:]
    [x + y     0]
    [    y  -3*x]
    sage: t[e,:]  # a shortcut of the above
    [x + y     0]
    [    y  -3*x]

All the components in some frame can be set at once, via the operator
[:]::

    sage: t[e,:] = [[x+y, 0], [y, -3*x]]
    sage: t[e,:]  # same as above:
    [x + y     0]
    [    y  -3*x]

To avoid any insconstency between the various components, the method
:meth:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal.set_comp`
clears the components in other frames.
To keep the other components, one must use the method
:meth:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal.add_comp`::

    sage: t = M.tensor_field(1, 1, 'T')  # Let us restart
    sage: t[:] = [[1, -x], [x*y, 2]]  # by first setting the components in the frame c_xy.frame()
    sage: # We now set the components in the frame e with add_comp:
    sage: t.add_comp(e)[:] = [[x+y, 0], [y, -3*x]]

The expansion of the tensor field in a given frame is obtained via the
method
:meth:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor.display` (the
symbol * stands for tensor product)::

    sage: t.display()  # expansion in the manifold's default frame
    T = d/dx*dx - x d/dx*dy + x*y d/dy*dx + 2 d/dy*dy
    sage: t.display(e)
    T = (x + y) e_1*e^1 + y e_2*e^1 - 3*x e_2*e^2

By definition, a tensor field acts as a multilinear map on 1-forms and vector
fields; in the present case, T being of type (1,1), it acts on pairs
(1-form, vector field)::

    sage: a = M.one_form('a')
    sage: a[:] = (1, x)
    sage: v = M.vector_field('V')
    sage: v[:] = (y, 2)
    sage: t(a,v)
    Scalar field T(a,V) on the 2-dimensional differentiable manifold M
    sage: t(a,v).display()
    T(a,V): M --> R
       (x, y) |--> x^2*y^2 + 2*x + y
       (u, v) |--> 1/16*u^4 - 1/8*u^2*v^2 + 1/16*v^4 + 3/2*u + 1/2*v
    sage: latex(t(a,v))
    T\left(a,V\right)

Check by means of the component expression of t(a,v)::

    sage: t(a,v).expr() - t[1,1]*a[1]*v[1] - t[1,2]*a[1]*v[2] \
    ....: - t[2,1]*a[2]*v[1] - t[2,2]*a[2]*v[2]
    0

A scalar field (rank-0 tensor field)::

    sage: f = M.scalar_field(x*y + 2, name='f') ; f
    Scalar field f on the 2-dimensional differentiable manifold M
    sage: f.tensor_type()
    (0, 0)

A scalar field acts on points on the manifold::

    sage: p = M.point((1,2))
    sage: f(p)
    4

See :class:`~sage.manifolds.differentiable.scalarfield.DiffScalarField` for
more details on scalar fields.

A vector field (rank-1 contravariant tensor field)::

    sage: v = M.vector_field('v') ; v
    Vector field v on the 2-dimensional differentiable manifold M
    sage: v.tensor_type()
    (1, 0)
    sage: v[1], v[2] = -x, y
    sage: v.display()
    v = -x d/dx + y d/dy

A field of symmetric bilinear forms::

    sage: q = M.sym_bilin_form_field('Q') ; q
    Field of symmetric bilinear forms Q on the 2-dimensional differentiable
     manifold M
    sage: q.tensor_type()
    (0, 2)

The components of a symmetric bilinear form are dealt by the subclass
:class:`~sage.tensor.modules.comp.CompFullySym` of the class
:class:`~sage.tensor.modules.comp.Components`, which takes into
account the symmetry between the two indices::

    sage: q[1,1], q[1,2], q[2,2] = (0, -x, y) # no need to set the component (2,1)
    sage: type(q.comp())
    <class 'sage.tensor.modules.comp.CompFullySym'>
    sage: q[:] # note that the component (2,1) is equal to the component (1,2)
    [ 0 -x]
    [-x  y]
    sage: q.display()
    Q = -x dx*dy - x dy*dx + y dy*dy

Internally (dictionary :attr:`_comp` of the class
:class:`~sage.tensor.modules.comp.Components`), only
the non-zero and non-redundant components are stored::

    sage: q.comp()._comp  # random (dictionary output)
    {(1, 2): Scalar field on the 2-dimensional differentiable manifold M,
     (2, 2): Scalar field on the 2-dimensional differentiable manifold M}
    sage: q.comp()._comp[(1,2)].expr()
    -x
    sage: q.comp()._comp[(2,2)].expr()
    y

More generally, tensor symmetries or antisymmetries can be specified via
the keywords ``sym`` and ``antisym``. For instance a rank-4 covariant
tensor symmetric with respect to its first two arguments (no. 0 and no. 1) and
antisymmetric with respect to its last two ones (no. 2 and no. 3) is declared
as follows::

    sage: t = M.tensor_field(0, 4, 'T', sym=(0,1), antisym=(2,3))
    sage: t[1,2,1,2] = 3
    sage: t[2,1,1,2] # check of the symmetry with respect to the first 2 indices
    3
    sage: t[1,2,2,1] # check of the antisymmetry with respect to the last 2 indices
    -3

"""

#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.rings.integer import Integer
from sage.structure.element import ModuleElement
from sage.tensor.modules.free_module_tensor import FreeModuleTensor
from sage.manifolds.differentiable.tensorfield import TensorField

class TensorFieldParal(FreeModuleTensor, TensorField):
    r"""
    Tensor field along a differentiable manifold, with values on a
    parallelizable manifold.

    An instance of this class is a tensor field along a differentiable
    manifold `U` with values on a parallelizable manifold `M`, via a
    differentiable map `\Phi: U \rightarrow M`. More precisely, given two
    non-negative integers `k` and `l` and a differentiable map

    .. MATH::

        \Phi:\ U \longrightarrow M,

    a *tensor field of type* `(k,l)` *along* `U` *with values on* `M` is
    a differentiable map

    .. MATH::

        t:\ U  \longrightarrow T^{(k,l)}M

    (where `T^{(k,l)}M` is the tensor bundle of type `(k,l)` over `M`) such
    that

    .. MATH::

        \forall p \in U,\ t(p) \in T^{(k,l)}(T_q M)

    i.e. `t(p)` is a tensor of type `(k,l)` on the tangent space `T_q M` at
    the point `q=\Phi(p)`, that is to say a multilinear map

    .. MATH::

        t(p):\ \underbrace{T_q^*M\times\cdots\times T_q^*M}_{k\ \; \mbox{times}}
        \times \underbrace{T_q M\times\cdots\times T_q M}_{l\ \; \mbox{times}}
        \longrightarrow K

    where `T_q^* M` is the dual vector space to `T_q M` and `K` is the
    topological field over which the manifold `M` is defined. The integer `k+l`
    is called the *tensor rank*.

    The standard case of a tensor
    field *on* a differentiable manifold corresponds to `U=M` and
    `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi` being an
    immersion and `\Phi` being a curve in `M` (`U` is then an open interval
    of `\RR`).

    If `M` is not parallelizable, the class
    :class:`~sage.manifolds.differentiable.tensorfield.TensorField`
    should be used instead.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.manifolds.differentiable.tensorfield_module.TensorFieldFreeModule`.

    INPUT:

    - ``vector_field_module`` -- free module `\mathcal{X}(U,\Phi)` of vector
      fields along `U` associated with the map `\Phi: U \rightarrow M` (cf.
      :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldFreeModule`)
    - ``tensor_type`` -- pair `(k,l)` with `k` being the contravariant rank
      and `l` the covariant rank
    - ``name`` -- (default: ``None``) name given to the tensor field
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the tensor
      field; if none is provided, the LaTeX symbol is set to ``name``
    - ``sym`` -- (default: ``None``) a symmetry or a list of symmetries among
      the tensor arguments: each symmetry is described by a tuple containing
      the positions of the involved arguments, with the convention position=0
      for the first argument. For instance:

      * ``sym=(0,1)`` for a symmetry between the 1st and 2nd arguments
      * ``sym=[(0,2),(1,3,4)]`` for a symmetry between the 1st and 3rd
        arguments and a symmetry between the 2nd, 4th and 5th arguments.

    - ``antisym`` -- (default: ``None``) antisymmetry or list of antisymmetries
      among the arguments, with the same convention as for ``sym``.


    EXAMPLES:

    A tensor field of type (2,0) on a 3-dimensional parallelizable manifold::

        sage: M = Manifold(3, 'M')
        sage: c_xyz.<x,y,z> = M.chart()  # makes M parallelizable
        sage: t = M.tensor_field(2, 0, 'T') ; t
        Tensor field T of type (2,0) on the 3-dimensional differentiable
         manifold M

    Tensor fields are considered as elements of a module over the ring
    `C^k(M)` of scalar fields on `M`::

        sage: t.parent()
        Free module T^(2,0)(M) of type-(2,0) tensors fields on the
         3-dimensional differentiable manifold M
        sage: t.parent().base_ring()
        Algebra of differentiable scalar fields on the 3-dimensional
         differentiable manifold M

    The components with respect to the manifold's default frame are set or read
    by means of square brackets::

        sage: e = M.vector_frame('e') ; M.set_default_frame(e)
        sage: for i in M.irange():
        ....:     for j in M.irange():
        ....:         t[i,j] = (i+1)**(j+1)
        ....:
        sage: [[ t[i,j] for j in M.irange()] for i in M.irange()]
        [[1, 1, 1], [2, 4, 8], [3, 9, 27]]

    A shortcut for the above is using [:]::

        sage: t[:]
        [ 1  1  1]
        [ 2  4  8]
        [ 3  9 27]

    The components with respect to another frame are set via the method
    :meth:`set_comp` and read via the method :meth:`comp`; both return an
    instance of :class:`~sage.tensor.modules.comp.Components`::

        sage: f = M.vector_frame('f')  # a new frame defined on M, in addition to e
        sage: t.set_comp(f)[0,0] = -3
        sage: t.comp(f)
        2-indices components w.r.t. Vector frame (M, (f_0,f_1,f_2))
        sage: t.comp(f)[0,0]
        -3
        sage: t.comp(f)[:]  # the full list of components
        [-3  0  0]
        [ 0  0  0]
        [ 0  0  0]

    To avoid any insconstency between the various components, the method
    :meth:`set_comp` deletes the components in other frames.
    Accordingly, the components in the frame e have been deleted::

        sage: t._components
        {Vector frame (M, (f_0,f_1,f_2)): 2-indices components w.r.t. Vector
         frame (M, (f_0,f_1,f_2))}

    To keep the other components, one must use the method :meth:`add_comp`::

        sage: t = M.tensor_field(2, 0, 'T')  # Let us restart
        sage: t[0,0] = 2                     # sets the components in the frame e
        sage: # We now set the components in the frame f with add_comp:
        sage: t.add_comp(f)[0,0] = -3
        sage: # The components w.r.t. frame e have been kept:
        sage: t._components  # random (dictionary output)
        {Vector frame (M, (e_0,e_1,e_2)): 2-indices components w.r.t. Vector frame (M, (e_0,e_1,e_2)),
         Vector frame (M, (f_0,f_1,f_2)): 2-indices components w.r.t. Vector frame (M, (f_0,f_1,f_2))}

    The basic properties of a tensor field are::

        sage: t.domain()
        3-dimensional differentiable manifold M
        sage: t.tensor_type()
        (2, 0)

    Internally, the components w.r.t. various vector frames are stored in the
    dictionary :attr:`_components`::

        sage: t._components  # random (dictionary output)
        {Vector frame (M, (e_0,e_1,e_2)): 2-indices components w.r.t. Vector frame (M, (e_0,e_1,e_2)),
         Vector frame (M, (f_0,f_1,f_2)): 2-indices components w.r.t. Vector frame (M, (f_0,f_1,f_2))}

    Symmetries and antisymmetries are declared via the keywords ``sym`` and
    ``antisym``. For instance, a rank-6 covariant tensor that is symmetric with
    respect to its 1st and 3rd arguments and antisymmetric with respect to the
    2nd, 5th and 6th arguments is set up as follows::

        sage: a = M.tensor_field(0, 6, 'T', sym=(0,2), antisym=(1,4,5))
        sage: a[0,0,1,0,1,2] = 3
        sage: a[1,0,0,0,1,2] # check of the symmetry
        3
        sage: a[0,1,1,0,0,2], a[0,1,1,0,2,0] # check of the antisymmetry
        (-3, 3)

    Multiple symmetries or antisymmetries are allowed; they must then be
    declared as a list. For instance, a rank-4 covariant tensor that is
    antisymmetric with respect to its 1st and 2nd arguments and with respect to
    its 3rd and 4th argument must be declared as::

        sage: r = M.tensor_field(0, 4, 'T', antisym=[(0,1), (2,3)])
        sage: r[0,1,2,0] = 3
        sage: r[1,0,2,0] # first antisymmetry
        -3
        sage: r[0,1,0,2] # second antisymmetry
        -3
        sage: r[1,0,0,2] # both antisymmetries acting
        3

    Tensor fields of the same type can be added and subtracted::

        sage: a = M.tensor_field(2, 0)
        sage: a[0,0], a[0,1], a[0,2] = (1,2,3)
        sage: b = M.tensor_field(2, 0)
        sage: b[0,0], b[1,1], b[2,2], b[0,2] = (4,5,6,7)
        sage: s = a + 2*b ; s
        Tensor field of type (2,0) on the 3-dimensional differentiable
         manifold M
        sage: a[:], (2*b)[:], s[:]
        (
        [1 2 3]  [ 8  0 14]  [ 9  2 17]
        [0 0 0]  [ 0 10  0]  [ 0 10  0]
        [0 0 0], [ 0  0 12], [ 0  0 12]
        )
        sage: s = a - b ; s
        Tensor field of type (2,0) on the 3-dimensional differentiable
         manifold M
        sage: a[:], b[:], s[:]
        (
        [1 2 3]  [4 0 7]  [-3  2 -4]
        [0 0 0]  [0 5 0]  [ 0 -5  0]
        [0 0 0], [0 0 6], [ 0  0 -6]
        )

    Symmetries are preserved by the addition whenever it is possible::

        sage: a = M.tensor_field(2, 0, sym=(0,1))
        sage: a[0,0], a[0,1], a[0,2] = (1,2,3)
        sage: s = a + b
        sage: a[:], b[:], s[:]
        (
        [1 2 3]  [4 0 7]  [ 5  2 10]
        [2 0 0]  [0 5 0]  [ 2  5  0]
        [3 0 0], [0 0 6], [ 3  0  6]
        )
        sage: a.symmetries()
        symmetry: (0, 1);  no antisymmetry
        sage: b.symmetries()
        no symmetry;  no antisymmetry
        sage: s.symmetries()
        no symmetry;  no antisymmetry
        sage: # let us now make b symmetric:
        sage: b = M.tensor_field(2, 0, sym=(0,1))
        sage: b[0,0], b[1,1], b[2,2], b[0,2] = (4,5,6,7)
        sage: s = a + b
        sage: a[:], b[:], s[:]
        (
        [1 2 3]  [4 0 7]  [ 5  2 10]
        [2 0 0]  [0 5 0]  [ 2  5  0]
        [3 0 0], [7 0 6], [10  0  6]
        )
        sage: s.symmetries()  # s is symmetric because both a and b are
        symmetry: (0, 1);  no antisymmetry

    The tensor product is taken with the operator \*::

        sage: c = a*b ; c
        Tensor field of type (4,0) on the 3-dimensional differentiable
         manifold M
        sage: c.symmetries()  # since a and b are both symmetric, a*b has two symmetries:
        symmetries: [(0, 1), (2, 3)];  no antisymmetry

    The tensor product of two fully contravariant tensors is not symmetric in
    general::

        sage: a*b == b*a
        False

    The tensor product of a fully contravariant tensor by a fully covariant one
    is symmetric::

        sage: d = M.diff_form(2)  # a fully covariant tensor field
        sage: d[0,1], d[0,2], d[1,2] = (3, 2, 1)
        sage: s = a*d ; s
        Tensor field of type (2,2) on the 3-dimensional differentiable
         manifold M
        sage: s.symmetries()
        symmetry: (0, 1);  antisymmetry: (2, 3)
        sage: s1 = d*a ; s1
        Tensor field of type (2,2) on the 3-dimensional differentiable
         manifold M
        sage: s1.symmetries()
        symmetry: (0, 1);  antisymmetry: (2, 3)
        sage: d*a == a*d
        True

    Example of tensor field associated with a non-trivial differentiable
    map `\Phi`: tensor field along a curve in `M`::

        sage: R = Manifold(1, 'R')  # R as a 1-dimensional manifold
        sage: T.<t> = R.chart()  # canonical chart on R
        sage: Phi = R.diff_map(M, [cos(t), sin(t), t], name='Phi') ; Phi
        Differentiable map Phi from the 1-dimensional differentiable manifold R
         to the 3-dimensional differentiable manifold M
        sage: h = R.tensor_field(2, 0, name='h', dest_map=Phi) ; h
        Tensor field h of type (2,0) along the 1-dimensional differentiable
         manifold R with values on the 3-dimensional differentiable manifold M
        sage: h.parent()
        Free module T^(2,0)(R,Phi) of type-(2,0) tensors fields along the
         1-dimensional differentiable manifold R mapped into the 3-dimensional
         differentiable manifold M
        sage: h[0,0], h[0,1], h[2,0] = 1+t, t^2, sin(t)
        sage: h.display()
        h = (t + 1) d/dx*d/dx + t^2 d/dx*d/dy + sin(t) d/dz*d/dx

    """
    def __init__(self, vector_field_module, tensor_type, name=None,
                 latex_name=None, sym=None, antisym=None):
        r"""
        Construct a tensor field.

        TESTS:

        Construction via ``parent.element_class``, and not via a direct call
        to ``TensorFieldParal``, to fit with the category framework::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: XM = M.vector_field_module()
            sage: T02 = M.tensor_field_module((0,2))
            sage: t = T02.element_class(XM, (0,2), name='t'); t
            Tensor field t of type (0,2) on the 2-dimensional differentiable
             manifold M
            sage: t[:] = [[1+x^2, x*y], [0, 1+y^2]]
            sage: t.display()
            t = (x^2 + 1) dx*dx + x*y dx*dy + (y^2 + 1) dy*dy
            sage: t.parent()
            Free module T^(0,2)(M) of type-(0,2) tensors fields on the
             2-dimensional differentiable manifold M
            sage: TestSuite(t).run()

        """
        FreeModuleTensor.__init__(self, vector_field_module, tensor_type,
                                  name=name, latex_name=latex_name,
                                  sym=sym, antisym=antisym)
        # TensorField attributes:
        self._vmodule = vector_field_module
        self._domain = vector_field_module._domain
        self._ambient_domain = vector_field_module._ambient_domain
        # NB: the TensorField attribute self._restrictions is considered as a
        #     derived quantity in the present case (the primary attribute
        #     being self._components, which is initialized by
        #     FreeModuleTensor.__init__ ); accordingly self._restrictions is
        #     initialized by _init_derived() and cleared by _del_derived().
        # Initialization of derived quantities:
        self._init_derived()

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: t = M.tensor_field(1,1, name='t')
            sage: t._repr_()
            'Tensor field t of type (1,1) on the 2-dimensional differentiable manifold M'
            sage: repr(t)  # indirect doctest
            'Tensor field t of type (1,1) on the 2-dimensional differentiable manifold M'
            sage: t  # indirect doctest
            Tensor field t of type (1,1) on the 2-dimensional differentiable
             manifold M

        """
        return TensorField._repr_(self)

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self`` on the same
        vector field module, with the same tensor type and same symmetries.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: t = M.tensor_field(1,1, name='t')
            sage: t._new_instance()
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: type(t._new_instance()) is type(t)
            True

        """
        return type(self)(self._fmodule, self._tensor_type, sym=self._sym,
                          antisym=self._antisym)

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TEST::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: t = M.tensor_field(1,1, name='t')
            sage: t._init_derived()

        """
        FreeModuleTensor._init_derived(self)
        TensorField._init_derived(self)
        self._restrictions = {} # dict. of restrictions of self on subdomains
                                # of self._domain, with the subdomains as keys

    def _del_derived(self, del_restrictions=True):
        r"""
        Delete the derived quantities

        INPUT:

        - ``del_restrictions`` -- (default: ``True``) determines whether the
          restrictions of ``self`` to subdomains are deleted.

        TEST::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: t = M.tensor_field(1,1, name='t')
            sage: t._del_derived()

        """
        FreeModuleTensor._del_derived(self)
        TensorField._del_derived(self)
        if del_restrictions:
            self._restrictions.clear()

    def set_comp(self, basis=None):
        r"""
        Return the components of the tensor field in a given vector frame
        for assignment.

        The components with respect to other frames on the same domain are
        deleted, in order to avoid any inconsistency. To keep them, use the
        method :meth:`add_comp` instead.

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame in which the components
          are defined; if none is provided, the components are assumed to refer
          to the tensor field domain's default frame.

        OUTPUT:

        - components in the given frame, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`; if such
          components did not exist previously, they are created.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: e_xy = X.frame()
            sage: t = M.tensor_field(1,1, name='t')
            sage: t.set_comp(e_xy)
            2-indices components w.r.t. Coordinate frame (M, (d/dx,d/dy))
            sage: t.set_comp(e_xy)[1,0] = 2
            sage: t.display(e_xy)
            t = 2 d/dy*dx

        Setting components in a new frame (``e``)::

            sage: e = M.vector_frame('e')
            sage: t.set_comp(e)
            2-indices components w.r.t. Vector frame (M, (e_0,e_1))
            sage: t.set_comp(e)[0,1] = x
            sage: t.display(e)
            t = x e_0*e^1

        The components w.r.t. the frame ``e_xy`` have be erased::

            sage: t.display(e_xy)
            Traceback (most recent call last):
            ...
            ValueError: no basis could be found for computing the components
             in the Coordinate frame (M, (d/dx,d/dy))

        Setting components in a frame defined on a subdomain deletes previously
        defined components as well::

            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: f = U.vector_frame('f')
            sage: t.set_comp(f)
            2-indices components w.r.t. Vector frame (U, (f_0,f_1))
            sage: t.set_comp(f)[0,1] = 1+y
            sage: t.display(f)
            t = (y + 1) f_0*f^1
            sage: t.display(e)
            Traceback (most recent call last):
            ...
            ValueError: no basis could be found for computing the components
             in the Vector frame (M, (e_0,e_1))

        """
        if basis is None:
            basis = self._fmodule._def_basis
        if basis._domain == self._domain:
            # Setting components on the tensor field domain:
            return FreeModuleTensor.set_comp(self, basis=basis)
        else:
            # Setting components on a subdomain:
            #
            # Creating or saving the restriction to the subdomain:
            rst = self.restrict(basis._domain, dest_map=basis._dest_map)
            # Deleting all the components on self._domain and the derived
            # quantities:
            self._components.clear()
            self._del_derived()
            # Restoring the restriction to the subdomain (which has been
            # deleted by _del_derived):
            self._restrictions[basis._domain] = rst
            # The set_comp operation is performed on the subdomain:
            return rst.set_comp(basis=basis)

    def add_comp(self, basis=None):
        r"""
        Return the components of the tensor field in a given vector frame
        for assignment.

        The components with respect to other frames on the same domain are
        kept. To delete them them, use the method :meth:`set_comp` instead.

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame in which the components
          are defined; if none is provided, the components are assumed to refer
          to the tensor field domain's default frame.

        OUTPUT:

        - components in the given frame, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`; if such
          components did not exist previously, they are created.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: e_xy = X.frame()
            sage: t = M.tensor_field(1,1, name='t')
            sage: t.add_comp(e_xy)
            2-indices components w.r.t. Coordinate frame (M, (d/dx,d/dy))
            sage: t.add_comp(e_xy)[1,0] = 2
            sage: t.display(e_xy)
            t = 2 d/dy*dx

        Adding components w.r.t. a new frame (``e``)::

            sage: e = M.vector_frame('e')
            sage: t.add_comp(e)
            2-indices components w.r.t. Vector frame (M, (e_0,e_1))
            sage: t.add_comp(e)[0,1] = x
            sage: t.display(e)
            t = x e_0*e^1

        The components w.r.t. frame ``e_xy`` are kept::

            sage: t.display(e_xy)
            t = 2 d/dy*dx

        Adding components in a frame defined on a subdomain::

            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: f = U.vector_frame('f')
            sage: t.add_comp(f)
            2-indices components w.r.t. Vector frame (U, (f_0,f_1))
            sage: t.add_comp(f)[0,1] = 1+y
            sage: t.display(f)
            t = (y + 1) f_0*f^1

        The components previously defined are kept::

            sage: t.display(e_xy)
            t = 2 d/dy*dx
            sage: t.display(e)
            t = x e_0*e^1

        """
        if basis is None:
            basis = self._fmodule._def_basis
        if basis._domain == self._domain:
            # Adding components on the tensor field domain:
            return FreeModuleTensor.add_comp(self, basis=basis)
        else:
            # Adding components on a subdomain:
            #
            # Creating or saving the restriction to the subdomain:
            rst = self.restrict(basis._domain, dest_map=basis._dest_map)
            # Deleting the derived quantities except for the restrictions to
            # subdomains:
            self._del_derived(del_restrictions=False)
            # The add_comp operation is performed on the subdomain:
            return rst.add_comp(basis=basis)

    def comp(self, basis=None, from_basis=None):
        r"""
        Return the components in a given vector frame.

        If the components are not known already, they are computed by the
        tensor change-of-basis formula from components in another vector frame.

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame in which the components
          are required; if none is provided, the components are assumed to
          refer to the tensor field domain's default frame
        - ``from_basis`` -- (default: ``None``) vector frame from which the
          required components are computed, via the tensor change-of-basis
          formula, if they are not known already in the basis ``basis``

        OUTPUT:

        - components in the vector frame ``basis``, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`

        EXAMPLES::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: t = M.tensor_field(1,2, name='t')
            sage: t[1,2,1] = x*y
            sage: t.comp(X.frame())
            3-indices components w.r.t. Coordinate frame (M, (d/dx,d/dy))
            sage: t.comp()  # the default frame is X.frame()
            3-indices components w.r.t. Coordinate frame (M, (d/dx,d/dy))
            sage: t.comp()[:]
            [[[0, 0], [x*y, 0]], [[0, 0], [0, 0]]]
            sage: e = M.vector_frame('e')
            sage: t[e, 2,1,1] = x-3
            sage: t.comp(e)
            3-indices components w.r.t. Vector frame (M, (e_1,e_2))
            sage: t.comp(e)[:]
            [[[0, 0], [0, 0]], [[x - 3, 0], [0, 0]]]

        """
        if basis is None:
            basis = self._fmodule._def_basis
        if basis._domain == self._domain:
            # components on the tensor field domain:
            return FreeModuleTensor.comp(self, basis=basis,
                                         from_basis=from_basis)
        else:
            # components on a subdomain:
            rst = self.restrict(basis._domain, dest_map=basis._dest_map)
            return rst.comp(basis=basis, from_basis=from_basis)


    def _common_coord_frame(self, other):
        r"""
        Find a common coordinate frame for the components of ``self`` and
        ``other``.

        In case of multiple common bases, the domain's default coordinate
        basis is privileged.
        If the current components of ``self`` and ``other`` are all relative to
        different frames, a common frame is searched by performing a component
        transformation, via the transformations listed in
        ``self._domain._frame_changes``, still privileging transformations to
        the domain's default frame.

        INPUT:

        - ``other`` -- a tensor field (instance of :class:`TensorFieldParal`)

        OUPUT:

        - common coordinate frame; if no common basis is found, None is
          returned.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: a = M.tensor_field(1,2, name='a')
            sage: a[0,1,0] = 2
            sage: b = M.vector_field(name='b')
            sage: b[:] = [-y, x]
            sage: a._common_coord_frame(b)
            Coordinate frame (M, (d/dx,d/dy))

        Vector field defined on a new chart::

            sage: Y.<u,v> = M.chart()
            sage: c = M.vector_field(name='c')
            sage: c[Y.frame(), :, Y] = (1+u, u*v)
            sage: c.display(Y.frame(), Y)
            c = (u + 1) d/du + u*v d/dv

        There is no common ccordinate frame::

            sage: a._common_coord_frame(c)

        Connecting the two coordinate charts enables to find a common frame::

            sage: X_to_Y = X.transition_map(Y, [x+y, x-y])
            sage: Y_to_X = X_to_Y.inverse()
            sage: a._common_coord_frame(c)
            Coordinate frame (M, (d/dx,d/dy))

        Indeed, the components of ``c`` w.r.t. the frame ``(M, (d/dx,d/dy))``
        have been computed via the change-of-coordinate formulas::

            sage: c.display(a._common_coord_frame(c))
            c = (1/2*x^2 - 1/2*y^2 + 1/2*x + 1/2*y + 1/2) d/dx
             + (-1/2*x^2 + 1/2*y^2 + 1/2*x + 1/2*y + 1/2) d/dy

        """
        from sage.manifolds.differentiable.vectorframe import CoordFrame
        # Compatibility checks:
        if not isinstance(other, TensorFieldParal):
            raise TypeError("the argument must be of type TensorFieldParal")
        dom = self._domain
        def_frame = dom._def_frame
        #
        # 1/ Search for a common frame among the existing components, i.e.
        #    without performing any component transformation.
        #    -------------------------------------------------------------
        # 1a/ Direct search
        if def_frame in self._components and \
           def_frame in other._components and \
           isinstance(dom._def_frame, CoordFrame):
            return def_frame # the domain's default frame is privileged
        for frame1 in self._components:
            if frame1 in other._components and \
               isinstance(frame1, CoordFrame):
                return frame1
        # 1b/ Search involving subframes
        dom2 = other._domain
        for frame1 in self._components:
            if not isinstance(frame1, CoordFrame):
                continue
            for frame2 in other._components:
                if not isinstance(frame2, CoordFrame):
                    continue
                if frame2 in frame1._subframes:
                    self.comp(frame2)
                    return frame2
                if frame1 in frame2._subframes:
                    other.comp(frame1)
                    return frame1
        #
        # 2/ Search for a common frame via one component transformation
        #    ----------------------------------------------------------
        # If this point is reached, it is indeed necessary to perform at least
        # one component transformation to get a common frame
        if isinstance(dom._def_frame, CoordFrame):
            if def_frame in self._components:
                for oframe in other._components:
                    if (oframe, def_frame) in dom._frame_changes:
                        other.comp(def_frame, from_basis=oframe)
                        return def_frame
            if def_frame in other._components:
                for sframe in self._components:
                    if (sframe, def_frame) in dom._frame_changes:
                        self.comp(def_frame, from_basis=sframe)
                        return def_frame
        # If this point is reached, then def_frame cannot be a common frame
        # via a single component transformation
        for sframe in self._components:
            if not isinstance(sframe, CoordFrame):
                continue
            for oframe in other._components:
                if not isinstance(oframe, CoordFrame):
                    continue
                if (oframe, sframe) in dom._frame_changes:
                    other.comp(sframe, from_basis=oframe)
                    return sframe
                if (sframe, oframe) in dom._frame_changes:
                    self.comp(oframe, from_basis=sframe)
                    return oframe
        #
        # 3/ Search for a common frame via two component transformations
        #    -----------------------------------------------------------
        # If this point is reached, it is indeed necessary to perform at least
        # two component transformations to get a common frame
        for sframe in self._components:
            for oframe in other._components:
                if (sframe, def_frame) in dom._frame_changes and \
                   (oframe, def_frame) in dom._frame_changes and \
                   isinstance(def_frame, CoordFrame):
                    self.comp(def_frame, from_basis=sframe)
                    other.comp(def_frame, from_basis=oframe)
                    return def_frame
                for frame in dom._frames:
                    if (sframe, frame) in dom._frame_changes and \
                       (oframe, frame) in dom._frame_changes and \
                       isinstance(frame, CoordFrame):
                        self.comp(frame, from_basis=sframe)
                        other.comp(frame, from_basis=oframe)
                        return frame
        #
        # If this point is reached, no common frame could be found, even at
        # the price of component transformations:
        return None


    def lie_der(self, vector):
        r"""
        Compute the Lie derivative with respect to a vector field.

        The Lie derivative is stored in the dictionary
        :attr:`_lie_derivatives`, so that there is no need to
        recompute it at the next call if neither ``self`` nor ``vector``
        have been modified meanwhile.

        INPUT:

        - ``vector`` -- vector field with respect to which the Lie derivative
          is to be taken

        OUTPUT:

        - the tensor field that is the Lie derivative of ``self`` with respect
          to ``vector``

        EXAMPLES:

        Lie derivative of a vector::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: c_xy.<x,y> = M.chart()
            sage: v = M.vector_field('v')
            sage: v[:] = (-y, x)
            sage: w = M.vector_field()
            sage: w[:] = (2*x+y, x*y)
            sage: w.lie_der(v)
            Vector field on the 2-dimensional differentiable manifold M
            sage: w.lie_der(v).display()
            ((x - 2)*y + x) d/dx + (x^2 - y^2 - 2*x - y) d/dy

        The Lie derivative is antisymmetric::

            sage: w.lie_der(v) == -v.lie_der(w)
            True

        For vectors, it coincides with the commutator::

            sage: f = M.scalar_field(x^3 + x*y^2)
            sage: w.lie_der(v)(f).display()
            M --> R
            (x, y) |--> -(x + 2)*y^3 + 3*x^3 - x*y^2 + 5*(x^3 - 2*x^2)*y
            sage: w.lie_der(v)(f) == v(w(f)) - w(v(f))  # rhs = commutator [v,w] acting on f
            True

        Lie derivative of a 1-form::

            sage: om = M.one_form()
            sage: om[:] = (y^2*sin(x), x^3*cos(y))
            sage: om.lie_der(v)
            1-form on the 2-dimensional differentiable manifold M
            sage: om.lie_der(v).display()
            (-y^3*cos(x) + x^3*cos(y) + 2*x*y*sin(x)) dx
             + (-x^4*sin(y) - 3*x^2*y*cos(y) - y^2*sin(x)) dy

        Check of Cartan identity::

            sage: om.lie_der(v) == v.contract(0,om.exterior_der(),0) + (om(v)).exterior_der()
            True

        """
        if vector._tensor_type != (1,0):
            raise TypeError("the argument must be a vector field")
        if id(vector) not in self._lie_derivatives:
            # A new computation must be performed
            #
            # 1/ Search for a common coordinate frame:
            coord_frame = self._common_coord_frame(vector)
            if coord_frame is None:
                raise ValueError("no common coordinate frame found")
            chart = coord_frame._chart
            #
            # 2/ Component computation:
            tc = self._components[coord_frame]
            vc = vector._components[coord_frame]
            # the result has the same tensor type and same symmetries as self:
            resc = self._new_comp(coord_frame)
            n_con = self._tensor_type[0]
            vf_module = vector._fmodule
            for ind in resc.non_redundant_index_generator():
                rsum = 0
                for i in vf_module.irange():
                    rsum += vc[[i]].coord_function(chart) * \
                           tc[[ind]].coord_function(chart).diff(i)
                # loop on contravariant indices:
                for k in range(n_con):
                    for i in vf_module.irange():
                        indk = list(ind)
                        indk[k] = i
                        rsum -= tc[[indk]].coord_function(chart) * \
                                vc[[ind[k]]].coord_function(chart).diff(i)
                # loop on covariant indices:
                for k in range(n_con, self._tensor_rank):
                    for i in vf_module.irange():
                        indk = list(ind)
                        indk[k] = i
                        rsum += tc[[indk]].coord_function(chart) * \
                                vc[[i]].coord_function(chart).diff(ind[k])
                resc[[ind]] = rsum.scalar_field()
            #
            # 3/ Final result (the tensor)
            res = vf_module.tensor_from_comp(self._tensor_type, resc)
            self._lie_derivatives[id(vector)] = (vector, res)
            vector._lie_der_along_self[id(self)] = self
        return self._lie_derivatives[id(vector)][1]

    def restrict(self, subdomain, dest_map=None):
        r"""
        Return the restriction of the tensor field to some subdomain.

        If the restriction has not been defined yet, it is constructed here.

        INPUT:

        - ``subdomain`` -- open subset `U` of the tensor field domain `S`
          (must be an instance of
          :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`)
        - ``dest_map`` -- (default: ``None``) destination map
          `\Psi:\ U \rightarrow V`, where `V` is an open subset of the manifold
          `M` where the tensor field takes it values
          (must be an instance of
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`).
          If ``None``, the restriction of `\Phi` to `U` is used, `\Phi` being
          the differentiable map `S\rightarrow M` associated with the tensor
          field.

        OUTPUT:

        - instance of :class:`TensorFieldParal` representing the restriction.

        EXAMPLES:

        Restriction of a vector field defined on `\RR^2` to a disk::

            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: v = M.vector_field('v')
            sage: v[:] = [x+y, -1+x^2]
            sage: D = M.open_subset('D') # the unit open disc
            sage: c_cart_D = c_cart.restrict(D, x^2+y^2<1)
            sage: v_D = v.restrict(D) ; v_D
            Vector field v on the Open subset D of the 2-dimensional
             differentiable manifold R^2
            sage: v_D.display()
            v = (x + y) d/dx + (x^2 - 1) d/dy

        The symbolic expressions of the components w.r.t. Cartesian coordinates
        are equal::

            sage: bool( v_D[1].expr() == v[1].expr() )
            True

        but neither the chart functions representing the components (they are
        defined on different charts)::

            sage: v_D[1] == v[1]
            False

        nor the scalar fields representing the components (they are
        defined on different open subsets)::

            sage: v_D[[1]] == v[[1]]
            False

        The restriction of the vector field to its own domain is of course
        itself::

            sage: v.restrict(M) is v
            True

        """
        if subdomain == self._domain and \
                    (dest_map is None or dest_map == self._vmodule._dest_map) :
            return self
        if subdomain not in self._restrictions:
            if not subdomain.is_subset(self._domain):
                raise ValueError("the provided domain is not a subset of " +
                                 "the field's domain.")
            if dest_map is None:
                dest_map = self._fmodule._dest_map.restrict(subdomain)
            elif not dest_map._codomain.is_subset(self._ambient_domain):
                raise ValueError("the argument 'dest_map' is not compatible " +
                                 "with the ambient domain of " +
                                 "the {}".format(self))
            #!# First one tries to derive the restriction from a tighter domain:
            #for dom, rst in self._restrictions.iteritems():
            #    if subdomain.is_subset(dom):
            #        self._restrictions[subdomain] = rst.restrict(subdomain)
            #        break
            # If this fails, the restriction is created from scratch:
            #else:
            smodule = subdomain.vector_field_module(dest_map=dest_map)
            resu = smodule.tensor(self._tensor_type, name=self._name,
                                  latex_name=self._latex_name, sym=self._sym,
                                  antisym=self._antisym,
                                  specific_type=type(self))
            for frame in self._components:
                for sframe in subdomain._covering_frames:
                    if sframe in frame._subframes:
                        comp_store = self._components[frame]._comp
                        scomp = resu._new_comp(sframe)
                        scomp_store = scomp._comp
                        # the components of the restriction are evaluated
                        # index by index:
                        for ind, value in comp_store.iteritems():
                            scomp_store[ind] = value.restrict(subdomain)
                        resu._components[sframe] = scomp
            self._restrictions[subdomain] = resu
        return self._restrictions[subdomain]

    def __call__(self, *args):
        r"""
        The tensor field acting on 1-forms and vector fields as a multilinear
        map.

        In the particular case of tensor field of type (1,1), the action can
        be on a single vector field, the tensor field being identified to a
        field of tangent-space endomorphisms. The output is then a vector
        field.

        INPUT:

        - ``*args`` -- list of k 1-forms and l vector fields, ``self``
          being a tensor of type (k,l).

        OUTPUT:

        - either the scalar field resulting from the action of ``self`` on
          the 1-forms and vector fields passed as arguments or the vector
          field resulting from the action of ``self`` as a field of
          tangent-space endomorphisms (case of a type-(1,1) tensor field)

        TESTS:

        Action of a tensor field of type (1,1)::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: t = M.tensor_field(1,1, name='t')
            sage: t[:] = [[1+x, 2], [y, -x^2]]
            sage: v = M.vector_field(name='v')
            sage: v[:] = [-y, x]
            sage: a = M.one_form(name='a')
            sage: a[:] = [3, 1-y]
            sage: s = t.__call__(a,v); s
            Scalar field t(a,v) on the 2-dimensional differentiable manifold M
            sage: s.display()
            t(a,v): M --> R
               (x, y) |--> -x^3 + y^3 + (x^3 - 3*x - 3)*y - y^2 + 6*x
            sage: s.coord_function() == sum(sum(t[i,j]*a[i]*v[j] for j in [0..1])
            ....:                           for i in [0..1])
            True
            sage: s == t(a,v)  # indirect doctest
            True

        The tensor field acting on vector field, as a field of tangent-space
        endomorphisms::

            sage: s = t.__call__(v); s
            Vector field t(v) on the 2-dimensional differentiable manifold M
            sage: s.display()
            t(v) = (-(x + 1)*y + 2*x) d/dx + (-x^3 - y^2) d/dy
            sage: s[0] == t[0,0]*v[0] + t[0,1]*v[1]
            True
            sage: s[1] == t[1,0]*v[0] + t[1,1]*v[1]
            True
            sage: s == t(v)  # indirect doctest
            True

        """
        from sage.categories.homset import End
        p = len(args)
        if p == 1 and self._tensor_type == (1,1):
            # type-(1,1) tensor acting as an endomorphism:
            vector = args[0]
            if vector._tensor_type != (1,0):
                raise TypeError("the argument must be a vector field")
            dom = self._domain.intersection(vector._domain)
            sd = self.restrict(dom)
            vd = vector.restrict(dom)
            endom = End(vd.parent())(sd)
            return endom(vd)
        # Generic case
        if p != self._tensor_rank:
            raise TypeError("{} arguments must be ".format(self._tensor_rank) +
                            "provided")
        # Domain of the result
        dom_resu = self._domain
        for arg in args:
            dom_resu = dom_resu.intersection(arg._domain)
        # Restriction to the result domain
        self_r = self.restrict(dom_resu)
        args_r = [args[i].restrict(dom_resu) for i in range(p)]
        # Call of the FreeModuleTensor version
        return FreeModuleTensor.__call__(self_r, *args_r)

    def contract(self, *args):
        r"""
        Contraction with another tensor field, on one or more indices.

        INPUT:

        - ``pos1`` -- positions of the indices in ``self`` involved in the
          contraction; ``pos1`` must be a sequence of integers, with 0 standing
          for the first index position, 1 for the second one, etc. If ``pos1``
          is not provided, a single contraction on the last index position of
          ``self`` is assumed
        - ``other`` -- the tensor field to contract with
        - ``pos2`` -- positions of the indices in ``other`` involved in the
          contraction, with the same conventions as for ``pos1``. If ``pos2``
          is not provided, a single contraction on the first index position of
          ``other`` is assumed

        OUTPUT:

        - tensor field resulting from the contraction at the positions
          ``pos1`` and ``pos2`` of ``self`` with ``other``

        EXAMPLES:

        Contraction of a tensor field of type (2,0) with a tensor field of
        type (1,1)::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: a = M.tensor_field(2,0, name='a')
            sage: a[:] = [[1+x, 2], [y, -x^2]]
            sage: b = M.tensor_field(1,1, name='b')
            sage: b[:] = [[-y, 1], [x, x+y]]
            sage: s = a.contract(0, b, 1); s
            Tensor field of type (2,0) on the 2-dimensional differentiable manifold M
            sage: s.display()
            -x*y d/dx*d/dx + (x^2 + x*y + y^2 + x) d/dx*d/dy
             + (-x^2 - 2*y) d/dy*d/dx + (-x^3 - x^2*y + 2*x) d/dy*d/dy

        Check::

            sage: all([s[ind] == sum(a[k, ind[0]]*b[ind[1], k] for k in [0..1])
            ....:      for ind in M.index_generator(2)])
            True

        The same contraction with repeated index notation::

            sage: s == a['^ki']*b['^j_k']
            True

        Contraction on the second index of ``a``::

            sage: s = a.contract(1, b, 1); s
            Tensor field of type (2,0) on the 2-dimensional differentiable manifold M
            sage: s.display()
            (-(x + 1)*y + 2) d/dx*d/dx + (x^2 + 3*x + 2*y) d/dx*d/dy
             + (-x^2 - y^2) d/dy*d/dx + (-x^3 - (x^2 - x)*y) d/dy*d/dy

        Check::

            sage: all([s[ind] == sum(a[ind[0], k]*b[ind[1], k] for k in [0..1])
            ....:      for ind in M.index_generator(2)])
            True

        The same contraction with repeated index notation::

            sage: s == a['^ik']*b['^j_k']
            True

        See
        :meth:`sage.manifolds.differentiable.tensorfield.TensorField.contract`
        for more examples.

        """
        # This is to ensure the call to the TensorField version instead of
        # the FreeModuleTensor one
        return TensorField.contract(self, *args)

    def __mul__(self, other):
        r"""
        Tensor product (or multiplication of the right by a scalar).

        INPUT:

        - ``other`` -- a tensor field, on the same manifold as ``self`` (or an
          object that can be coerced to a scalar field on the same manifold
          as ``self``)

        OUPUT:

        - the tensor field resulting from the tensor product of ``self`` with
          ``other`` (or from the product ``other * self`` if ``other`` is a
          scalar)

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: a = M.tensor_field(0,2, name='a')
            sage: a[:] = [[1+x, 2], [y, -x^2]]

        Tensor product with another tensor field::

            sage: v = M.vector_field(name='v')
            sage: v[:] = [-y, x]
            sage: s = a.__mul__(v); s
            Tensor field a*v of type (1,2) on the 2-dimensional differentiable
             manifold M
            sage: s.display()
            a*v = -(x + 1)*y d/dx*dx*dx - 2*y d/dx*dx*dy - y^2 d/dx*dy*dx
             + x^2*y d/dx*dy*dy + (x^2 + x) d/dy*dx*dx + 2*x d/dy*dx*dy
             + x*y d/dy*dy*dx - x^3 d/dy*dy*dy
            sage: all([s[ind] == v[ind[0]] * a[ind[1],ind[2]]
            ....:      for ind in M.index_generator(3)])
            True

        Multiplication on the right by a scalar field::

            sage: f = M.scalar_field({X: x+y}, name='f')
            sage: s = a.__mul__(f); s
            Tensor field of type (0,2) on the 2-dimensional differentiable
             manifold M
            sage: s.display()
            (x^2 + (x + 1)*y + x) dx*dx + (2*x + 2*y) dx*dy + (x*y + y^2) dy*dx
             + (-x^3 - x^2*y) dy*dy
            sage: s == f*a
            True

        """
        # This is to ensure the call to the TensorField version instead of
        # the FreeModuleTensor one
        return TensorField.__mul__(self, other)


    def display_comp(self, frame=None, chart=None, coordinate_labels=True,
                     only_nonzero=True, only_nonredundant=False):
        r"""
        Display the tensor components w.r.t. a given frame, one per
        line.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``frame`` -- (default: ``None``) vector frame with respect to which
          the tensor field components are defined; if ``None``, then

          - if ``chart`` is not ``None``, the coordinate frame associated to
            ``chart`` is used
          - otherwise, the default basis of the vector field module on which
            the tensor field is defined is used

        - ``chart`` -- (default: ``None``) chart specifying the coordinate
          expression of the components; if ``None``, the default chart of the
          tensor field domain is used
        - ``coordinate_labels`` -- (default: ``True``) boolean; if ``True``,
          coordinate symbols are used by default (instead of integers) as
          index labels whenever ``frame`` is a coordinate frame
        - ``only_nonzero`` -- (default: ``True``) boolean; if ``True``, only
          nonzero components are displayed
        - ``only_nonredundant`` -- (default: ``False``) boolean; if ``True``,
          only nonredundant components are displayed in case of symmetries

        EXAMPLES:

        Display of the components of a type-(2,1) tensor field on a
        2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: t = M.tensor_field(2, 1, name='t', sym=(0,1))
            sage: t[0,0,0], t[0,1,0], t[1,1,1] = x+y, x*y, -3
            sage: t.display_comp()
            t^xx_x = x + y
            t^xy_x = x*y
            t^yx_x = x*y
            t^yy_y = -3

        By default, only the non-vanishing components are displayed; to see
        all the components, the argument ``only_nonzero`` must be set to
        ``False``::

            sage: t.display_comp(only_nonzero=False)
            t^xx_x = x + y
            t^xx_y = 0
            t^xy_x = x*y
            t^xy_y = 0
            t^yx_x = x*y
            t^yx_y = 0
            t^yy_x = 0
            t^yy_y = -3

        ``t`` being symmetric w.r.t. to its first two indices, one may ask to
        skip the components that can be deduced by symmetry::

            sage: t.display_comp(only_nonredundant=True)
            t^xx_x = x + y
            t^xy_x = x*y
            t^yy_y = -3

        Instead of coordinate labels, one may ask for integers::

            sage: t.display_comp(coordinate_labels=False)
            t^00_0 = x + y
            t^01_0 = x*y
            t^10_0 = x*y
            t^11_1 = -3

        Display in a frame different from the default one (note that since
        ``f`` is not a coordinate frame, integer are used to label the
        indices)::

            sage: a = M.automorphism_field()
            sage: a[:] = [[1+y^2, 0], [0, 2+x^2]]
            sage: f = X.frame().new_frame(a, 'f')
            sage: t.display_comp(frame=f)
            t^00_0 = (x + y)/(y^2 + 1)
            t^01_0 = x*y/(x^2 + 2)
            t^10_0 = x*y/(x^2 + 2)
            t^11_1 = -3/(x^2 + 2)

        Display with respect to a chart different from the default one::

            sage: Y.<u,v> = M.chart()
            sage: X_to_Y = X.transition_map(Y, [x+y, x-y])
            sage: Y_to_X = X_to_Y.inverse()
            sage: t.display_comp(chart=Y)
            t^uu_u = 1/4*u^2 - 1/4*v^2 + 1/2*u - 3/2
            t^uu_v = 1/4*u^2 - 1/4*v^2 + 1/2*u + 3/2
            t^uv_u = 1/2*u + 3/2
            t^uv_v = 1/2*u - 3/2
            t^vu_u = 1/2*u + 3/2
            t^vu_v = 1/2*u - 3/2
            t^vv_u = -1/4*u^2 + 1/4*v^2 + 1/2*u - 3/2
            t^vv_v = -1/4*u^2 + 1/4*v^2 + 1/2*u + 3/2

        Note that the frame defining the components is the coordinate frame
        associated with chart ``Y``, i.e. we have::

            sage: str(t.display_comp(chart=Y)) == str(t.display_comp(frame=Y.frame(), chart=Y))
            True

        Display of the components with respect to a specific frame, expressed
        in terms of a specific chart::

            sage: t.display_comp(frame=f, chart=Y)
            t^00_0 = 4*u/(u^2 - 2*u*v + v^2 + 4)
            t^01_0 = (u^2 - v^2)/(u^2 + 2*u*v + v^2 + 8)
            t^10_0 = (u^2 - v^2)/(u^2 + 2*u*v + v^2 + 8)
            t^11_1 = -12/(u^2 + 2*u*v + v^2 + 8)

        """
        from sage.misc.latex import latex
        from sage.manifolds.differentiable.vectorframe import CoordFrame
        if frame is None:
            if chart is not None:
                frame = chart.frame()
            else:
                frame = self._fmodule.default_basis()
        if chart is None:
            chart = self._domain.default_chart()
        index_labels = None
        index_latex_labels = None
        if isinstance(frame, CoordFrame) and coordinate_labels:
            ch = frame.chart()
            index_labels = map(str, ch[:])
            index_latex_labels = map(latex, ch[:])
        return FreeModuleTensor.display_comp(self, basis=frame,
                                  format_spec=chart, index_labels=index_labels,
                                  index_latex_labels=index_latex_labels,
                                  only_nonzero=only_nonzero,
                                  only_nonredundant=only_nonredundant)

    def at(self, point):
        r"""
        Value of the tensor field at a point of its domain.

        If the current tensor field is

        .. MATH::

            t:\ U  \longrightarrow T^{(k,l)} M

        associated with the differentiable map

        .. MATH::

            \Phi:\ U \longrightarrow M

        where `U` and `M` are two manifolds (possibly `U=M` and
        `\Phi=\mathrm{Id}_M`), then for any point `p\in U`, `t(p)` is a tensor
        on the tangent space to `M` at the point `\Phi(p)`.

        INPUT:

        - ``point`` -- (instance of
          :class:`~sage.manifolds.point.TopologicalManifoldPoint`) point `p` in the
          domain of the tensor field (`U`)

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`
          representing the tensor `t(p)` on the tangent vector space
          `T_{\Phi(p)} M`

        EXAMPLES:

        Vector in a tangent space of a 2-dimensional manifold::

            sage: from sage.manifolds.differentiable.tangent_space import TangentSpace  # for doctests only
            sage: TangentSpace._clear_cache_()  # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: p = M.point((-2,3), name='p')
            sage: v = M.vector_field('v')
            sage: v[:] = [y, x^2] ; v.display()
            v = y d/dx + x^2 d/dy
            sage: vp = v.at(p) ; vp
            Tangent vector v at Point p on the 2-dimensional differentiable
             manifold M
            sage: vp.parent()
            Tangent space at Point p on the 2-dimensional differentiable
             manifold M
            sage: vp.display()
            v = 3 d/dx + 4 d/dy

        A 1-form gives birth to a linear form in the tangent space::

            sage: w = M.one_form('w')
            sage: w[:] = [-x, 1+y] ; w.display()
            w = -x dx + (y + 1) dy
            sage: wp = w.at(p) ; wp
            Linear form w on the Tangent space at Point p on the 2-dimensional
             differentiable manifold M
            sage: wp.parent()
            Dual of the Tangent space at Point p on the 2-dimensional
             differentiable manifold M
            sage: wp.display()
            w = 2 dx + 4 dy

        A tensor field of type (1,1) yields a tensor of type(1,1) in the
        tangent space::

            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[0,0], t[0,1], t[1,1] = 1+x, x*y, 1-y
            sage: t.display()
            t = (x + 1) d/dx*dx + x*y d/dx*dy + (-y + 1) d/dy*dy
            sage: tp = t.at(p) ; tp
            Type-(1,1) tensor t on the Tangent space at Point p on the
             2-dimensional differentiable manifold M
            sage: tp.parent()
            Free module of type-(1,1) tensors on the Tangent space at Point p
             on the 2-dimensional differentiable manifold M
            sage: tp.display()
            t = -d/dx*dx - 6 d/dx*dy - 2 d/dy*dy

        A 2-form yields an alternating form of degree 2 in the tangent space::

            sage: a = M.diff_form(2, name='a')
            sage: a[0,1] = x*y
            sage: a.display()
            a = x*y dx/\dy
            sage: ap = a.at(p) ; ap
            Alternating form a of degree 2 on the Tangent space at Point p on
             the 2-dimensional differentiable manifold M
            sage: ap.parent()
            2nd exterior power of the dual of the Tangent space at Point p on
             the 2-dimensional differentiable manifold M
            sage: ap.display()
            a = -6 dx/\dy

        Example with a non trivial map `\Phi`::

            sage: U = Manifold(1, 'U')  # (0,2*pi) as a 1-dimensional manifold
            sage: T.<t> = U.chart(r't:(0,2*pi)')  # canonical chart on U
            sage: Phi = U.diff_map(M, [cos(t), sin(t)], name='Phi',
            ....:                  latex_name=r'\Phi')
            sage: v = U.vector_field(name='v', dest_map=Phi) ; v
            Vector field v along the 1-dimensional differentiable manifold U
             with values on the 2-dimensional differentiable manifold M
            sage: v[:] = [1+t, t^2]
            sage: v.display()
            v = (t + 1) d/dx + t^2 d/dy
            sage: p = U((pi/6,))
            sage: vp = v.at(p) ; vp
            Tangent vector v at Point on the 2-dimensional differentiable
             manifold M
            sage: vp.parent() is M.tangent_space(Phi(p))
            True
            sage: vp.display()
            v = (1/6*pi + 1) d/dx + 1/36*pi^2 d/dy

        """
        if point not in self._domain:
            raise ValueError("the {} is not in the domain of ".format(point) +
                             "the {}".format(self))
        dest_map = self._fmodule._dest_map
        if dest_map.is_identity():
            amb_point = point
        else:
            amb_point = dest_map(point)  #  "ambient" point
        ts = amb_point._manifold.tangent_space(amb_point)
        resu = ts.tensor(self._tensor_type, name=self._name,
                         latex_name=self._latex_name, sym=self._sym,
                         antisym=self._antisym)
        for frame, comp in self._components.iteritems():
            comp_resu = resu.add_comp(frame.at(point))
            for ind, val in comp._comp.iteritems():
                comp_resu._comp[ind] = val(point)
        return resu