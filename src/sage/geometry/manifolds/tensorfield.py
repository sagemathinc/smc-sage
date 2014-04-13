r"""
Tensor fields

The class :class:`TensorField` implements tensor fields on differentiable 
manifolds over `\RR`. 

A tensor field of type `(k,\ell)` is a field of multilinear maps:

.. MATH::

    \underbrace{T_p^*M\times\cdots\times T_p^*M}_{k\ \; \mbox{times}}
    \times \underbrace{T_p M\times\cdots\times T_p M}_{\ell\ \; \mbox{times}}
    \longrightarrow \RR
    
where `T_p M` stands for the tangent space at the point `p` on the
manifold `M` and `T_p^*M` for its dual vector space. The integer `k+\ell`
is called the tensor rank. 

Various derived classes of :class:`TensorField` are devoted to specific tensor
fields:

* :class:`~sage.geometry.manifolds.vectorfield.VectorField` for vector fields 
  (rank-1 contravariant tensor fields)
* :class:`~sage.geometry.manifolds.diffform.OneForm` for 1-forms (rank-1 
  covariant tensor fields)
* :class:`~sage.geometry.manifolds.rank2field.EndomorphismField` for fields of 
  endomorphisms (type (1,1) tensor fields)
* :class:`~sage.geometry.manifolds.rank2field.SymBilinFormField` for fields of 
  symmetric bilinear forms (rank-2 symmetric covariant tensor fields)
* :class:`~sage.geometry.manifolds.diffform.DiffForm` for differential forms 
  (fully antisymmetric covariant tensor fields)


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013, 2014) : initial version

EXAMPLES:

    A tensor field of type (1,1) on a 2-dimensional manifold::

        sage: M = Manifold(2, 'M', start_index=1)
        sage: c_xy.<x,y> = M.chart('x y')
        sage: t = M.tensor_field(1, 1, 'T') ; t
        tensor field 'T' of type (1,1) on the 2-dimensional manifold 'M'
        sage: t.tensor_rank
        2

    A just-created tensor field has no components::

        sage: t.components
        {}

    Components w.r.t. the manifold's default frame are created by providing the
    relevant indices inside square brackets::

        sage: t[1,1] = x^2

    Unset components are initialized to zero::

        sage: t[:]  # list of components w.r.t. the manifold's default vector frame
        [x^2   0]
        [  0   0]

    The full set of components w.r.t. a given vector frame is returned by the 
    method :meth:`comp`; it is an instance of the class 
    :class:`~sage.tensor.modules.comp.Components`::
    
        sage: t.comp(c_xy.frame)
        2-indices components w.r.t. coordinate frame (M, (d/dx,d/dy)) 
        sage: print type(t.comp(c_xy.frame))
        <class 'sage.tensor.modules.comp.Components'>

    The vector frame can be skipped, it is then assumed to be the
    manifold's default frame::
    
        sage: M.default_frame()
        coordinate frame (M, (d/dx,d/dy))
        sage: t.comp() is t.comp(c_xy.frame)
        True

    Individual components w.r.t. the manifold's default frame are accessed by 
    listing their indices inside double square brackets; they are scalar
    fields on the manifold, and therefore instances of the class 
    :class:`~sage.geometry.manifolds.scalarfield.ScalarField`::

        sage: t[[1,1]]
        scalar field on the 2-dimensional manifold 'M'
        sage: t[[1,1]].expr()
        x^2
        sage: t[[1,2]]
        zero scalar field on the 2-dimensional manifold 'M'
        sage: t[[1,2]].expr()
        0
        
    A direct access to the coordinate expression of some component is obtained
    via the single square brackets::
    
        sage: t[1,1] 
        x^2
        sage: t[1,1] is t[[1,1]].function_chart() # the coordinate function
        True
        sage: t[1,1] is t[[1,1]].function_chart(c_xy)
        True
        sage: t[1,1] == t[[1,1]].expr() # check the value of the coordinate function
        True
        sage: t[1,1].expr() is t[[1,1]].expr() # the symbolic expression
        True

    In other words, the single square brackets return an instance of 
    :class:`~sage.geometry.manifolds.chart.FunctionChart` that is the 
    coordinate function representing the component in some chart (by default, 
    the manifold's default chart)::
    
        sage: print type(t[1,1])    # single bracket --> FunctionChart
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        sage: print type(t[[1,1]])  # double bracket --> ScalarField
        <class 'sage.geometry.manifolds.scalarfield.ScalarFieldRing_with_category.element_class'>
    
    Expressions in a chart different from the manifold's default one are 
    obtained by specifying the chart as the last argument inside the
    single square brackets::
    
        sage: c_uv.<u,v> = M.chart('u v')
        sage: xy_to_uv = c_xy.coord_change(c_uv, x+y, x-y)  
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: t[1,1, c_uv] 
        1/4*u^2 + 1/2*u*v + 1/4*v^2

    Note that ``t[1,1, c_uv]`` is the component of the tensor t w.r.t. to 
    the coordinate frame associated to the chart (x,y) expressed in terms of 
    the coordinates (u,v). Indeed, ``t[1,1, c_uv]`` is a shortcut for 
    ``t.comp(c_xy.frame)[[1,1]].function_chart(c_uv)``::
        
        sage: t[1,1, c_uv] is t.comp(c_xy.frame)[[1,1]].function_chart(c_uv)
        True
    
    Similarly, ``t[1,1]`` is a shortcut for 
    ``t.comp(c_xy.frame)[[1,1]].function_chart(c_xy)``::
    
        sage: t[1,1] is t.comp(c_xy.frame)[[1,1]].function_chart(c_xy)            
        True
        sage: t[1,1] is t.comp()[[1,1]].function_chart()  # since c_xy.frame and c_xy are the manifold's default values                
        True

    Internally, the components are stored as a dictionary (attribute 
    :attr:`_comp` of the class 
    :class:`~sage.tensor.modules.comp.Components`) whose
    keys are the indices. Only the non-zero components and non-redundant
    components (in case of symmetries) are stored::

        sage: t.comp()._comp
        {(1, 1): scalar field on the 2-dimensional manifold 'M'}

    All the components can be set at once via [:]::
    
        sage: t[:] = [[1, -x], [x*y, 2]]
        sage: t[:]
        [  1  -x]
        [x*y   2]
        
    The different sets of components, corresponding to representations of the
    tensor in different vector frames, are stored in the dictionary 
    :attr:`components`, each item being an instance of the class 
    :class:`~sage.tensor.modules.comp.Components`::
    
        sage: t.components
        {coordinate frame (M, (d/dx,d/dy)): 2-indices components w.r.t. coordinate frame (M, (d/dx,d/dy))}
        sage: print type(t.components[c_xy.frame])
        <class 'sage.tensor.modules.comp.Components'>
        sage: print type(t.comp())
        <class 'sage.tensor.modules.comp.Components'>
        sage: t.comp() is t.components[c_xy.frame]
        True

    To set the components in a vector frame different from the manifold's 
    default one, the method :meth:`set_comp` must be employed::

        sage: e = M.vector_frame('e')
        sage: t.set_comp(e)[1,1], t.set_comp(e)[1,2] = (x+y, 0)
        sage: t.set_comp(e)[2,1], t.set_comp(e)[2,2] = (y, -3*x)
        sage: t.comp(e)
        2-indices components w.r.t. vector frame (M, (e_1,e_2))
        sage: t.comp(e)[:]
        [x + y     0]
        [    y  -3*x]

    All the components in some frame can be set at once, via the operator
    [:]::

        sage: t.set_comp(e)[:] = [[x+y, 0], [y, -3*x]]
        sage: t.comp(e)[:]  # same as above:
        [x + y     0]
        [    y  -3*x]
    
    To avoid any insconstency between the various components, the method 
    :meth:`set_comp` clears the components in other frames. 
    Accordingly, the components in the frame c_xy.frame have been deleted::
    
        sage: t.components
        {vector frame (M, (e_1,e_2)): 2-indices components w.r.t. vector frame (M, (e_1,e_2))}

    To keep the other components, one must use the method :meth:`add_comp`::
    
        sage: t = M.tensor_field(1, 1, 'T')  # Let us restart 
        sage: t[:] = [[1, -x], [x*y, 2]]  # by first setting the components in the frame c_xy.frame
        sage: # We now set the components in the frame e with add_comp:
        sage: t.add_comp(e)[:] = [[x+y, 0], [y, -3*x]]
        sage: t.components  # Both set of components are present:
        {coordinate frame (M, (d/dx,d/dy)): 2-indices components w.r.t. coordinate frame (M, (d/dx,d/dy)), vector frame (M, (e_1,e_2)): 2-indices components w.r.t. vector frame (M, (e_1,e_2))}

    The expansion of the tensor field in a given frame is displayed via the 
    method :meth:`view` (the symbol * stands for tensor product)::
    
        sage: t.view()  # expansion in the manifold's default frame
        T = d/dx*dx - x d/dx*dy + x*y d/dy*dx + 2 d/dy*dy
        sage: t.view(e)
        T = (x + y) e_1*e^1 + y e_2*e^1 - 3*x e_2*e^2

    A tensor field acts as a multilinear map on 1-forms and vector fields; 
    in the present case, T being of type (1,1), it acts on pairs 
    (1-form, vector)::
    
        sage: a = M.one_form('a')
        sage: a[:] = (1, x)
        sage: v = M.vector_field('V')
        sage: v[:] = (y, 2)
        sage: t(a,v)
        scalar field 'T(a,V)' on the 2-dimensional manifold 'M'
        sage: t(a,v).expr()
        x^2*y^2 + 2*x + y
        sage: latex(t(a,v))
        T\left(a,V\right)
    
    Check by means of the component expression of t(a,v)::
    
        sage: t[1,1]*a[1]*v[1] + t[1,2]*a[1]*v[2] + t[2,1]*a[2]*v[1] + t[2,2]*a[2]*v[2] - t(a,v).expr()
        0

    A scalar field (rank-0 tensor field)::
    
        sage: f = M.scalar_field(x*y + 2, name='f') ; f 
        scalar field 'f' on the 2-dimensional manifold 'M'
        sage: f.tensor_type
        (0, 0)
        
    A scalar field acts on points on the manifold::
    
        sage: p = M.point((1,2))
        sage: f(p)
        4
        
    A vector field (rank-1 contravariant tensor field)::
    
        sage: v = M.vector_field('v') ; v
        vector field 'v' on the 2-dimensional manifold 'M'
        sage: v.tensor_type
        (1, 0)
        sage: v[1], v[2] = -x, y
        sage: v.view()
        v = -x d/dx + y d/dy        

    A field of symmetric bilinear forms::
    
        sage: q = M.sym_bilin_form_field('Q') ; q
        field of symmetric bilinear forms 'Q' on the 2-dimensional manifold 'M'
        sage: q.tensor_type
        (0, 2)

    The components of a symmetric bilinear form are dealt by the subclass 
    :class:`~sage.tensor.modules.comp.CompFullySym` of the class 
    :class:`~sage.tensor.modules.comp.Components`, which takes into 
    account the symmetry between the two indices::
    
        sage: q[1,1], q[1,2], q[2,2] = (0, -x, y) # no need to set the component (2,1)
        sage: print type(q.comp())
        <class 'sage.tensor.modules.comp.CompFullySym'>
        sage: q[:] # note that the component (2,1) is equal to the component (1,2)
        [ 0 -x]
        [-x  y]
        sage: q.view()
        Q = -x dx*dy - x dy*dx + y dy*dy
    
    Internally (dictionary :attr:`_comp` of the class 
    :class:`~sage.tensor.modules.comp.Components`), only
    the non-zero and non-redundant components are stored::
    
        sage: q.comp()._comp
        {(1, 2): scalar field on the 2-dimensional manifold 'M',
        (2, 2): scalar field on the 2-dimensional manifold 'M'}
        sage: q.comp()._comp[(1,2)].expr()
        -x
        sage: q.comp()._comp[(2,2)].expr()
        y

    More generally, tensor symmetries or antisymmetries can be specified via
    the keywords ``sym`` and ``antisym``. For instance a rank-4 covariant 
    tensor symmetric with respect to its first two arguments and 
    antisymmetric with respect to its last two ones is declared as follows::
    
        sage: t = M.tensor_field(0, 4, 'T', sym=(0,1), antisym=(2,3))
        sage: t[1,2,1,2] = 3
        sage: t[2,1,1,2] # check of the symmetry with respect to the first 2 indices
        3
        sage: t[1,2,2,1] # check of the antisymmetry with respect to the last 2 indices
        -3

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

from sage.tensor.modules.free_module_tensor import FreeModuleTensor

class TensorFieldParal(FreeModuleTensor):
    r"""
    Base class for tensor fields on an open set of a differentiable manifold, 
    with values on parallelizable open subset of a differentiable manifold. 
    
    An instance of this class is a tensor field along an open subset `U` 
    of some manifold `S` with values in a parallelizable open subset `V` 
    of a manifold `M`, via a differentiable mapping `\Phi: U \rightarrow V`. 
    The standard case of a tensor field *on* a manifold corresponds to `S=M`, 
    `U=V` and `\Phi = \mathrm{Id}`. Another common case is `\Phi` being an
    immersion.

    A tensor field of type `(k,\ell)` is a field `t` on `U`, so that at each 
    point `p\in U`, `t(p)` is a multilinear map of the type:

    .. MATH::

        t(p):\ \underbrace{T_p^*M\times\cdots\times T_p^*M}_{k\ \; \mbox{times}}
        \times \underbrace{T_p M\times\cdots\times T_p M}_{\ell\ \; \mbox{times}}
        \longrightarrow \RR
    
    where `T_p M` stands for the tangent space at the point `p` on the
    manifold `M` and `T_p^* M` for its dual vector space. The integer `k+\ell`
    is called the tensor rank. 
    
    INPUT:
    
    - ``vector_field_module`` -- free module `\mathcal{X}(U,\Phi)` of vector 
      fields along `U` with values on `\Phi(U)\subset V \subset M`
    - ``tensor_type`` -- pair (k,l) with k being the contravariant rank and l 
      the covariant rank
    - ``name`` -- (default: None) name given to the tensor field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the tensor field; 
      if none is provided, the LaTeX symbol is set to ``name``
    - ``sym`` -- (default: None) a symmetry or a list of symmetries among the 
      tensor arguments: each symmetry is described by a tuple containing 
      the positions of the involved arguments, with the convention position=0
      for the first argument. For instance:
        * sym=(0,1) for a symmetry between the 1st and 2nd arguments 
        * sym=[(0,2),(1,3,4)] for a symmetry between the 1st and 3rd
          arguments and a symmetry between the 2nd, 4th and 5th arguments.
    - ``antisym`` -- (default: None) antisymmetry or list of antisymmetries 
      among the arguments, with the same convention as for ``sym``. 


    EXAMPLES:

    A tensor field of type (2,0) on a 3-dimensional manifold::
    
        sage: M = Manifold(3, 'M')
        sage: c_xyz.<x,y,z> = M.chart('x y z')
        sage: t = M.tensor_field(2, 0, 'T') ; t
        tensor field 'T' of type (2,0) on the 3-dimensional manifold 'M'

    Tensor fields are considered as elements of a module over the ring 
    `C^\infty(M)` of scalar fields on `M`::
    
        sage: t.parent()
        free module TF^(2,0)(M) of type-(2,0) tensors fields on the 3-dimensional manifold 'M'
        sage: t.parent().base_ring()
        ring of scalar fields on the 3-dimensional manifold 'M'

    The components with respect to the manifold's default frame are set or read
    by means of square brackets::
    
        sage: e = M.vector_frame('e') ; M.set_default_frame(e)
        sage: for i in M.irange():
        ...       for j in M.irange():
        ...           t[i,j] = (i+1)**(j+1)
        ...
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
        2-indices components w.r.t. vector frame (M, (f_0,f_1,f_2))
        sage: t.comp(f)[0,0]
        -3
        sage: t.comp(f)[:]  # the full list of components
        [-3  0  0]
        [ 0  0  0]
        [ 0  0  0]

    To avoid any insconstency between the various components, the method 
    :meth:`set_comp` deletes the components in other frames. 
    Accordingly, the components in the frame e have been deleted::
    
        sage: t.components
        {vector frame (M, (f_0,f_1,f_2)): 2-indices components w.r.t. vector frame (M, (f_0,f_1,f_2))}

    To keep the other components, one must use the method :meth:`add_comp`::
    
        sage: t = M.tensor_field(2, 0, 'T')  # Let us restart
        sage: t[0,0] = 2                     # sets the components in the frame e
        sage: # We now set the components in the frame f with add_comp:
        sage: t.add_comp(f)[0,0] = -3
        sage: # The components w.r.t. frame e have been kept: 
        sage: t.components
        {vector frame (M, (e_0,e_1,e_2)): 2-indices components w.r.t. vector frame (M, (e_0,e_1,e_2)), vector frame (M, (f_0,f_1,f_2)): 2-indices components w.r.t. vector frame (M, (f_0,f_1,f_2))}

    The basic attributes of :class:`TensorField` are :attr:`domain`, 
    :attr:`tensor_type` (the pair (k,l)), :attr:`tensor_rank` (the integer k+l),  
    and :attr:`components` (the dictionary of the components w.r.t. various 
    frames)::

        sage: t.domain
        3-dimensional manifold 'M'
        sage: t.tensor_type
        (2, 0)
        sage: t.tensor_rank
        2
        sage: t.components
        {vector frame (M, (e_0,e_1,e_2)): 2-indices components w.r.t. vector frame (M, (e_0,e_1,e_2)), vector frame (M, (f_0,f_1,f_2)): 2-indices components w.r.t. vector frame (M, (f_0,f_1,f_2))}

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
        tensor field of type (2,0) on the 3-dimensional manifold 'M'
        sage: a[:], (2*b)[:], s[:]
        (
        [1 2 3]  [ 8  0 14]  [ 9  2 17]
        [0 0 0]  [ 0 10  0]  [ 0 10  0]
        [0 0 0], [ 0  0 12], [ 0  0 12]
        )
        sage: s = a - b ; s
        tensor field of type (2,0) on the 3-dimensional manifold 'M'
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
        tensor field of type (4,0) on the 3-dimensional manifold 'M'
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
        tensor field of type (2,2) on the 3-dimensional manifold 'M'
        sage: s.symmetries()
        symmetry: (0, 1);  antisymmetry: (2, 3)
        sage: s1 = d*a ; s1 
        tensor field of type (2,2) on the 3-dimensional manifold 'M'
        sage: s1.symmetries()
        symmetry: (0, 1);  antisymmetry: (2, 3)
        sage: d*a == a*d
        True

    """
    def __init__(self, vector_field_module, tensor_type, name=None, 
                 latex_name=None, sym=None, antisym=None):
        FreeModuleTensor.__init__(self, vector_field_module, tensor_type,
                                  name=name, latex_name=latex_name,
                                  sym=sym, antisym=antisym)
        self.domain = vector_field_module.domain
        self.ambient_domain = vector_field_module.ambient_domain
        # Initialization of derived quantities:
        self._init_derived() 

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "tensor field"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " of type (%s,%s) " % (str(self.tensor_type[0]), 
                                             str(self.tensor_type[1]))
        return self._final_repr(description)

    def _final_repr(self, description):
        r"""
        Part of string representation common to all tensor fields.
        """
        if self.domain == self.ambient_domain:
            description += "on the " + str(self.domain)
        else:
            description += "along the " + str(self.domain) + \
                           " with values on the " + str(self.ambient_domain)
        return description
        
    def _new_instance(self):
        r"""
        Create a :class:`TensorFieldParal` instance of the same tensor type, 
        with the same symmetries and on the same domain.

        This method must be redefined by derived classes of 
        :class:`TensorFieldParal`.
        
        """
        return TensorFieldParal(self.fmodule, self.tensor_type, sym=self.sym, 
                                antisym=self.antisym)

    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        self._lie_derivatives = {} # collection of Lie derivatives of self

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        # First deletes any reference to self in the vectors' dictionary:
        for vid, val in self._lie_derivatives.items():
            del val[0]._lie_der_along_self[id(self)]
        self._lie_derivatives.clear()


    def common_coord_frame(self, other):
        r"""
        Find a common coordinate frame for the components of ``self`` and 
        ``other``. 
        
        In case of multiple common bases, the domain's default coordinate
        basis is privileged. 
        If the current components of ``self`` and ``other`` are all relative to
        different frames, a common frame is searched by performing a component
        transformation, via the transformations listed in 
        ``self.domain.frame_changes``, still privileging transformations to 
        the domain's default frame.
        
        INPUT:
        
        - ``other`` -- a tensor field
        
        OUPUT:
        
        - common coordinate frame; if no common basis is found, None is 
          returned. 
        
        """
        from vectorframe import CoordFrame
        # Compatibility checks:
        if not isinstance(other, TensorFieldParal):
            raise TypeError("The argument must be of type TensorFieldParal.")
        dom = self.domain
        def_frame = dom.def_frame
        #
        # 1/ Search for a common frame among the existing components, i.e. 
        #    without performing any component transformation. 
        #    -------------------------------------------------------------
        # 1a/ Direct search
        if def_frame in self.components and \
           def_frame in other.components and \
           isinstance(dom.def_frame, CoordFrame):
            return def_frame # the domain's default frame is privileged
        for frame1 in self.components:
            if frame1 in other.components and \
               isinstance(frame1, CoordFrame):
                return frame1
        # 1b/ Search involving subframes
        dom2 = other.domain
        for frame1 in self.components:
            if not isinstance(frame1, CoordFrame):
                continue
            for frame2 in other.components:
                if not isinstance(frame2, CoordFrame):
                    continue
                if frame2 in frame1.subframes:
                    self.comp(frame2)
                    return frame2
                if frame1 in frame2.subframes:
                    other.comp(frame1)
                    return frame1
        #
        # 2/ Search for a common frame via one component transformation
        #    ----------------------------------------------------------
        # If this point is reached, it is indeed necessary to perform at least 
        # one component transformation to get a common frame
        if isinstance(dom.def_frame, CoordFrame):
            if def_frame in self.components:
                for oframe in other.components:
                    if (oframe, def_frame) in dom.frame_changes:
                        other.comp(def_frame, from_frame=oframe)
                        return def_frame
            if def_frame in other.components:
                for sframe in self.components:
                    if (sframe, def_frame) in dom.frame_changes:
                        self.comp(def_frame, from_frame=sframe)
                        return def_frame
        # If this point is reached, then def_frame cannot be a common frame
        # via a single component transformation
        for sframe in self.components:
            if not instance(sframe, CoordFrame):
                continue
            for oframe in other.components:
                if not instance(oframe, CoordFrame):
                    continue
                if (oframe, sframe) in dom.frame_changes:
                    other.comp(sframe, from_frame=oframe)
                    return sframe
                if (sframe, oframe) in dom.frame_changes:
                    self.comp(oframe, from_frame=sframe)
                    return oframe
        #
        # 3/ Search for a common frame via two component transformations
        #    -----------------------------------------------------------
        # If this point is reached, it is indeed necessary to perform at two
        # component transformation to get a common frame
        for sframe in self.components:
            for oframe in other.components:
                if (sframe, def_frame) in dom.frame_changes and \
                   (oframe, def_frame) in dom.frame_changes and \
                   isinstance(def_frame, CoordFrame):
                    self.comp(def_frame, from_frame=sframe)
                    other.comp(def_frame, from_frame=oframe)
                    return def_frame
                for frame in dom.frames:
                    if (sframe, frame) in dom.frame_changes and \
                       (oframe, frame) in dom.frame_changes and \
                       isinstance(frame, CoordFrame):
                        self.comp(frame, from_frame=sframe)
                        other.comp(frame, from_frame=oframe)
                        return frame
        #
        # If this point is reached, no common frame could be found, even at 
        # the price of component transformations:
        return None
    

    def lie_der(self, vector):
        r"""
        Computes the Lie derivative with respect to a vector field.
        
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
            sage: c_xy.<x,y> = M.chart('x y')
            sage: v = M.vector_field('v')
            sage: v[:] = (-y, x)
            sage: w = M.vector_field()
            sage: w[:] = (2*x+y, x*y)
            sage: w.lie_der(v)
            vector field on the 2-dimensional manifold 'M'
            sage: w.lie_der(v).view()
            ((x - 2)*y + x) d/dx + (x^2 - y^2 - 2*x - y) d/dy

        The Lie derivative is antisymmetric::
        
            sage: w.lie_der(v) == -v.lie_der(w)
            True
            
        For vectors, it coincides with the commutator::

            sage: f = M.scalar_field(x^3 + x*y^2)
            sage: w.lie_der(v)(f).view()
            (x, y) |--> -(x + 2)*y^3 + 3*x^3 - x*y^2 + 5*(x^3 - 2*x^2)*y
            sage: w.lie_der(v)(f) == v(w(f)) - w(v(f))  # rhs = commutator [v,w] acting on f
            True
            
        Lie derivative of a 1-form::
        
            sage: om = M.one_form()
            sage: om[:] = (y^2*sin(x), x^3*cos(y))
            sage: om.lie_der(v)
            1-form on the 2-dimensional manifold 'M'
            sage: om.lie_der(v).view()
            (-y^3*cos(x) + x^3*cos(y) + 2*x*y*sin(x)) dx + (-x^4*sin(y) - 3*x^2*y*cos(y) - y^2*sin(x)) dy
            
        Check of Cartan identity::
        
            sage: om.lie_der(v) == v.contract(0,om.exterior_der(),0) + (om(v)).exterior_der()
            True
        
        """
        from vectorfield import VectorFieldParal
        if not isinstance(vector, VectorFieldParal):
            raise TypeError("The argument must be of type VectorFieldParal.")
        if id(vector) not in self._lie_derivatives:
            # A new computation must be performed
            #
            # 1/ Search for a common coordinate frame:
            coord_frame = self.common_coord_frame(vector)
            if coord_frame is None:
                raise TypeError("No common coordinate frame found.")
            chart = coord_frame.chart
            #
            # 2/ Component computation:
            tc = self.components[coord_frame]
            vc = vector.components[coord_frame]
            # the result has the same tensor type and same symmetries as self:
            resc = self._new_comp(coord_frame) 
            n_con = self.tensor_type[0]
            vf_module = vector.fmodule
            for ind in resc.non_redundant_index_generator():
                rsum = 0
                for i in vf_module.irange():
                    rsum += vc[[i]].function_chart(chart) * \
                           tc[[ind]].function_chart(chart).diff(i)
                # loop on contravariant indices:
                for k in range(n_con): 
                    for i in vf_module.irange():
                        indk = list(ind)
                        indk[k] = i  
                        rsum -= tc[[indk]].function_chart(chart) * \
                                vc[[ind[k]]].function_chart(chart).diff(i)
                # loop on covariant indices:
                for k in range(n_con, self.tensor_rank): 
                    for i in vf_module.irange():
                        indk = list(ind)
                        indk[k] = i  
                        rsum += tc[[indk]].function_chart(chart) * \
                                vc[[i]].function_chart(chart).diff(ind[k])
                resc[[ind]] = rsum.scalar_field()
            #
            # 3/ Final result (the tensor)
            res = vf_module.tensor_from_comp(self.tensor_type, resc)
            self._lie_derivatives[id(vector)] = (vector, res)
            vector._lie_der_along_self[id(self)] = self
        return self._lie_derivatives[id(vector)][1]
